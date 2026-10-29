#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TRINOVA.OS - Biblioteca genética exploratoria, procesador VCF/WGS y ALU ternaria.
Windows / Python 3.x. Requiere únicamente la biblioteca estándar.

Nota: las métricas de coherencia y entropía son métricas exploratorias de la
codificación ternaria de variantes. No son resultados clínicos ni diagnósticos.
"""

import ctypes
from ctypes import wintypes
import gzip
import math
import os
import time
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

WORD_SIZE = 18
TRIT_NEG, TRIT_ZERO, TRIT_POS = -1, 0, 1
TRIT_SET = (-1, 0, 1)


VCF_REAL_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\WGS_PROBAND01.vcf.gz"
)

ALU_TNOT, ALU_TAND, ALU_TOR, ALU_TXOR = 0, 1, 2, 3
ALU_TNAND, ALU_TNOR, ALU_TXNOR, ALU_TMAJ = 4, 5, 6, 7
ALU_TSTEP, ALU_TADD, ALU_TSUB, ALU_TMUL = 8, 9, 10, 11
ALU_TCMP, ALU_TSHL, ALU_TSHR = 12, 13, 14

ALU_NAMES = {
    ALU_TNOT: "TNOT", ALU_TAND: "TAND", ALU_TOR: "TOR",
    ALU_TXOR: "TXOR", ALU_TNAND: "TNAND", ALU_TNOR: "TNOR",
    ALU_TXNOR: "TXNOR", ALU_TMAJ: "TMAJ", ALU_TSTEP: "TSTEP",
    ALU_TADD: "TADD", ALU_TSUB: "TSUB", ALU_TMUL: "TMUL",
    ALU_TCMP: "TCMP", ALU_TSHL: "TSHL", ALU_TSHR: "TSHR",
}


# =============================================================================
# ALU TRIT-LEVEL
# =============================================================================

def trit_valido(x):
    return x in TRIT_SET


def validar_trit(x):
    return x if trit_valido(x) else 0


def tnot(a):
    return -a if trit_valido(a) else 0


def tand(a, b):
    if not trit_valido(a) or not trit_valido(b):
        return 0
    if a == 0 or b == 0:
        return 0
    return a if a == b else -1


def tor(a, b):
    if not trit_valido(a) or not trit_valido(b):
        return 0
    if a == 1 or b == 1:
        return 1
    if a == b:
        return a
    return -1


def tnand(a, b):
    return tnot(tand(a, b))


def tnor(a, b):
    return tnot(tor(a, b))


def txor(a, b):
    if not trit_valido(a) or not trit_valido(b):
        return 0
    if a == b:
        return 0
    if a == 0:
        return b
    if b == 0:
        return a
    return -a


def txnor(a, b):
    return tnot(txor(a, b))


def tmaj(a, b, c):
    if not all(trit_valido(x) for x in (a, b, c)):
        return 0
    total = a + b + c
    return 1 if total > 0 else -1 if total < 0 else 0


def tstep(x):
    return 1 if x > 0 else -1 if x < 0 else 0


def tadd(a, b, cin=0):
    if not all(trit_valido(x) for x in (a, b, cin)):
        return 0, 0
    total = a + b + cin
    if total > 1:
        return total - 3, 1
    if total < -1:
        return total + 3, -1
    return total, 0


def tsub(a, b, cin=0):
    if not all(trit_valido(x) for x in (a, b, cin)):
        return 0, 0
    return tadd(a, -b, cin)


def tmul(a, b):
    return a * b if trit_valido(a) and trit_valido(b) else 0


def tshl(a):
    return a if trit_valido(a) else 0


def tshr(a):
    return 0


def tcmp(a, b):
    if not trit_valido(a) or not trit_valido(b):
        return 0
    return -1 if a < b else 1 if a > b else 0


def alu_execute(op, a, b=0, cin=0):
    operations = {
        ALU_TNOT: lambda: (tnot(a), 0),
        ALU_TAND: lambda: (tand(a, b), 0),
        ALU_TOR: lambda: (tor(a, b), 0),
        ALU_TXOR: lambda: (txor(a, b), 0),
        ALU_TNAND: lambda: (tnand(a, b), 0),
        ALU_TNOR: lambda: (tnor(a, b), 0),
        ALU_TXNOR: lambda: (txnor(a, b), 0),
        ALU_TMAJ: lambda: (tmaj(a, b, cin), 0),
        ALU_TSTEP: lambda: (tstep(a), 0),
        ALU_TADD: lambda: tadd(a, b, cin),
        ALU_TSUB: lambda: tsub(a, b, cin),
        ALU_TMUL: lambda: (tmul(a, b), 0),
        ALU_TCMP: lambda: (tcmp(a, b), 0),
        ALU_TSHL: lambda: (tshl(a), 0),
        ALU_TSHR: lambda: (tshr(a), 0),
    }
    return operations.get(op, lambda: (0, 0))()


# =============================================================================
# WORD18 / ALU18
# =============================================================================

def word18_cero():
    return [0] * WORD_SIZE


def word18_valido(word):
    return (
        isinstance(word, (list, tuple))
        and len(word) == WORD_SIZE
        and all(trit_valido(x) for x in word)
    )


def word18_from_string(text):
    values = [{"+" : 1, "-" : -1, "0": 0}[c] for c in text if c in "+-0"]
    return (values + [0] * WORD_SIZE)[:WORD_SIZE]


def word18_to_string(word):
    if not word18_valido(word):
        return "?" * WORD_SIZE
    return "".join("+" if x == 1 else "-" if x == -1 else "0" for x in word)


def word18_to_integer(word):
    if not word18_valido(word):
        return 0
    value = 0
    for trit in word:
        value = value * 3 + trit
    return value


def integer_to_word18(value):
    limit = (3 ** WORD_SIZE - 1) // 2
    value = max(-limit, min(limit, int(value)))
    result = [0] * WORD_SIZE
    for i in range(WORD_SIZE - 1, -1, -1):
        remainder = value % 3
        if remainder == 0:
            trit = 0
        elif remainder == 1:
            trit = 1
        else:
            trit = -1
            value += 1
        result[i] = trit
        value //= 3
    return result


class ALU18:
    @staticmethod
    def tnot(word):
        return [tnot(x) for x in word]

    @staticmethod
    def tand(a, b):
        return [tand(x, y) for x, y in zip(a, b)]

    @staticmethod
    def tor(a, b):
        return [tor(x, y) for x, y in zip(a, b)]

    @staticmethod
    def txor(a, b):
        return [txor(x, y) for x, y in zip(a, b)]

    @staticmethod
    def tnand(a, b):
        return [tnand(x, y) for x, y in zip(a, b)]

    @staticmethod
    def tnor(a, b):
        return [tnor(x, y) for x, y in zip(a, b)]

    @staticmethod
    def txnor(a, b):
        return [txnor(x, y) for x, y in zip(a, b)]

    @staticmethod
    def tmaj(a, b, c):
        return [tmaj(x, y, z) for x, y, z in zip(a, b, c)]

    @staticmethod
    def tstep(word):
        return [tstep(x) for x in word]

    @staticmethod
    def tmul(a, b):
        return [tmul(x, y) for x, y in zip(a, b)]

    @staticmethod
    def tcmp(a, b):
        for x, y in zip(a, b):
            result = tcmp(x, y)
            if result:
                return result
        return 0

    @staticmethod
    def tadd(a, b, cin=0):
        result, carry = [0] * WORD_SIZE, cin
        for i in range(WORD_SIZE - 1, -1, -1):
            result[i], carry = tadd(a[i], b[i], carry)
        return result, carry

    @staticmethod
    def tsub(a, b, cin=0):
        result, carry = [0] * WORD_SIZE, cin
        for i in range(WORD_SIZE - 1, -1, -1):
            result[i], carry = tsub(a[i], b[i], carry)
        return result, carry

    @staticmethod
    def tshl(word):
        return word[1:] + [0] if word18_valido(word) else word18_cero()

    @staticmethod
    def tshr(word):
        return [0] + word[:-1] if word18_valido(word) else word18_cero()

    @staticmethod
    def execute(op, a, b=None, c=None, cin=0):
        b = word18_cero() if b is None else b
        c = word18_cero() if c is None else c
        if not word18_valido(a) or not word18_valido(b) or not word18_valido(c):
            return word18_cero(), 0
        if op == ALU_TNOT:
            return ALU18.tnot(a), 0
        if op == ALU_TAND:
            return ALU18.tand(a, b), 0
        if op == ALU_TOR:
            return ALU18.tor(a, b), 0
        if op == ALU_TXOR:
            return ALU18.txor(a, b), 0
        if op == ALU_TNAND:
            return ALU18.tnand(a, b), 0
        if op == ALU_TNOR:
            return ALU18.tnor(a, b), 0
        if op == ALU_TXNOR:
            return ALU18.txnor(a, b), 0
        if op == ALU_TMAJ:
            return ALU18.tmaj(a, b, c), 0
        if op == ALU_TSTEP:
            return ALU18.tstep(a), 0
        if op == ALU_TADD:
            return ALU18.tadd(a, b, cin)
        if op == ALU_TSUB:
            return ALU18.tsub(a, b, cin)
        if op == ALU_TMUL:
            return ALU18.tmul(a, b), 0
        if op == ALU_TCMP:
            return ALU18.tcmp(a, b), 0
        if op == ALU_TSHL:
            return ALU18.tshl(a), 0
        if op == ALU_TSHR:
            return ALU18.tshr(a), 0
        return word18_cero(), 0


# =============================================================================
# CÓDIGO GENÉTICO Y CONVERSIÓN DNA
# =============================================================================

GENETIC_CODE = {
    "TTT":"Fenilalanina","TTC":"Fenilalanina","TTA":"Leucina","TTG":"Leucina",
    "CTT":"Leucina","CTC":"Leucina","CTA":"Leucina","CTG":"Leucina",
    "ATT":"Isoleucina","ATC":"Isoleucina","ATA":"Isoleucina","ATG":"Metionina",
    "GTT":"Valina","GTC":"Valina","GTA":"Valina","GTG":"Valina",
    "TCT":"Serina","TCC":"Serina","TCA":"Serina","TCG":"Serina",
    "CCT":"Prolina","CCC":"Prolina","CCA":"Prolina","CCG":"Prolina",
    "ACT":"Treonina","ACC":"Treonina","ACA":"Treonina","ACG":"Treonina",
    "GCT":"Alanina","GCC":"Alanina","GCA":"Alanina","GCG":"Alanina",
    "TAT":"Tirosina","TAC":"Tirosina","TAA":"STOP","TAG":"STOP",
    "CAT":"Histidina","CAC":"Histidina","CAA":"Glutamina","CAG":"Glutamina",
    "AAT":"Asparagina","AAC":"Asparagina","AAA":"Lisina","AAG":"Lisina",
    "GAT":"Ácido aspártico","GAC":"Ácido aspártico",
    "GAA":"Ácido glutámico","GAG":"Ácido glutámico",
    "TGT":"Cisteína","TGC":"Cisteína","TGA":"STOP","TGG":"Triptófano",
    "CGT":"Arginina","CGC":"Arginina","CGA":"Arginina","CGG":"Arginina",
    "AGT":"Serina","AGC":"Serina","AGA":"Arginina","AGG":"Arginina",
    "GGT":"Glicina","GGC":"Glicina","GGA":"Glicina","GGG":"Glicina",
}


def nucleotide_to_trit(base):
    return {"A": 1, "T": -1, "U": -1, "C": 0, "G": 0}.get(base.upper(), 0)


def dna_to_trits(sequence):
    return [nucleotide_to_trit(b) for b in sequence if b.upper() in "ATCGU"]


def dna_to_word18(sequence):
    return (dna_to_trits(sequence) + [0] * WORD_SIZE)[:WORD_SIZE]


def codon_to_aminoacid(codon):
    return GENETIC_CODE.get(codon.upper().replace("U", "T"), "DESCONOCIDO")


def translate_dna(sequence):
    sequence = sequence.upper().replace("U", "T")
    return [(sequence[i:i+3], codon_to_aminoacid(sequence[i:i+3]))
            for i in range(0, len(sequence)-2, 3)]


# =============================================================================
# NEXUSISA
# =============================================================================

class NexusInstruction:
    OPCODE_SIZE, DEST_SIZE, SRC_SIZE, IMM_SIZE = 3, 3, 3, 9

    def __init__(self, opcode=0, dest=0, src=0, imm=0):
        self.opcode, self.dest, self.src, self.imm = opcode, dest, src, imm

    @staticmethod
    def encode_field(value, size):
        limit = (3 ** size - 1) // 2
        value = max(-limit, min(limit, int(value)))
        return integer_to_word18(value)[-size:]

    @staticmethod
    def decode_field(trits):
        value = 0
        for trit in trits:
            value = value * 3 + trit
        return value

    def encode(self):
        return (
            self.encode_field(self.opcode, self.OPCODE_SIZE)
            + self.encode_field(self.dest, self.DEST_SIZE)
            + self.encode_field(self.src, self.SRC_SIZE)
            + self.encode_field(self.imm, self.IMM_SIZE)
        )

    @classmethod
    def decode(cls, word):
        if not word18_valido(word):
            return cls()
        p = 0
        op = cls.decode_field(word[p:p+cls.OPCODE_SIZE]); p += cls.OPCODE_SIZE
        dest = cls.decode_field(word[p:p+cls.DEST_SIZE]); p += cls.DEST_SIZE
        src = cls.decode_field(word[p:p+cls.SRC_SIZE]); p += cls.SRC_SIZE
        imm = cls.decode_field(word[p:p+cls.IMM_SIZE])
        return cls(op, dest, src, imm)


class NexusCPU:
    def __init__(self):
        self.registers = {i: word18_cero() for i in range(3)}
        self.cout = 0
        self.last_result = word18_cero()

    def set_register(self, index, word):
        if index in self.registers and word18_valido(word):
            self.registers[index] = list(word)

    def get_register(self, index):
        return list(self.registers.get(index, word18_cero()))

    def execute(self, instruction):
        op, dest, src = instruction.opcode, instruction.dest, instruction.src
        a, b = self.get_register(dest), self.get_register(src)
        immediate = integer_to_word18(instruction.imm)
        if op in (ALU_TNOT, ALU_TSTEP, ALU_TSHL, ALU_TSHR):
            result, cout = ALU18.execute(op, a)
        elif op in (ALU_TADD, ALU_TSUB):
            result, cout = ALU18.execute(op, a, b, cin=0)
        elif op == ALU_TMAJ:
            result, cout = ALU18.execute(op, a, b, immediate)
        else:
            result, cout = ALU18.execute(op, a, b)
        self.cout, self.last_result = cout, list(result)
        if dest in self.registers:
            self.registers[dest] = list(result)
        return result, cout


# =============================================================================
# VCF / WGS REAL
# =============================================================================

def check_vcf(path):
    """Comprueba existencia, gzip, cabecera VCF y al menos una variante válida."""
    if not os.path.isfile(path):
        return False, "archivo no encontrado"
    try:
        with open(path, "rb") as raw:
            if raw.read(2) != b"\x1f\x8b":
                return False, "el archivo no tiene cabecera gzip"
        has_fileformat = False
        has_columns = False
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as f:
            for line in f:
                line = line.rstrip("\r\n")
                if line.startswith("##fileformat=VCF"):
                    has_fileformat = True
                elif line.startswith("#CHROM"):
                    columns = line.split("\t")
                    has_columns = len(columns) >= 8 and columns[:5] == [
                        "#CHROM", "POS", "ID", "REF", "ALT"
                    ]
                elif line and not line.startswith("#"):
                    fields = line.split("\t")
                    if has_fileformat and has_columns and len(fields) >= 8:
                        try:
                            if int(fields[1]) > 0 and fields[0] and fields[3] and fields[4]:
                                return True, "cabecera y primera variante válidas"
                        except ValueError:
                            continue
        return False, "cabecera o variantes válidas no encontradas"
    except (OSError, EOFError, ValueError, UnicodeError) as exc:
        return False, f"error leyendo VCF: {exc}"


class WholeGenomeCodonProcessor:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.chromosomes_data = {}
        self.ranked_coherence = []
        self.ranked_entropy = []
        self.headers = []
        self.variant_count = 0
        self.vcf_ok, self.vcf_message = check_vcf(vcf_path)
        self.load_and_process()

    def load_and_process(self):
        target_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
        temp_data = {c: [] for c in target_chroms}
        if not self.vcf_ok:
            self.headers.append("VCF NO: " + self.vcf_message)
            self.chromosomes_data = temp_data
            self.compute_metrics()
            return

        try:
            with gzip.open(self.vcf_path, "rt", encoding="utf-8", errors="strict") as f:
                for line in f:
                    if line.startswith("#"):
                        if line.startswith("##fileformat=") or line.startswith("#CHROM"):
                            self.headers.append(line.strip())
                        continue
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) < 8:
                        continue
                    chrom = parts[0].lower().replace("chr", "").strip().upper()
                    if chrom not in temp_data:
                        continue
                    ref, alt = parts[3], parts[4]
                    if len(ref) == 1 and len(alt) == 1 and ref.upper() in "ACGT" and alt.upper() in "ACGT":
                        if len(temp_data[chrom]) < 300:
                            temp_data[chrom].extend([
                                nucleotide_to_trit(ref), nucleotide_to_trit(alt)
                            ])
                        self.variant_count += 1
        except (OSError, EOFError, ValueError, UnicodeError) as exc:
            self.vcf_ok = False
            self.vcf_message = f"error leyendo VCF: {exc}"
            self.headers.append("VCF NO: " + self.vcf_message)
            temp_data = {c: [] for c in target_chroms}

        self.chromosomes_data = temp_data
        self.compute_metrics()

    def compute_metrics(self):
        results = []
        for chrom, trits in self.chromosomes_data.items():
            total = len(trits)
            if not total:
                continue
            agreements = sum(1 for i in range(total - 1) if trits[i] == trits[i + 1])
            coherence = min(0.9850, round(0.5 + agreements / max(1, total - 1) * 0.4127, 4))
            probs = [trits.count(v) / total for v in (0, 1, -1)]
            entropy = round(-sum(p * math.log2(p) for p in probs if p > 0), 4)
            results.append({"chrom": chrom, "coherence": coherence, "entropy": entropy, "total": total})
        self.ranked_coherence = sorted(results, key=lambda x: x["coherence"], reverse=True)
        self.ranked_entropy = sorted(results, key=lambda x: x["entropy"], reverse=True)


# =============================================================================
# INFORME DE SALIDA
# =============================================================================

def write_report(processor, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "trinova_vcf_report.txt")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "w", encoding="utf-8") as f:
        f.write("TRINOVA.OS - INFORME DE PROCESAMIENTO VCF\n")
        f.write(f"Fecha: {now}\n")
        f.write(f"VCF: {processor.vcf_path}\n")
        f.write(f"Estado VCF: {'VCF OK' if processor.vcf_ok else 'VCF NO'}\n")
        f.write(f"Detalle: {processor.vcf_message}\n")
        f.write(f"Variantes SNV de cromosomas 1-22/X/Y contabilizadas: {processor.variant_count}\n")
        f.write("\nNOTA: métricas experimentales, no diagnóstico clínico.\n\n")
        f.write("TOP COHERENCIA\n")
        for row in processor.ranked_coherence[:10]:
            f.write(f"Chr {row['chrom']}: C={row['coherence']} trits={row['total']}\n")
        f.write("\nTOP ENTROPIA\n")
        for row in processor.ranked_entropy[:10]:
            f.write(f"Chr {row['chrom']}: S={row['entropy']} trits={row['total']}\n")
    return path


# =============================================================================
# WIN32 GUI
# =============================================================================

if os.name == "nt":
    HCURSOR = ctypes.c_void_p
    HICON = ctypes.c_void_p
    HBRUSH = ctypes.c_void_p
    HMENU = ctypes.c_void_p
    HWND = ctypes.c_void_p
    HINSTANCE = ctypes.c_void_p
    LRESULT = ctypes.c_int64 if ctypes.sizeof(ctypes.c_void_p) == 8 else ctypes.c_int32

    class BITMAPINFOHEADER(ctypes.Structure):
        _fields_ = [
            ("biSize", wintypes.DWORD), ("biWidth", wintypes.LONG),
            ("biHeight", wintypes.LONG), ("biPlanes", wintypes.WORD),
            ("biBitCount", wintypes.WORD), ("biCompression", wintypes.DWORD),
            ("biSizeImage", wintypes.DWORD), ("biXPelsPerMeter", wintypes.LONG),
            ("biYPelsPerMeter", wintypes.LONG), ("biClrUsed", wintypes.DWORD),
            ("biClrImportant", wintypes.DWORD),
        ]

    WNDPROC = ctypes.WINFUNCTYPE(LRESULT, HWND, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM)

    class WNDCLASSW(ctypes.Structure):
        _fields_ = [
            ("style", wintypes.UINT), ("lpfnWndProc", WNDPROC),
            ("cbClsExtra", ctypes.c_int), ("cbWndExtra", ctypes.c_int),
            ("hInstance", HINSTANCE), ("hIcon", HICON), ("hCursor", HCURSOR),
            ("hbrBackground", HBRUSH), ("lpszMenuName", wintypes.LPCWSTR),
            ("lpszClassName", wintypes.LPCWSTR),
        ]


class TriNovaCanvas:
    def __init__(self, width=1280, height=800):
        self.width, self.height = width, height
        self.stride = (width * 3 + 3) & ~3
        self.vram = bytearray(self.stride * height)
        self.colors = {
            "desktop": bytes((12, 8, 6)), "win_bg": bytes((18, 14, 10)),
            "win_header": bytes((35, 24, 14)), "win_header_active": bytes((150, 80, 10)),
            "border_cyan": bytes((160, 110, 0)), "border_amber": bytes((0, 120, 220)),
        }

    def clear(self, key="desktop"):
        color = self.colors[key]
        row = color * self.width + b"\x00" * (self.stride - self.width * 3)
        self.vram[:] = row * self.height

    def set_pixel(self, x, y, color):
        if 0 <= x < self.width and 0 <= y < self.height:
            offset = int(y) * self.stride + int(x) * 3
            self.vram[offset:offset+3] = color

    def draw_line(self, x0, y0, x1, y1, color):
        x0, y0, x1, y1 = map(int, (x0, y0, x1, y1))
        dx, dy = abs(x1-x0), -abs(y1-y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            self.set_pixel(x0, y0, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def draw_rect_filled(self, x, y, w, h, color):
        x0, y0 = max(0, int(x)), max(0, int(y))
        x1, y1 = min(self.width, int(x+w)), min(self.height, int(y+h))
        if x0 >= x1 or y0 >= y1:
            return
        segment = color * (x1-x0)
        for py in range(y0, y1):
            offset = py*self.stride + x0*3
            self.vram[offset:offset+len(segment)] = segment

    def draw_rect_outline(self, x, y, w, h, color):
        self.draw_line(x, y, x+w, y, color)
        self.draw_line(x, y+h, x+w, y+h, color)
        self.draw_line(x, y, x, y+h, color)
        self.draw_line(x+w, y, x+w, y+h, color)


class Window:
    def __init__(self, title, x, y, w, h, app_type):
        self.title, self.x, self.y, self.w, self.h, self.app_type = title, x, y, w, h, app_type
        self.minimized = False
        self.focused = False

    def is_inside_window(self, mx, my):
        return self.x <= mx <= self.x+self.w and self.y <= my <= self.y+self.h

    def is_inside_header(self, mx, my):
        return self.x <= mx <= self.x+self.w and self.y <= my <= self.y+25

    def is_inside_close_btn(self, mx, my):
        return self.x+self.w-22 <= mx <= self.x+self.w-5 and self.y+4 <= my <= self.y+21


class TriNovaDesktopOS:
    def __init__(self, width=1280, height=800):
        if os.name != "nt":
            raise OSError("TriNova.OS GUI requiere Windows (Win32).")
        self.width, self.height = width, height
        self.canvas = TriNovaCanvas(width, height)
        self.mouse_x, self.mouse_y = 640, 400
        self.dragging_win = None
        self.drag_offset_x = self.drag_offset_y = 0

        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.base_dir = os.path.dirname(VCF_REAL_PATH)
        self.vcf_path = VCF_REAL_PATH
        self.gtf_path = os.path.join(self.base_dir, "gencode.v50.annotation.gtf.gz")
        self.fasta_path = os.path.join(self.base_dir, "GRCh38.primary_assembly.genome.fa.gz")
        self.output_dir = self.script_dir

        self.genome_processor = WholeGenomeCodonProcessor(self.vcf_path)
        self.report_path = write_report(self.genome_processor, self.output_dir)
        print("[VCF] " + ("VCF OK" if self.genome_processor.vcf_ok else "VCF NO"))
        print("[VCF] " + self.genome_processor.vcf_message)
        print("[VCF] Ruta: " + self.vcf_path)
        print("[OUTPUT] Informe: " + self.report_path)
        for label, path in (("GTF", self.gtf_path), ("FASTA", self.fasta_path)):
            print(f"[{'RUTA OK' if os.path.isfile(path) else 'AVISO'}] {label}: {path}")

        self.cpu = NexusCPU()
        self.demo_dna = "ATGGCTTTTTGA"
        self.demo_word18 = dna_to_word18(self.demo_dna)
        self.cpu.set_register(1, self.demo_word18)
        self.cpu.set_register(0, integer_to_word18(1))
        self.demo_instruction = NexusInstruction(ALU_TADD, 1, 0, 0)
        self.demo_instruction_word = self.demo_instruction.encode()
        self.demo_result, self.demo_cout = self.cpu.execute(self.demo_instruction)

        self.windows = [
            Window("BIBLIOTECA GENETICA & ALU18", 20, 40, 400, 360, "LIBRARY_VIEW"),
            Window("TOP 10: COHERENCIA GENOMICA (C)", 430, 40, 410, 360, "TOP_COHERENCE"),
            Window("TOP 10: ENTROPIA TERNARIA (S)", 850, 40, 410, 360, "TOP_ENTROPY"),
            Window("TRINOVA.OS - CONSOLA Y FLUJO VCF WGS", 20, 410, 1240, 330, "TERMINAL"),
        ]
        self.focus_window(self.windows[0])
        self._init_win32()

    def focus_window(self, target):
        for win in self.windows:
            win.focused = False
        target.focused = True
        self.windows.remove(target)
        self.windows.append(target)

    def _init_win32(self):
        self.user32, self.gdi32, self.kernel32 = ctypes.windll.user32, ctypes.windll.gdi32, ctypes.windll.kernel32
        self.user32.DefWindowProcW.argtypes = [HWND, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM]
        self.user32.DefWindowProcW.restype = LRESULT
        self.hInstance = self.kernel32.GetModuleHandleW(None)
        self.wndProc = WNDPROC(self._wnd_proc)
        wc = WNDCLASSW()
        wc.style, wc.lpfnWndProc, wc.hInstance = 3, self.wndProc, self.hInstance
        wc.hIcon = self.user32.LoadIconW(0, ctypes.c_wchar_p(32512))
        wc.hCursor = self.user32.LoadCursorW(0, ctypes.c_wchar_p(32512))
        wc.lpszClassName = "TriNovaUnifiedClass"
        self.user32.RegisterClassW(ctypes.byref(wc))
        rect = wintypes.RECT(0, 0, self.width, self.height)
        self.user32.AdjustWindowRect(ctypes.byref(rect), 0x00CF0000, False)
        self.hwnd = self.user32.CreateWindowExW(
            0, "TriNovaUnifiedClass", "TriNova.OS - Biblioteca Genetica y Genoma",
            0x00CF0000 | 0x10000000, 50, 50, rect.right-rect.left, rect.bottom-rect.top,
            0, 0, self.hInstance, None
        )
        if not self.hwnd:
            raise ctypes.WinError()
        self.hdc = self.user32.GetDC(self.hwnd)
        self.bmi = BITMAPINFOHEADER()
        self.bmi.biSize = ctypes.sizeof(BITMAPINFOHEADER)
        self.bmi.biWidth, self.bmi.biHeight = self.width, -self.height
        self.bmi.biPlanes, self.bmi.biBitCount, self.bmi.biCompression = 1, 24, 0
        self.hfont = self.gdi32.CreateFontW(16, 0, 0, 0, 400, False, False, False, 1, 0, 0, 0, 0, "Consolas")
        self.gdi32.SelectObject(self.hdc, self.hfont)
        self.gdi32.SetBkMode(self.hdc, 1)

    def _wnd_proc(self, hwnd, msg, wparam, lparam):
        if msg == 0x0014:  # WM_ERASEBKGND: no borrar antes del siguiente frame
            return 1
        if msg == 0x000F:  # WM_PAINT: validar región para evitar repintado pendiente
            ps = ctypes.create_string_buffer(72)
            self.user32.BeginPaint(hwnd, ps)
            self.user32.EndPaint(hwnd, ps)
            return 0
        if msg == 0x0200:  # WM_MOUSEMOVE
            self.mouse_x, self.mouse_y = lparam & 0xFFFF, (lparam >> 16) & 0xFFFF
            if self.dragging_win:
                self.dragging_win.x = self.mouse_x - self.drag_offset_x
                self.dragging_win.y = self.mouse_y - self.drag_offset_y
        elif msg == 0x0201:  # WM_LBUTTONDOWN
            self.mouse_x, self.mouse_y = lparam & 0xFFFF, (lparam >> 16) & 0xFFFF
            for win in reversed(self.windows):
                if not win.minimized and win.is_inside_window(self.mouse_x, self.mouse_y):
                    self.focus_window(win)
                    if win.is_inside_close_btn(self.mouse_x, self.mouse_y):
                        win.minimized = True
                    elif win.is_inside_header(self.mouse_x, self.mouse_y):
                        self.dragging_win = win
                        self.drag_offset_x, self.drag_offset_y = self.mouse_x-win.x, self.mouse_y-win.y
                    break
        elif msg == 0x0202:
            self.dragging_win = None
        elif msg == 0x0002:
            self.user32.PostQuitMessage(0)
            return 0
        return self.user32.DefWindowProcW(hwnd, msg, wparam, lparam)

    def draw_text(self, text, x, y, color_rgb=(255, 255, 255)):
        self.gdi32.SetTextColor(self.hdc, color_rgb[0] | (color_rgb[1] << 8) | (color_rgb[2] << 16))
        self.gdi32.TextOutW(self.hdc, int(x), int(y), text, len(text))

    def render_desktop(self):
        self.canvas.clear("desktop")
        for win in self.windows:
            if win.minimized:
                continue
            bg = self.canvas.colors["win_bg"]
            header = self.canvas.colors["win_header_active"] if win.focused else self.canvas.colors["win_header"]
            border = self.canvas.colors["border_amber"] if win.focused else self.canvas.colors["border_cyan"]
            self.canvas.draw_rect_filled(win.x, win.y, win.w, win.h, bg)
            self.canvas.draw_rect_filled(win.x, win.y, win.w, 25, header)
            self.canvas.draw_rect_outline(win.x, win.y, win.w, win.h, border)

        buffer_type = ctypes.c_ubyte * len(self.canvas.vram)
        self.gdi32.StretchDIBits(
            self.hdc, 0, 0, self.width, self.height, 0, 0, self.width, self.height,
            buffer_type.from_buffer(self.canvas.vram), ctypes.byref(self.bmi), 0, 0x00CC0020
        )

        for win in self.windows:
            if win.minimized:
                continue
            self.draw_text(win.title, win.x+8, win.y+4, (255, 255, 255))
            if win.app_type == "LIBRARY_VIEW":
                lines = [
                    "NIVEL 1: DNA['A','C','G','T'] RNA['A','C','G','U']",
                    "NIVEL 2: Codon: AUG -> Anticodon: UAC",
                    "NIVEL 3: AUG -> Metionina",
                    f"GEN: SISTEMA ({self.demo_dna})",
                    f"WORD18: {word18_to_string(self.demo_word18)}",
                    f"ALU18 TNOT: {word18_to_string(ALU18.tnot(self.demo_word18))}",
                    f"TADD18 COUT: {self.demo_cout}",
                ]
                y = win.y + 38
                for line in lines:
                    self.draw_text(line, win.x+10, y, (255, 255, 0))
                    y += 24
            elif win.app_type in ("TOP_COHERENCE", "TOP_ENTROPY"):
                coherence = win.app_type == "TOP_COHERENCE"
                data = self.genome_processor.ranked_coherence if coherence else self.genome_processor.ranked_entropy
                y = win.y + 40
                for idx, item in enumerate(data[:10]):
                    metric = f"C:{item['coherence']}" if coherence else f"S:{item['entropy']}"
                    color = (239, 204, 34) if coherence else (200, 68, 239)
                    self.draw_text(f"#{idx+1} Chr {item['chrom']} | {metric}", win.x+10, y, color)
                    y += 25
                if not data:
                    self.draw_text("Sin datos VCF reales disponibles.", win.x+10, y, (255, 180, 100))
            elif win.app_type == "TERMINAL":
                p = self.genome_processor
                top = p.ranked_coherence[0]["chrom"] if p.ranked_coherence else "N/A"
                logs = [
                    "[OK] Biblioteca genética TriNova y ALU18 cargadas.",
                    f"[VCF] {'VCF OK' if p.vcf_ok else 'VCF NO'}",
                    f"[INFO] {p.vcf_message}",
                    f"[INFO] WORD18 DNA: {word18_to_string(self.demo_word18)}",
                    f"[INFO] NexusISA: {word18_to_string(self.demo_instruction_word)}",
                    f"[INFO] OPCODE: {self.demo_instruction.opcode} ({ALU_NAMES.get(self.demo_instruction.opcode, 'UNKNOWN')})",
                    f"[ALU18] RESULT: {word18_to_string(self.demo_result)}",
                    f"[ALU18] COUT: {self.demo_cout}",
                    f"[INFO] Cromosoma con mayor coherencia: Chr {top}",
                    f"[INFO] Variantes SNV contadas: {p.variant_count}",
                    f"[OUTPUT] {os.path.basename(self.report_path)}",
                    "[READY] Sistema + WGS + ALU18 sincronizados.",
                ]
                y = win.y + 45
                for log in logs:
                    self.draw_text(log, win.x+10, y, (74, 222, 128))
                    y += 22

    def run(self):
        msg = wintypes.MSG()
        self.user32.ShowWindow(self.hwnd, 5)
        self.user32.UpdateWindow(self.hwnd)
        while True:
            while self.user32.PeekMessageW(ctypes.byref(msg), 0, 0, 0, 1):
                if msg.message == 0x0012:
                    return
                self.user32.TranslateMessage(ctypes.byref(msg))
                self.user32.DispatchMessageW(ctypes.byref(msg))
            self.render_desktop()
            time.sleep(0.016)


# =============================================================================
# SELF TEST
# =============================================================================

def self_test():
    print("=" * 70)
    print("TRINOVA.OS - SELF TEST")
    print("=" * 70)
    print("[ALU-TRIT]")
    print("TNOT(+1)       =", tnot(1))
    print("TAND(+1,-1)    =", tand(1, -1))
    print("TOR(+1,-1)     =", tor(1, -1))
    print("TMAJ(-1,+1,+1) =", tmaj(-1, 1, 1))
    print("TADD(+1,+1,0)  =", tadd(1, 1, 0))
    dna = "ATGGCTTTTTGA"
    word = dna_to_word18(dna)
    print("\n[BIOLOGIA]")
    print("DNA            :", dna)
    for codon, amino in translate_dna(dna):
        print(f"CODON          : {codon} -> {amino}")
    print("\n[WORD18]")
    print("WORD18         :", word18_to_string(word))
    print("LONGITUD       :", len(word))
    print("VALOR TERNARIO :", word18_to_integer(word))
    print("\n[ALU18]")
    result, cout = ALU18.tadd(word, integer_to_word18(1), 0)
    print("TADD18(A,B)    :", word18_to_string(result))
    print("COUT           :", cout)
    instruction = NexusInstruction(ALU_TADD, 1, 0, 0)
    cpu = NexusCPU()
    cpu.set_register(1, word)
    cpu.set_register(0, integer_to_word18(1))
    result, cout = cpu.execute(instruction)
    print("\n[NEXUSISA]")
    print("INSTRUCTION    :", word18_to_string(instruction.encode()))
    print("RESULT         :", word18_to_string(result))
    print("COUT           :", cout)
    print("\n[OK] SELF TEST FINALIZADO")
    print("=" * 70)


if __name__ == "__main__":
    self_test()
    print("\nIniciando TriNova.OS...")
    try:
        app = TriNovaDesktopOS()
        app.run()
    except Exception as exc:
        print(f"[ERROR] No se pudo iniciar TriNova.OS: {exc}")
        raise
