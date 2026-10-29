# -*- coding: utf-8 -*-
"""
TRINOVA / MVA-DATA — Visor genómico
Windows 11 + Python 3.13
Solo biblioteca estándar de Python.
"""

import gzip
import json
import os
import queue
import random
import re
import threading
import tkinter as tk

from tkinter import ttk, messagebox
from collections import Counter


# ============================================================
# 1. CONFIGURACIÓN Y RUTAS
# ============================================================

BASE_DIR = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data"
)

CLINICAL_PATH = os.path.join(BASE_DIR, "clinico.json")
GTF_PATH = os.path.join(BASE_DIR, "gencode.v50.annotation.gtf.gz")
FASTA_PATH = os.path.join(
    BASE_DIR, "GRCh38.primary_assembly.genome.fa.gz"
)
VCF_PATH = os.path.join(
    BASE_DIR, "WGS_proband01.vcf.gz"
)
VCF_INDEX_PATH = VCF_PATH + ".tbi"

TARGET_GENE = "BUB1B"
MAX_PLOT_POINTS = 12000

PLOT_BACKGROUND = "#101722"
PLOT_TEXT_COLOR = "#dbe6f5"
OTHER_VARIANT_COLOR = "#8b98a8"

GENE_COLORS = {
    "BUB1B": "#ff5c7a",
    "GOLGA8IP": "#42d6c5",
    "NIPA1": "#ffc857",
    "SNRPN": "#b69cff",
    "ARHGAP11A-DT": "#70a7ff",
}

FALLBACK_GENE = {
    "chrom": "15",
    "start": 40160984,
    "end": 40221137,
    "strand": "+",
    "gene": "BUB1B",
}

CANDIDATE_VARIANTS = [
    {
        "gene": "BUB1B", "chrom": "15",
        "pos": 40209701, "ref": "T", "alt": "G",
        "epcr": 0.98, "finding_type": "primary",
        "notes": "Candidato principal; consecuencia por comprobar.",
    },
    {
        "gene": "BUB1B", "chrom": "15",
        "pos": 40220612, "ref": "T", "alt": "G",
        "epcr": 0.98, "finding_type": "primary",
        "notes": "Candidato principal; genotipo y fase por comprobar.",
    },
    {
        "gene": "GOLGA8IP", "chrom": "15",
        "pos": 22587629, "ref": "T", "alt": "C",
        "epcr": 0.85, "finding_type": "secondary",
        "notes": "Candidato secundario; genotipo y fase por comprobar.",
    },
    {
        "gene": "GOLGA8IP", "chrom": "15",
        "pos": 22608113, "ref": "C", "alt": "T",
        "epcr": 0.85, "finding_type": "secondary",
        "notes": "Candidato secundario; genotipo y fase por comprobar.",
    },
    {
        "gene": "NIPA1", "chrom": "15",
        "pos": 22814504, "ref": "A", "alt": "G",
        "epcr": 0.80, "finding_type": "secondary",
        "notes": "Consecuencia y evidencia clínica por comprobar.",
    },
    {
        "gene": "NIPA1", "chrom": "15",
        "pos": 22826077, "ref": "G", "alt": "GA",
        "epcr": 0.80, "finding_type": "secondary",
        "notes": "Inserción candidata; representación por comprobar.",
    },
    {
        "gene": "SNRPN", "chrom": "15",
        "pos": 24958486, "ref": "GTTTTTTT", "alt": "G",
        "epcr": 0.75, "finding_type": "secondary",
        "notes": "Revisar región y representación de la variante.",
    },
    {
        "gene": "SNRPN", "chrom": "15",
        "pos": 24974365, "ref": "T", "alt": "C",
        "epcr": 0.75, "finding_type": "secondary",
        "notes": "Consecuencia por comprobar.",
    },
    {
        "gene": "ARHGAP11A-DT", "chrom": "15",
        "pos": 32573465, "ref": "T", "alt": "C",
        "epcr": 0.65, "finding_type": "secondary",
        "notes": "Anotación funcional por comprobar.",
    },
    {
        "gene": "ARHGAP11A-DT", "chrom": "15",
        "pos": 32574480, "ref": "G", "alt": "C",
        "epcr": 0.65, "finding_type": "secondary",
        "notes": "Anotación funcional por comprobar.",
    },
]

CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

CHROM_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559,
    "4": 190214555, "5": 181538259, "6": 170805979,
    "7": 159345973, "8": 145138636, "9": 138394717,
    "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285,
    "19": 58617616, "20": 64444167, "21": 46709983,
    "22": 50818468, "X": 156040895, "Y": 57227415,
    "MT": 16569,
}


# ============================================================
# 2. UTILIDADES
# ============================================================

def open_text_auto(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8",
                         errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def normalize_chrom(chrom):
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in ("M", "MT"):
        return "MT"
    if chrom.upper() in ("X", "Y"):
        return chrom.upper()
    return chrom


def safe_int(value):
    try:
        if value in ("", ".", None):
            return None
        return int(value)
    except (ValueError, TypeError):
        return None


def safe_float(value):
    try:
        if value in ("", ".", None):
            return None
        return float(value)
    except (ValueError, TypeError):
        return None


def format_number(value):
    if value is None:
        return "N/D"
    try:
        return f"{int(value):,}"
    except (ValueError, TypeError):
        return str(value)


def parse_gtf_attributes(text):
    return {
        match.group(1): match.group(2)
        for match in re.finditer(r'(\S+)\s+"([^"]*)"', text)
    }


def parse_info(text):
    result = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = True
    return result


def parse_sample(format_text, sample_text):
    keys = format_text.split(":") if format_text else []
    values = sample_text.split(":") if sample_text else []
    sample = dict(zip(keys, values))
    return {
        "GT": sample.get("GT", "N/D"),
        "DP": safe_int(sample.get("DP")),
        "GQ": safe_int(sample.get("GQ")),
    }


def parse_vcf_record(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return None

    pos = safe_int(fields[1])
    if pos is None:
        return None

    info = parse_info(fields[7])
    sample = parse_sample(
        fields[8] if len(fields) > 8 else "",
        fields[9] if len(fields) > 9 else "",
    )

    af = None
    if info.get("AF") not in (None, True):
        af = safe_float(str(info["AF"]).split(",")[0])

    return {
        "chrom": normalize_chrom(fields[0]),
        "pos": pos,
        "id": fields[2],
        "ref": fields[3].upper(),
        "alt": fields[4].upper(),
        "qual": safe_float(fields[5]),
        "filter": fields[6],
        "af": af,
        "dp": sample["DP"],
        "gq": sample["GQ"],
        "gt": sample["GT"],
        "info": info,
    }


def candidate_key(candidate):
    return (
        normalize_chrom(candidate["chrom"]),
        candidate["pos"],
        candidate["ref"].upper(),
        candidate["alt"].upper(),
    )


def record_matches_candidate(record, candidate):
    """Compara cromosoma, posición, REF y cada ALT del registro."""
    alts = record["alt"].upper().split(",")

    return (
        normalize_chrom(record["chrom"])
        == normalize_chrom(candidate["chrom"])
        and record["pos"] == candidate["pos"]
        and record["ref"].upper() == candidate["ref"].upper()
        and candidate["alt"].upper() in alts
    )


def candidate_gene_for_record(record):
    for candidate in CANDIDATE_VARIANTS:
        if record_matches_candidate(record, candidate):
            return candidate["gene"]
    return None


# ============================================================
# 3. JSON CLÍNICO Y HPO
# ============================================================

def load_clinical_json(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No se encuentra el JSON:\n{path}")

    with open(path, "rt", encoding="utf-8-sig") as file:
        data = json.load(file)

    entries = []
    seen = set()

    def walk(obj):
        if isinstance(obj, dict):
            possible_id = (
                obj.get("hpo_id") or obj.get("HPO")
                or obj.get("hpo") or obj.get("term_id")
                or obj.get("id")
            )

            if isinstance(possible_id, str):
                match = re.search(r"HP:\d{7}", possible_id)
                if match:
                    hpo_id = match.group(0)
                    label = (
                        obj.get("label") or obj.get("name")
                        or obj.get("term") or obj.get("description")
                        or "Etiqueta no disponible"
                    )
                    key = (hpo_id, str(label))
                    if key not in seen:
                        seen.add(key)
                        entries.append({
                            "hpo_id": hpo_id,
                            "label": str(label),
                        })

            for value in obj.values():
                walk(value)

        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(data)
    return data, entries


# ============================================================
# 4. BÚSQUEDA DE BUB1B EN GENCODE
# ============================================================

def load_gene_from_gtf(path, gene_name="BUB1B"):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No se encuentra el GTF:\n{path}")

    found = []

    with open_text_auto(path) as file:
        for line in file:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attrs = parse_gtf_attributes(fields[8])
            if attrs.get("gene_name") != gene_name:
                continue

            start, end = safe_int(fields[3]), safe_int(fields[4])
            if start is None or end is None:
                continue

            found.append({
                "chrom": normalize_chrom(fields[0]),
                "start": start,
                "end": end,
                "strand": fields[6],
                "gene": gene_name,
                "gene_id": attrs.get("gene_id", "N/D"),
                "gene_type": (
                    attrs.get("gene_type")
                    or attrs.get("gene_biotype", "N/D")
                ),
            })

    if not found:
        return None

    return max(found, key=lambda item: item["end"] - item["start"])


# ============================================================
# 5. ESCANEO DEL VCF
# ============================================================

def scan_vcf(path, target_chrom, start_pos, end_pos,
             progress_callback=None):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No se encuentra el VCF:\n{path}")

    if start_pos is not None and start_pos < 1:
        raise ValueError("La posición inicial debe ser >= 1.")
    if end_pos is not None and end_pos < 1:
        raise ValueError("La posición final debe ser >= 1.")
    if start_pos is not None and end_pos is not None:
        if start_pos > end_pos:
            raise ValueError("La posición inicial supera la final.")

    rng = random.Random(42)
    reservoir = []
    candidate_matches = {}
    total_records = 0
    filtered_records = 0
    chromosome_counts = Counter()

    with open_text_auto(path) as file:
        for line in file:
            if not line or line.startswith("#"):
                continue

            total_records += 1
            record = parse_vcf_record(line)
            if record is None:
                continue

            # Busca candidatos en TODO el VCF, no solo en el filtro.
            for candidate in CANDIDATE_VARIANTS:
                if record_matches_candidate(record, candidate):
                    candidate_matches[candidate_key(candidate)] = (
                        record.copy()
                    )

            chrom, pos = record["chrom"], record["pos"]

            if target_chrom != "all" and chrom != target_chrom:
                continue
            if start_pos is not None and pos < start_pos:
                continue
            if end_pos is not None and pos > end_pos:
                continue

            filtered_records += 1
            chromosome_counts[chrom] += 1

            if len(reservoir) < MAX_PLOT_POINTS:
                reservoir.append(record)
            else:
                replacement = rng.randrange(filtered_records)
                if replacement < MAX_PLOT_POINTS:
                    reservoir[replacement] = record

            if progress_callback and total_records % 100000 == 0:
                progress_callback(
                    f"Leídos {total_records:,} registros; "
                    f"{filtered_records:,} cumplen el filtro."
                )

    reservoir.sort(key=lambda r: (
        CHROM_ORDER.index(r["chrom"])
        if r["chrom"] in CHROM_ORDER else 99,
        r["pos"],
    ))

    return {
        "total_records": total_records,
        "filtered_records": filtered_records,
        "chromosome_counts": chromosome_counts,
        "records": reservoir,
        "plotted_records": len(reservoir),
        "target_chrom": target_chrom,
        "start_pos": start_pos,
        "end_pos": end_pos,
        "candidate_matches": candidate_matches,
    }


# ============================================================
# 6. LECTURA DE REGIÓN FASTA
# ============================================================

def read_fasta_region(path, target_chrom, start, end,
                      progress_callback=None):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No se encuentra el FASTA:\n{path}")

    target_chrom = normalize_chrom(target_chrom)
    sequence_parts = []
    current_chrom = None
    current_pos = 1
    found_chrom = False

    with open_text_auto(path) as file:
        for line in file:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                current_chrom = normalize_chrom(header)
                current_pos = 1

                if current_chrom == target_chrom:
                    found_chrom = True
                    if progress_callback:
                        progress_callback(
                            f"FASTA: buscando {target_chrom}:"
                            f"{start}-{end}..."
                        )
                continue

            if current_chrom != target_chrom:
                continue

            seq_line = line.strip().upper()
            if not seq_line:
                continue

            line_start = current_pos
            line_end = current_pos + len(seq_line) - 1
            current_pos = line_end + 1

            overlap_start = max(start, line_start)
            overlap_end = min(end, line_end)

            if overlap_start <= overlap_end:
                a = overlap_start - line_start
                b = overlap_end - line_start + 1
                sequence_parts.append(seq_line[a:b])

            if line_end >= end:
                break

    if not found_chrom:
        raise ValueError(f"No se encontró el cromosoma {target_chrom}.")

    sequence = "".join(sequence_parts)
    expected_length = end - start + 1

    if len(sequence) != expected_length:
        raise ValueError(
            f"Longitud FASTA inesperada: esperada {expected_length:,}; "
            f"obtenida {len(sequence):,}. Comprueba referencia y coordenadas."
        )

    return sequence


def validate_bub1b_ref(vcf_path, fasta_path, gene,
                        progress_callback=None):
    chrom = normalize_chrom(gene["chrom"])
    start, end = int(gene["start"]), int(gene["end"])

    reference = read_fasta_region(
        fasta_path, chrom, start, end, progress_callback
    )

    results = []

    with open_text_auto(vcf_path) as file:
        for line in file:
            if not line or line.startswith("#"):
                continue

            record = parse_vcf_record(line)
            if record is None or record["chrom"] != chrom:
                continue
            if not start <= record["pos"] <= end:
                continue

            pos, ref = record["pos"], record["ref"].upper()
            offset = pos - start
            expected_ref = reference[offset:offset + len(ref)]

            if len(expected_ref) != len(ref):
                status = "NO EVALUABLE"
            elif expected_ref == ref:
                status = "COINCIDE"
            else:
                status = "NO COINCIDE"

            results.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": record["alt"],
                "expected_ref": expected_ref,
                "status": status,
                "gt": record["gt"],
                "dp": record["dp"],
                "gq": record["gq"],
            })

    return {
        "gene": gene,
        "variants": results,
        "counts": Counter(item["status"] for item in results),
        "region_length": len(reference),
    }


# ============================================================
# 7. INTERFAZ
# ============================================================

class TriNovaVisualizer:
    def __init__(self, root):
        self.root = root
        self.root.title("TRINOVA | Visor genómico y candidatos")
        self.root.geometry("1350x900")
        self.root.minsize(1000, 700)

        self.result_queue = queue.Queue()
        self.busy = False
        self.results = None
        self.ref_results = None
        self.clinical_data = None
        self.hpo_entries = []
        self.bub1b = FALLBACK_GENE.copy()
        self.gtf_gene_found = False

        self.candidate_matches = {}
        self.scan_completed = False
        self.candidate_tree = None
        self.candidate_window = None

        self._build_ui()
        self._check_paths()

        threading.Thread(
            target=self._load_metadata_worker, daemon=True
        ).start()

        self.root.after(150, self._poll_queue)

    def _build_ui(self):
        main = ttk.Frame(self.root, padding=10)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main,
            text="TRINOVA — Visor de variantes genómicas",
            font=("Segoe UI", 16, "bold"),
        ).pack(anchor="w", pady=(0, 8))

        filters = ttk.LabelFrame(main, text="Filtro de consulta", padding=8)
        filters.pack(fill="x", pady=(0, 8))

        ttk.Label(filters, text="Cromosoma:").grid(
            row=0, column=0, padx=4, pady=4, sticky="w"
        )

        self.chrom_var = tk.StringVar(value="all")
        self.chrom_combo = ttk.Combobox(
            filters, textvariable=self.chrom_var,
            values=["all"] + CHROM_ORDER, state="readonly", width=10
        )
        self.chrom_combo.grid(row=0, column=1, padx=4, pady=4)

        ttk.Label(filters, text="Posición inicial:").grid(
            row=0, column=2, padx=4, pady=4, sticky="w"
        )
        self.start_var = tk.StringVar(value="")
        ttk.Entry(
            filters, textvariable=self.start_var, width=15
        ).grid(row=0, column=3, padx=4, pady=4)

        ttk.Label(filters, text="Posición final:").grid(
            row=0, column=4, padx=4, pady=4, sticky="w"
        )
        self.end_var = tk.StringVar(value="")
        ttk.Entry(
            filters, textvariable=self.end_var, width=15
        ).grid(row=0, column=5, padx=4, pady=4)

        self.scan_button = ttk.Button(
            filters, text="Cargar / filtrar VCF", command=self.start_scan
        )
        self.scan_button.grid(row=0, column=6, padx=5, pady=4)

        self.gene_button = ttk.Button(
            filters, text="Ir a BUB1B", command=self.go_to_bub1b
        )
        self.gene_button.grid(row=0, column=7, padx=5, pady=4)

        self.validate_button = ttk.Button(
            filters, text="Validar REF contra FASTA",
            command=self.start_ref_validation
        )
        self.validate_button.grid(row=0, column=8, padx=5, pady=4)

        self.candidates_button = ttk.Button(
            filters, text="Ver genes candidatos",
            command=self.show_candidate_table
        )
        self.candidates_button.grid(
            row=1, column=0, columnspan=3,
            padx=4, pady=5, sticky="w"
        )

        self.status_var = tk.StringVar(value="Preparando la aplicación...")
        ttk.Label(
            main, textvariable=self.status_var, foreground="#1d5f91"
        ).pack(anchor="w", pady=(0, 5))

        summary = ttk.LabelFrame(main, text="Resumen", padding=8)
        summary.pack(fill="x", pady=(0, 8))

        self.summary_var = tk.StringVar(value="Esperando la carga de datos.")
        ttk.Label(
            summary, textvariable=self.summary_var, justify="left"
        ).pack(anchor="w")

        plot_frame = ttk.LabelFrame(
            main, text="Distribución de variantes (muestra visual)", padding=5
        )
        plot_frame.pack(fill="both", expand=True, pady=(0, 8))

        self.canvas = tk.Canvas(
            plot_frame, background=PLOT_BACKGROUND,
            highlightthickness=0, height=330
        )
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", lambda _event: self.draw_plot())

        legend_frame = ttk.Frame(plot_frame)
        legend_frame.pack(fill="x", pady=(5, 0))

        for gene_name, color in GENE_COLORS.items():
            item = ttk.Frame(legend_frame)
            item.pack(side="left", padx=(3, 12))

            swatch = tk.Canvas(
                item, width=13, height=13, highlightthickness=0,
                background=PLOT_BACKGROUND
            )
            swatch.pack(side="left", padx=(0, 4))
            swatch.create_rectangle(1, 1, 12, 12, fill=color, outline=color)
            ttk.Label(item, text=gene_name).pack(side="left")

        other = ttk.Frame(legend_frame)
        other.pack(side="left", padx=(3, 12))
        swatch = tk.Canvas(
            other, width=13, height=13, highlightthickness=0,
            background=PLOT_BACKGROUND
        )
        swatch.pack(side="left", padx=(0, 4))
        swatch.create_rectangle(
            1, 1, 12, 12,
            fill=OTHER_VARIANT_COLOR, outline=OTHER_VARIANT_COLOR
        )
        ttk.Label(other, text="Otras variantes").pack(side="left")

        bottom = ttk.Panedwindow(main, orient="horizontal")
        bottom.pack(fill="both", expand=True)

        hpo_frame = ttk.LabelFrame(bottom, text="Fenotipos HPO", padding=5)
        log_frame = ttk.LabelFrame(bottom, text="Registro de ejecución", padding=5)
        bottom.add(hpo_frame, weight=1)
        bottom.add(log_frame, weight=2)

        self.hpo_text = tk.Text(
            hpo_frame, height=8, width=42, wrap="word", state="disabled"
        )
        self.hpo_text.pack(fill="both", expand=True)

        self.log_text = tk.Text(
            log_frame, height=8, wrap="word", state="disabled"
        )
        self.log_text.pack(fill="both", expand=True)

    def _check_paths(self):
        self.log_message(f"Directorio de datos: {BASE_DIR}")
        for label, path in [
            ("JSON clínico", CLINICAL_PATH),
            ("GENCODE GTF", GTF_PATH),
            ("FASTA GRCh38", FASTA_PATH),
            ("VCF", VCF_PATH),
            ("Índice VCF", VCF_INDEX_PATH),
        ]:
            self.log_message(
                f"{'OK' if os.path.isfile(path) else 'NO ENCONTRADO'}"
                f" | {label}: {path}"
            )

        self.log_message(
            "El VCF se lee secuencialmente; el índice .tbi no se utiliza."
        )

    def log_message(self, message):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", str(message) + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def set_hpo_text(self, entries):
        self.hpo_text.configure(state="normal")
        self.hpo_text.delete("1.0", "end")

        if not entries:
            self.hpo_text.insert(
                "end", "No se detectaron entradas HPO en el JSON."
            )
        else:
            for item in entries:
                self.hpo_text.insert(
                    "end", f"{item['hpo_id']} | {item['label']}\n"
                )

        self.hpo_text.configure(state="disabled")

    def _load_metadata_worker(self):
        try:
            clinical, hpo, gene, errors = None, [], None, []

            if os.path.isfile(CLINICAL_PATH):
                try:
                    clinical, hpo = load_clinical_json(CLINICAL_PATH)
                except Exception as exc:
                    errors.append(f"JSON: {type(exc).__name__}: {exc}")
            else:
                errors.append(f"No se encuentra JSON: {CLINICAL_PATH}")

            if os.path.isfile(GTF_PATH):
                try:
                    gene = load_gene_from_gtf(GTF_PATH, TARGET_GENE)
                except Exception as exc:
                    errors.append(f"GTF: {type(exc).__name__}: {exc}")
            else:
                errors.append(f"No se encuentra GTF: {GTF_PATH}")

            self.result_queue.put((
                "metadata_ok",
                {"clinical": clinical, "hpo": hpo,
                 "gene": gene, "errors": errors}
            ))
        except Exception as exc:
            self.result_queue.put((
                "metadata_error", f"{type(exc).__name__}: {exc}"
            ))

    def _parse_optional_position(self, value, field_name):
        value = value.strip()
        if not value:
            return None
        try:
            number = int(value)
        except ValueError:
            raise ValueError(f"{field_name} debe ser entero o quedar vacío.")
        if number < 1:
            raise ValueError(f"{field_name} debe ser >= 1.")
        return number

    def start_scan(self):
        if self.busy:
            messagebox.showinfo("Proceso en curso", "Espera a que termine.")
            return

        try:
            chrom = normalize_chrom(self.chrom_var.get())
            start = self._parse_optional_position(
                self.start_var.get(), "Posición inicial"
            )
            end = self._parse_optional_position(
                self.end_var.get(), "Posición final"
            )
            if start is not None and end is not None and start > end:
                raise ValueError("La posición inicial supera la final.")
        except ValueError as exc:
            messagebox.showerror("Filtro inválido", str(exc))
            return

        if not os.path.isfile(VCF_PATH):
            messagebox.showerror("Falta el VCF", VCF_PATH)
            return

        self.busy = True
        self.scan_button.configure(state="disabled")
        self.validate_button.configure(state="disabled")
        self.status_var.set("Escaneando el VCF en segundo plano...")
        self.log_message(
            f"Consulta: cromosoma={chrom}, inicio={start}, final={end}"
        )

        def worker():
            try:
                result = scan_vcf(
                    VCF_PATH, chrom, start, end,
                    progress_callback=lambda msg:
                    self.result_queue.put(("progress", msg))
                )
                self.result_queue.put(("scan_ok", result))
            except Exception as exc:
                self.result_queue.put((
                    "scan_error", f"{type(exc).__name__}: {exc}"
                ))

        threading.Thread(target=worker, daemon=True).start()

    def go_to_bub1b(self):
        gene = self.bub1b or FALLBACK_GENE
        self.chrom_var.set(normalize_chrom(gene["chrom"]))
        self.start_var.set(str(gene["start"]))
        self.end_var.set(str(gene["end"]))
        self.log_message(
            f"Filtro BUB1B: chr{gene['chrom']}:{gene['start']}-{gene['end']}"
        )
        self.start_scan()

    # --------------------------------------------------------
    # TABLA DE CANDIDATOS
    # --------------------------------------------------------

    def show_candidate_table(self):
        if (
            self.candidate_window is not None
            and self.candidate_window.winfo_exists()
        ):
            self.candidate_window.lift()
            self._refresh_candidate_table()
            return

        window = tk.Toplevel(self.root)
        self.candidate_window = window
        window.title("TRINOVA — Variantes candidatas")
        window.geometry("1250x490")
        window.minsize(900, 350)

        ttk.Label(
            window, text="Variantes candidatas configuradas",
            font=("Segoe UI", 12, "bold")
        ).pack(anchor="w", padx=10, pady=(10, 3))

        ttk.Label(
            window,
            text=(
                "La presencia en el VCF no confirma patogenicidad. "
                "EPCR se presenta como dato de entrada, no como "
                "probabilidad de patogenicidad."
            ),
            wraplength=1200
        ).pack(anchor="w", padx=10, pady=(0, 8))

        columns = (
            "gene", "chrom", "pos", "ref", "alt",
            "gt", "dp", "gq", "epcr", "type", "status"
        )

        container = ttk.Frame(window, padding=(10, 0, 10, 8))
        container.pack(fill="both", expand=True)

        tree = ttk.Treeview(
            container, columns=columns, show="headings", height=14
        )
        self.candidate_tree = tree

        headings = {
            "gene": "Gen", "chrom": "Cromosoma", "pos": "Posición",
            "ref": "REF", "alt": "ALT", "gt": "GT", "dp": "DP",
            "gq": "GQ", "epcr": "EPCR", "type": "Tipo",
            "status": "Estado en VCF"
        }
        widths = {
            "gene": 145, "chrom": 85, "pos": 100, "ref": 90,
            "alt": 90, "gt": 65, "dp": 55, "gq": 55,
            "epcr": 65, "type": 90, "status": 240
        }

        for col in columns:
            tree.heading(col, text=headings[col])
            tree.column(
                col, width=widths[col], minwidth=45,
                anchor="center", stretch=True
            )

        vertical = ttk.Scrollbar(
            container, orient="vertical", command=tree.yview
        )
        horizontal = ttk.Scrollbar(
            container, orient="horizontal", command=tree.xview
        )
        tree.configure(
            yscrollcommand=vertical.set, xscrollcommand=horizontal.set
        )
        tree.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal.grid(row=1, column=0, sticky="ew")
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)

        for gene_name, color in GENE_COLORS.items():
            tree.tag_configure(gene_name, foreground=color)

        self._refresh_candidate_table()

        ttk.Button(
            window, text="Actualizar tabla",
            command=self._refresh_candidate_table
        ).pack(anchor="e", padx=10, pady=(0, 10))

    def _refresh_candidate_table(self):
        tree = self.candidate_tree
        if tree is None or not tree.winfo_exists():
            return

        for item_id in tree.get_children():
            tree.delete(item_id)

        for candidate in CANDIDATE_VARIANTS:
            record = self.candidate_matches.get(candidate_key(candidate))

            if record is not None:
                status = "Encontrada en VCF"
                gt = record.get("gt", "N/D")
                dp = format_number(record.get("dp"))
                gq = format_number(record.get("gq"))
            elif self.scan_completed:
                status = "No encontrada en VCF"
                gt = dp = gq = "N/D"
            else:
                status = "Pendiente de escaneo"
                gt = dp = gq = "N/D"

            tree.insert(
                "", "end",
                values=(
                    candidate["gene"],
                    "chr" + normalize_chrom(candidate["chrom"]),
                    candidate["pos"], candidate["ref"], candidate["alt"],
                    gt, dp, gq, f"{candidate['epcr']:.2f}",
                    candidate["finding_type"], status
                ),
                tags=(candidate["gene"],)
            )

    # --------------------------------------------------------
    # VALIDACIÓN REF
    # --------------------------------------------------------

    def start_ref_validation(self):
        if self.busy:
            messagebox.showinfo("Proceso en curso", "Espera a que termine.")
            return
        if not os.path.isfile(VCF_PATH):
            messagebox.showerror("Falta el VCF", VCF_PATH)
            return
        if not os.path.isfile(FASTA_PATH):
            messagebox.showerror("Falta el FASTA", FASTA_PATH)
            return

        gene = (self.bub1b or FALLBACK_GENE).copy()
        self.busy = True
        self.scan_button.configure(state="disabled")
        self.validate_button.configure(state="disabled")
        self.status_var.set("Validando REF de BUB1B contra FASTA...")
        self.log_message("Inicio de validación REF de BUB1B.")

        def worker():
            try:
                result = validate_bub1b_ref(
                    VCF_PATH, FASTA_PATH, gene,
                    progress_callback=lambda msg:
                    self.result_queue.put(("progress", msg))
                )
                self.result_queue.put(("ref_ok", result))
            except Exception as exc:
                self.result_queue.put((
                    "ref_error", f"{type(exc).__name__}: {exc}"
                ))

        threading.Thread(target=worker, daemon=True).start()

    # --------------------------------------------------------
    # FILTRO ACTIVO Y CANDIDATOS VISIBLES
    # --------------------------------------------------------

    def _record_is_in_current_filter(self, record):
        if not self.results:
            return False

        target = self.results["target_chrom"]
        start = self.results["start_pos"]
        end = self.results["end_pos"]

        if target != "all" and record["chrom"] != target:
            return False
        if start is not None and record["pos"] < start:
            return False
        if end is not None and record["pos"] > end:
            return False
        return True

    def _records_for_plot(self):
        """
        Devuelve la muestra visual más TODOS los candidatos encontrados
        que pertenezcan al filtro actual. Evita duplicados.
        """
        if not self.results:
            return []

        records = list(self.results["records"])
        seen = {
            (
                r["chrom"], r["pos"], r["ref"].upper(),
                r["alt"].upper()
            )
            for r in records
        }

        for candidate in CANDIDATE_VARIANTS:
            record = self.candidate_matches.get(candidate_key(candidate))
            if record is None:
                continue
            if not self._record_is_in_current_filter(record):
                continue

            # Un registro multialélico puede contener varios ALT.
            # Guardamos una copia con el ALT específico del candidato
            # para evitar que las distintas variantes se confundan.
            plot_record = record.copy()
            plot_record["alt"] = candidate["alt"].upper()

            identity = (
                plot_record["chrom"], plot_record["pos"],
                plot_record["ref"].upper(), plot_record["alt"].upper()
            )
            if identity not in seen:
                records.append(plot_record)
                seen.add(identity)

        records.sort(key=lambda r: (
            CHROM_ORDER.index(r["chrom"])
            if r["chrom"] in CHROM_ORDER else 99,
            r["pos"]
        ))
        return records

    # --------------------------------------------------------
    # RESUMEN
    # --------------------------------------------------------

    def update_summary(self):
        lines = []

        if self.results:
            r = self.results
            lines.extend([
                f"Registros VCF recorridos: {r['total_records']:,}",
                f"Registros que cumplen el filtro: {r['filtered_records']:,}",
                f"Puntos de la muestra: {r['plotted_records']:,}",
                f"Filtro: cromosoma={r['target_chrom']}, "
                f"inicio={r['start_pos']}, final={r['end_pos']}",
            ])

        if self.bub1b:
            g = self.bub1b
            length = g["end"] - g["start"] + 1
            source = "GENCODE" if self.gtf_gene_found else "coordenadas de respaldo"
            lines.append(
                f"BUB1B: chr{g['chrom']}:{g['start']}-{g['end']} "
                f"({g['strand']}); {length:,} bp; fuente={source}"
            )

        if self.scan_completed:
            lines.append(
                f"Candidatos encontrados en el VCF: "
                f"{len(self.candidate_matches)}/{len(CANDIDATE_VARIANTS)}"
            )

        if self.ref_results:
            counts = self.ref_results["counts"]
            lines.append(
                "Validación REF BUB1B: "
                f"coinciden={counts.get('COINCIDE', 0)}, "
                f"no coinciden={counts.get('NO COINCIDE', 0)}, "
                f"no evaluables={counts.get('NO EVALUABLE', 0)}"
            )

        self.summary_var.set(
            "\n".join(lines) if lines else "Esperando la carga de datos."
        )

    # --------------------------------------------------------
    # COORDENADAS DEL GRÁFICO
    # --------------------------------------------------------

    def _plot_coordinates(self, records, width, height):
        left, right = 65, max(70, width - 20)
        top, bottom = 20, max(50, height - 45)
        plot_w = max(1, right - left)
        plot_h = max(1, bottom - top)

        selected_chrom = self.results["target_chrom"]
        values = []

        for record in records:
            if selected_chrom == "all":
                chrom = record["chrom"]
                try:
                    chrom_idx = CHROM_ORDER.index(chrom)
                except ValueError:
                    chrom_idx = 0

                length = CHROM_LENGTHS.get(chrom, 1)
                fraction = min(1.0, max(0.0, record["pos"] / length))
                x_value = chrom_idx + fraction
            else:
                x_value = record["pos"]

            y_value = record["qual"]
            if y_value is None:
                y_value = record["dp"]
            if y_value is None:
                y_value = 0.0

            values.append((x_value, float(y_value), record))

        if not values:
            return [], (left, right, top, bottom)

        xs = [item[0] for item in values]
        ys = [item[1] for item in values]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)

        if xmin == xmax:
            xmax = xmin + 1
        if ymin == ymax:
            ymax = ymin + 1

        points = []
        for x_value, y_value, record in values:
            x = left + (x_value - xmin) / (xmax - xmin) * plot_w
            y = bottom - (y_value - ymin) / (ymax - ymin) * plot_h
            points.append((x, y, record))

        return points, (left, right, top, bottom)

    # --------------------------------------------------------
    # DIBUJO DEL GRÁFICO
    # --------------------------------------------------------

    def draw_plot(self):
        if not hasattr(self, "canvas"):
            return

        canvas = self.canvas
        canvas.delete("all")
        width = max(100, canvas.winfo_width())
        height = max(100, canvas.winfo_height())

        canvas.create_text(
            12, 8, anchor="nw",
            text="QUAL (si existe) o DP; sin dato = línea base",
            fill=PLOT_TEXT_COLOR, font=("Segoe UI", 9)
        )

        if not self.results:
            canvas.create_text(
                width // 2, height // 2,
                text="Carga el VCF para visualizar variantes",
                fill="#e9eff8", font=("Segoe UI", 12)
            )
            return

        # CLAVE: se añaden candidatos encontrados aunque el muestreo
        # aleatorio de 12.000 puntos no los haya seleccionado.
        records = self._records_for_plot()

        if not records:
            canvas.create_text(
                width // 2, height // 2,
                text="No hay variantes dentro del filtro actual",
                fill="#e9eff8", font=("Segoe UI", 12)
            )
            return

        points, bounds = self._plot_coordinates(records, width, height)
        left, right, top, bottom = bounds

        canvas.create_line(left, top, left, bottom, fill="#9aa9bd")
        canvas.create_line(left, bottom, right, bottom, fill="#9aa9bd")

        # Dibujar primero las variantes generales y después los candidatos,
        # para que los puntos de colores queden por encima.
        for x, y, record in points:
            if candidate_gene_for_record(record) is not None:
                continue
            canvas.create_oval(
                x - 2, y - 2, x + 2, y + 2,
                outline="", fill=OTHER_VARIANT_COLOR
            )

        for x, y, record in points:
            gene_name = candidate_gene_for_record(record)
            if gene_name is None:
                continue

            color = GENE_COLORS[gene_name]
            radius = 5

            canvas.create_oval(
                x - radius, y - radius, x + radius, y + radius,
                outline="#ffffff", width=1, fill=color
            )

            canvas.create_text(
                x + 7, y - 8, text=gene_name,
                anchor="w", fill=color,
                font=("Segoe UI", 8, "bold")
            )

        if self.results["target_chrom"] == "all":
            canvas.create_text(
                (left + right) / 2, height - 18,
                text="Cromosomas (agrupados por cromosoma)",
                fill=PLOT_TEXT_COLOR, font=("Segoe UI", 9)
            )

            span = right - left
            n = len(CHROM_ORDER)
            for i, chrom in enumerate(CHROM_ORDER):
                x = left + (i + 0.5) / n * span
                canvas.create_text(
                    x, bottom + 12, text=chrom,
                    fill=PLOT_TEXT_COLOR, font=("Segoe UI", 7)
                )
        else:
            canvas.create_text(
                (left + right) / 2, height - 18,
                text=f"Posición genómica — cromosoma "
                     f"{self.results['target_chrom']}",
                fill=PLOT_TEXT_COLOR, font=("Segoe UI", 9)
            )

        canvas.create_text(
            8, (top + bottom) / 2,
            text="QUAL / DP", fill=PLOT_TEXT_COLOR,
            angle=90, font=("Segoe UI", 9)
        )

    # --------------------------------------------------------
    # RESULTADOS DE LOS HILOS
    # --------------------------------------------------------

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.result_queue.get_nowait()

                if kind == "metadata_ok":
                    self.clinical_data = payload["clinical"]
                    self.hpo_entries = payload["hpo"]

                    if payload["gene"]:
                        self.bub1b = payload["gene"]
                        self.gtf_gene_found = True
                    else:
                        self.bub1b = FALLBACK_GENE.copy()
                        self.gtf_gene_found = False
                        self.log_message(
                            "Aviso: BUB1B no está en el GTF; "
                            "se usan coordenadas de respaldo."
                        )

                    self.set_hpo_text(self.hpo_entries)
                    self.log_message(
                        f"Fenotipos HPO detectados: {len(self.hpo_entries)}"
                    )

                    g = self.bub1b
                    self.log_message(
                        f"BUB1B: chr{g['chrom']}:{g['start']}-"
                        f"{g['end']} ({g['strand']})"
                    )

                    for error in payload["errors"]:
                        self.log_message("AVISO: " + error)

                    self.update_summary()
                    self.status_var.set("Metadatos cargados.")

                elif kind == "metadata_error":
                    self.log_message(f"Error de metadatos: {payload}")
                    self.status_var.set("Error al cargar metadatos.")

                elif kind == "progress":
                    self.status_var.set(payload)

                elif kind == "scan_ok":
                    self.results = payload
                    self.candidate_matches = payload["candidate_matches"]
                    self.scan_completed = True
                    self.busy = False

                    self.scan_button.configure(state="normal")
                    self.validate_button.configure(state="normal")

                    self.log_message("Escaneo VCF completado.")
                    self.log_message(
                        f"Registros recorridos: {payload['total_records']:,}"
                    )
                    self.log_message(
                        f"Registros filtrados: {payload['filtered_records']:,}"
                    )
                    self.log_message(
                        f"Puntos de muestra: {payload['plotted_records']:,}"
                    )
                    self.log_message(
                        f"Candidatos exactos encontrados: "
                        f"{len(self.candidate_matches)}/"
                        f"{len(CANDIDATE_VARIANTS)}"
                    )

                    for candidate in CANDIDATE_VARIANTS:
                        key = candidate_key(candidate)
                        if key in self.candidate_matches:
                            self.log_message(
                                f"CANDIDATO ENCONTRADO: {candidate['gene']} "
                                f"chr{candidate['chrom']}:{candidate['pos']} "
                                f"{candidate['ref']}>{candidate['alt']}"
                            )
                        else:
                            self.log_message(
                                f"CANDIDATO NO ENCONTRADO: "
                                f"{candidate['gene']} "
                                f"chr{candidate['chrom']}:{candidate['pos']} "
                                f"{candidate['ref']}>{candidate['alt']}"
                            )

                    counts = payload["chromosome_counts"]
                    if counts:
                        summary = ", ".join(
                            f"{chrom}: {count:,}"
                            for chrom, count in sorted(
                                counts.items(),
                                key=lambda item: (
                                    CHROM_ORDER.index(item[0])
                                    if item[0] in CHROM_ORDER else 99
                                )
                            )
                        )
                        self.log_message("Conteo por cromosoma: " + summary)

                    self.update_summary()
                    self.draw_plot()
                    self._refresh_candidate_table()
                    self.status_var.set("Consulta completada.")

                elif kind == "scan_error":
                    self.busy = False
                    self.scan_button.configure(state="normal")
                    self.validate_button.configure(state="normal")
                    self.log_message("ERROR VCF: " + payload)
                    self.status_var.set("Error durante el escaneo.")
                    messagebox.showerror("Error VCF", payload)

                elif kind == "ref_ok":
                    self.ref_results = payload
                    self.busy = False
                    self.scan_button.configure(state="normal")
                    self.validate_button.configure(state="normal")

                    counts = payload["counts"]
                    self.log_message("Validación REF BUB1B terminada.")
                    self.log_message(
                        f"Variantes examinadas: {len(payload['variants']):,}"
                    )
                    self.log_message(
                        f"REF coincide: {counts.get('COINCIDE', 0):,}"
                    )
                    self.log_message(
                        f"REF no coincide: {counts.get('NO COINCIDE', 0):,}"
                    )
                    self.log_message(
                        f"No evaluables: {counts.get('NO EVALUABLE', 0):,}"
                    )

                    for item in payload["variants"]:
                        self.log_message(
                            f"chr{item['chrom']}:{item['pos']} "
                            f"{item['ref']}>{item['alt']} | "
                            f"FASTA={item['expected_ref']} | "
                            f"{item['status']} | GT={item['gt']} | "
                            f"DP={item['dp']} | GQ={item['gq']}"
                        )

                    self.update_summary()
                    self.status_var.set("Validación REF completada.")

                elif kind == "ref_error":
                    self.busy = False
                    self.scan_button.configure(state="normal")
                    self.validate_button.configure(state="normal")
                    self.log_message("ERROR FASTA/VCF: " + payload)
                    self.status_var.set("Error en validación REF.")
                    messagebox.showerror("Error de validación", payload)

        except queue.Empty:
            pass

        self.root.after(150, self._poll_queue)


# ============================================================
# 8. INICIO
# ============================================================

def main():
    root = tk.Tk()
    TriNovaVisualizer(root)
    root.mainloop()


if __name__ == "__main__":
    main()
