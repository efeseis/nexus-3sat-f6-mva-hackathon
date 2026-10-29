# -*- coding: utf-8 -*-

"""
===============================================================================
TRINOVA.OS
ZOOM FUNCIONAL BUB1B
===============================================================================

OBJETIVO
--------
Segunda capa sobre el resultado de ZOOM PROFUNDO BUB1B.

NO modifica el motor TriNova original.

Este script:
    1. Carga VCF, GENCODE v50, FASTA GRCh38 y clinico.json.
    2. Localiza BUB1B usando la anotación GENCODE v50.
    3. Recupera las variantes dentro de BUB1B.
    4. Separa SNV / INDEL.
    5. Valida REF contra FASTA GRCh38.
    6. Determina solapamiento:
           - CDS
           - UTR
           - exon no-CDS
           - intron
           - gene span
    7. Para SNV dentro de CDS intenta calcular:
           - transcript
           - codón
           - posición del codón
           - codón REF
           - codón ALT
           - aminoácido REF
           - aminoácido ALT
           - cambio proteico
    8. Marca variantes cercanas a límites exon/intron.
    9. Conserva GT / AD / DP / GQ / QUAL / FILTER / AF.
   10. Conserva la representación ternaria TriNova para SNV.
   11. Genera archivos TSV y un resumen legible.

IMPORTANTE
----------
No se asigna:
    - patogenicidad
    - benignidad
    - causalidad
    - ACMG
    - score clínico

El objetivo es caracterizar molecularmente las variantes encontradas.
===============================================================================
"""

import gzip
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict


# =============================================================================
# RUTAS
# =============================================================================

BASE = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data"
)

VCF_PATH = os.path.join(
    BASE,
    "WGS_proband01.vcf.gz"
)

VCF_INDEX = os.path.join(
    BASE,
    "WGS_proband01.vcf.gz.tbi"
)

GTF_PATH = os.path.join(
    BASE,
    "gencode.v50.annotation.gtf.gz"
)

FASTA_PATH = os.path.join(
    BASE,
    "GRCh38.primary_assembly.genome.fa.gz"
)

CLINICAL_PATH = os.path.join(
    BASE,
    "clinico.json"
)


# =============================================================================
# CONFIGURACIÓN BUB1B
# =============================================================================

BUB1B_GENE_NAME = "BUB1B"
BUB1B_GENE_ID = "ENSG00000156970.15"

BUB1B_CHROM_GTF = "chr15"

BUB1B_START = 40160984
BUB1B_END = 40221137

# Ventana para detectar proximidad a splice junction.
# No significa que una variante sea spliceogenic.
SPLICE_WINDOW = 10


# =============================================================================
# UTILIDADES
# =============================================================================
# =============================================================================
# NORMALIZACIÓN DE CROMOSOMAS
# =============================================================================

def normalize_chrom(chrom):
    """
    Normaliza nombres cromosómicos para poder comparar VCF,
    GENCODE y FASTA aunque utilicen convenciones distintas.

    Ejemplos:

        chr15 -> 15
        15    -> 15

        chrX  -> X
        X     -> X

        chrY  -> Y
        Y     -> Y

        chrM  -> MT
        M     -> MT
        MT    -> MT
        chrMT -> MT
    """

    if chrom is None:
        return ""

    chrom = str(chrom).strip()

    if not chrom:
        return ""

    # Eliminar prefijo chr
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    chrom = chrom.upper()

    # Unificar cromosoma mitocondrial
    if chrom in ("M", "MT"):
        chrom = "MT"

    return chrom


def open_text(path):
    """
    Abre archivos .gz como texto o archivos normales como texto.
    """
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def safe_int(value, default=None):
    try:
        return int(value)
    except Exception:
        return default


def safe_float(value, default=None):
    try:
        return float(value)
    except Exception:
        return default


def reverse_complement(seq):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(table)[::-1]


def nucleotide_to_trit(base):
    """
    MISMO MAPEO DEL MOTOR TRINOVA ORIGINAL:

        A -> +1
        T/U -> -1
        C/G -> 0
    """
    base = base.upper()

    if base == "A":
        return 1

    if base in ("T", "U"):
        return -1

    if base in ("C", "G"):
        return 0

    return None


def shannon_entropy(values):
    """
    Shannon base 2.
    """
    if not values:
        return 0.0

    counts = Counter(values)
    total = len(values)

    entropy = 0.0

    for count in counts.values():
        p = count / total

        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def parse_gt(gt):
    """
    Acepta:
        0/0
        0/1
        1/1
        1|0
        ./.
    """
    if not gt:
        return []

    gt = gt.replace("|", "/")

    parts = gt.split("/")

    result = []

    for part in parts:
        if part == ".":
            result.append(None)
        else:
            try:
                result.append(int(part))
            except Exception:
                result.append(None)

    return result


def genotype_observed(gt):
    alleles = parse_gt(gt)

    return any(
        allele is not None and allele > 0
        for allele in alleles
    )


# =============================================================================
# GTF
# =============================================================================

def parse_gtf_attributes(text):
    attrs = {}

    pattern = re.compile(
        r'(\S+)\s+"([^"]*)"'
    )

    for key, value in pattern.findall(text):
        attrs[key] = value

    return attrs


def interval_contains(start, end, pos):
    return start <= pos <= end


def interval_overlap(start1, end1, start2, end2):
    return not (
        end1 < start2 or
        end2 < start1
    )


def load_bub1b_annotation():

    print()
    print("=" * 80)
    print("CARGANDO GENCODE v50 — ANOTACIÓN FUNCIONAL BUB1B")
    print("=" * 80)

    transcripts = {}

    exon_intervals = []
    cds_intervals = []
    utr_intervals = []

    # transcript_id -> listas
    transcript_exons = defaultdict(list)
    transcript_cds = defaultdict(list)
    transcript_utrs = defaultdict(list)

    gene_records = 0
    annotation_records = 0

    with open_text(GTF_PATH) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chrom = fields[0]
            source = fields[1]
            feature = fields[2]

            start = safe_int(fields[3])
            end = safe_int(fields[4])

            score = fields[5]
            strand = fields[6]
            frame = fields[7]

            attrs = parse_gtf_attributes(fields[8])

            gene_name = attrs.get("gene_name", "")
            gene_id = attrs.get("gene_id", "")

            if gene_name != BUB1B_GENE_NAME:
                continue

            if gene_id != BUB1B_GENE_ID:
                continue

            annotation_records += 1

            transcript_id = attrs.get("transcript_id")

            if feature == "gene":
                gene_records += 1

            if transcript_id:

                transcripts[transcript_id] = {
                    "transcript_id": transcript_id,
                    "transcript_name": attrs.get(
                        "transcript_name",
                        ""
                    ),
                    "transcript_type": attrs.get(
                        "transcript_type",
                        ""
                    ),
                    "strand": strand,
                    "chrom": chrom,
                }

            record = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "frame": frame,
                "source": source,
                "feature": feature,
                "attrs": attrs,
                "transcript_id": transcript_id,
            }

            if feature == "exon":

                exon_intervals.append(record)

                if transcript_id:
                    transcript_exons[
                        transcript_id
                    ].append(record)

            elif feature == "CDS":

                cds_intervals.append(record)

                if transcript_id:
                    transcript_cds[
                        transcript_id
                    ].append(record)

            elif feature in (
                "UTR",
                "five_prime_UTR",
                "three_prime_UTR",
            ):

                utr_intervals.append(record)

                if transcript_id:
                    transcript_utrs[
                        transcript_id
                    ].append(record)

    print()
    print("BUB1B")
    print("  Gene ID:", BUB1B_GENE_ID)
    print("  Chrom:", BUB1B_CHROM_GTF)
    print("  Start:", BUB1B_START)
    print("  End:", BUB1B_END)

    print()
    print("Transcripts:", len(transcripts))
    print("Exones:", len(exon_intervals))
    print("CDS:", len(cds_intervals))
    print("UTR:", len(utr_intervals))

    return {
        "transcripts": transcripts,
        "exons": exon_intervals,
        "cds": cds_intervals,
        "utrs": utr_intervals,
        "transcript_exons": transcript_exons,
        "transcript_cds": transcript_cds,
        "transcript_utrs": transcript_utrs,
    }


# =============================================================================
# FASTA
# =============================================================================

class FastaReader:

    def __init__(self, path):
        self.path = path
        self.sequences = {}

    def load_chromosome(self, chromosome):

        if chromosome in self.sequences:
            return self.sequences[chromosome]

        print()
        print("Cargando FASTA:", chromosome)

        wanted = normalize_chrom(chromosome)

        sequence = []
        found = False
        fasta_header = None

        with open_text(self.path) as fh:

            for line in fh:

                line = line.rstrip("\n")

                if line.startswith(">"):

                    header = line[1:].split()[0]

                    header_norm = normalize_chrom(header)

                    if found:
                        break

                    if header_norm == wanted:
                        found = True
                        fasta_header = header

                    continue

                if found:
                    sequence.append(line.strip())

        if not sequence:

            raise RuntimeError(
                f"No se encontró {chromosome} en FASTA"
            )

        seq = "".join(sequence).upper()

        self.sequences[chromosome] = seq

        print(
            f"[OK] FASTA {fasta_header or chromosome} "
            f"(solicitado: {chromosome}): "
            f"{len(seq):,} bases cargadas"
        )

        return seq

    def get_base(self, chromosome, pos):

        seq = self.load_chromosome(chromosome)

        if pos < 1 or pos > len(seq):
            return None

        return seq[pos - 1]

    def get_sequence(self, chromosome, start, end):

        seq = self.load_chromosome(chromosome)

        if start < 1:
            return ""

        if end > len(seq):
            return ""

        return seq[start - 1:end]


# =============================================================================
# CLASIFICACIÓN DE POSICIÓN
# =============================================================================

def find_overlaps(intervals, pos):

    hits = []

    for record in intervals:

        if interval_contains(
            record["start"],
            record["end"],
            pos
        ):
            hits.append(record)

    return hits


def classify_position(pos, annotation):

    cds_hits = find_overlaps(
        annotation["cds"],
        pos
    )

    utr_hits = find_overlaps(
        annotation["utrs"],
        pos
    )

    exon_hits = find_overlaps(
        annotation["exons"],
        pos
    )

    # Importante:
    # CDS tiene prioridad sobre exon.
    if cds_hits:
        region = "CDS"

    elif utr_hits:
        region = "UTR"

    elif exon_hits:
        region = "EXON_NON_CDS"

    else:
        region = "INTRON_OR_GENE_SPAN"

    transcripts = set()

    for hit in (
        cds_hits +
        utr_hits +
        exon_hits
    ):

        tid = hit.get("transcript_id")

        if tid:
            transcripts.add(tid)

    return {
        "region": region,
        "cds_hits": cds_hits,
        "utr_hits": utr_hits,
        "exon_hits": exon_hits,
        "transcripts": sorted(transcripts),
    }


# =============================================================================
# SPLICE JUNCTION
# =============================================================================

def nearest_exon_boundary(pos, annotation):

    nearest_distance = None
    nearest_boundary = None
    nearest_exon = None

    for exon in annotation["exons"]:

        for boundary_name, boundary in (
            ("exon_start", exon["start"]),
            ("exon_end", exon["end"]),
        ):

            distance = abs(pos - boundary)

            if (
                nearest_distance is None or
                distance < nearest_distance
            ):
                nearest_distance = distance
                nearest_boundary = boundary_name
                nearest_exon = exon

    if nearest_distance is None:
        return {
            "splice_distance": None,
            "splice_flag": False,
            "splice_boundary": "",
            "splice_transcript": "",
        }

    tid = nearest_exon.get("transcript_id", "")

    return {
        "splice_distance": nearest_distance,
        "splice_flag": nearest_distance <= SPLICE_WINDOW,
        "splice_boundary": nearest_boundary,
        "splice_transcript": tid or "",
    }


# =============================================================================
# CONSTRUCCIÓN DE SECUENCIA CODIFICANTE
# =============================================================================

def build_transcript_cds(annotation, transcript_id, fasta):

    cds_records = annotation["transcript_cds"].get(
        transcript_id,
        []
    )

    if not cds_records:
        return None

    strand = cds_records[0]["strand"]

    # Orden genómico
    cds_records = sorted(
        cds_records,
        key=lambda x: (
            x["start"],
            x["end"]
        )
    )

    pieces = []

    genomic_positions = []

    for record in cds_records:

        seq = fasta.get_sequence(
            record["chrom"],
            record["start"],
            record["end"]
        )

        if not seq:
            continue

        positions = list(
            range(
                record["start"],
                record["end"] + 1
            )
        )

        if strand == "-":
            seq = reverse_complement(seq)
            positions = positions[::-1]

        pieces.append(seq)
        genomic_positions.extend(positions)

    if not pieces:
        return None

    cds_sequence = "".join(pieces)

    return {
        "transcript_id": transcript_id,
        "strand": strand,
        "sequence": cds_sequence,
        "genomic_positions": genomic_positions,
    }


# =============================================================================
# TABLA CODONES
# =============================================================================

CODON_TABLE = {

    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",

    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",

    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",

    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",

    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",

    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",

    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",

    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",

    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",

    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",

    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",

    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",

    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",

    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",

    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",

    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def protein_from_cds(cds_sequence):

    protein = []

    for i in range(0, len(cds_sequence) - 2, 3):

        codon = cds_sequence[i:i + 3]

        aa = CODON_TABLE.get(
            codon,
            "?"
        )

        protein.append(aa)

    return "".join(protein)


# =============================================================================
# CONSECUENCIA DE SNV EN CDS
# =============================================================================

def calculate_cds_effect(
    chrom,
    pos,
    ref,
    alt,
    transcript_id,
    annotation,
    fasta
):

    cds_info = build_transcript_cds(
        annotation,
        transcript_id,
        fasta
    )

    if not cds_info:
        return None

    positions = cds_info["genomic_positions"]
    sequence = cds_info["sequence"]
    strand = cds_info["strand"]

    if pos not in positions:
        return None

    index = positions.index(pos)

    ref_base_transcript = sequence[index]

    if strand == "+":
        alt_base_transcript = alt.upper()

    else:
        alt_base_transcript = reverse_complement(
            alt.upper()
        )

    if len(ref) != 1 or len(alt) != 1:
        return None

    if len(ref_base_transcript) != 1:
        return None

    codon_start = (index // 3) * 3

    codon_position = (index % 3) + 1

    if codon_start + 3 > len(sequence):
        return None

    ref_codon = sequence[
        codon_start:codon_start + 3
    ]

    alt_codon_list = list(ref_codon)

    alt_codon_list[
        codon_position - 1
    ] = alt_base_transcript

    alt_codon = "".join(
        alt_codon_list
    )

    ref_aa = CODON_TABLE.get(
        ref_codon,
        "?"
    )

    alt_aa = CODON_TABLE.get(
        alt_codon,
        "?"
    )

    if ref_aa == alt_aa:
        consequence = "SYNONYMOUS"

    else:

        if alt_aa == "*":
            consequence = "STOP_GAINED"

        elif ref_aa == "*":
            consequence = "STOP_LOST"

        else:
            consequence = "MISSENSE"

    return {
        "transcript_id": transcript_id,
        "strand": strand,
        "cds_index": index + 1,
        "codon_position": codon_position,
        "ref_base_transcript": ref_base_transcript,
        "alt_base_transcript": alt_base_transcript,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "protein_change": f"{ref_aa}>{alt_aa}",
        "consequence": consequence,
    }


# =============================================================================
# ANOTACIÓN FUNCIONAL DE UNA VARIANTE
# =============================================================================

def annotate_variant(
    variant,
    annotation,
    fasta
):

    chrom = variant["chrom"]
    pos = variant["pos"]
    ref = variant["ref"]
    alt = variant["alt"]

    classification = classify_position(
        pos,
        annotation
    )

    splice = nearest_exon_boundary(
        pos,
        annotation
    )

    variant["REGION"] = classification["region"]

    variant["TRANSCRIPTS"] = ",".join(
        classification["transcripts"]
    )

    variant["SPLICE_DISTANCE"] = (
        splice["splice_distance"]
    )

    variant["NEAR_SPLICE"] = (
        "YES"
        if splice["splice_flag"]
        else "NO"
    )

    variant["SPLICE_BOUNDARY"] = (
        splice["splice_boundary"]
    )

    variant["SPLICE_TRANSCRIPT"] = (
        splice["splice_transcript"]
    )

    # -------------------------------------------------------------------------
    # Validación FASTA
    # -------------------------------------------------------------------------

    fasta_ref = fasta.get_base(
        chrom,
        pos
    )

    variant["FASTA_REF"] = (
        fasta_ref
        if fasta_ref is not None
        else ""
    )

    variant["REF_MATCH_FASTA"] = (
        "YES"
        if fasta_ref is not None
        and fasta_ref.upper() == ref.upper()
        else "NO"
    )

    # -------------------------------------------------------------------------
    # Defaults
    # -------------------------------------------------------------------------

    variant["CDS_TRANSCRIPT"] = ""
    variant["CODON_POSITION"] = ""
    variant["REF_CODON"] = ""
    variant["ALT_CODON"] = ""
    variant["REF_AA"] = ""
    variant["ALT_AA"] = ""
    variant["PROTEIN_CHANGE"] = ""
    variant["CONSEQUENCE"] = ""

    # -------------------------------------------------------------------------
    # Consecuencia CDS
    # -------------------------------------------------------------------------

    if (
        len(ref) == 1
        and len(alt) == 1
        and classification["cds_hits"]
    ):

        for transcript_id in (
            classification["transcripts"]
        ):

            effect = calculate_cds_effect(
                chrom,
                pos,
                ref,
                alt,
                transcript_id,
                annotation,
                fasta
            )

            if effect:

                variant["CDS_TRANSCRIPT"] = (
                    effect["transcript_id"]
                )

                variant["CODON_POSITION"] = (
                    effect["codon_position"]
                )

                variant["REF_CODON"] = (
                    effect["ref_codon"]
                )

                variant["ALT_CODON"] = (
                    effect["alt_codon"]
                )

                variant["REF_AA"] = (
                    effect["ref_aa"]
                )

                variant["ALT_AA"] = (
                    effect["alt_aa"]
                )

                variant["PROTEIN_CHANGE"] = (
                    effect["protein_change"]
                )

                variant["CONSEQUENCE"] = (
                    effect["consequence"]
                )

                # Para este análisis basta con la primera
                # anotación CDS válida.
                break

    return variant


# =============================================================================
# PARSEO DE INFO
# =============================================================================

def parse_info(info):

    result = {}

    if info == ".":
        return result

    for item in info.split(";"):

        if "=" in item:

            key, value = item.split(
                "=",
                1
            )

            result[key] = value

        else:

            result[item] = True

    return result


# =============================================================================
# PARSEO FORMAT / SAMPLE
# =============================================================================

def parse_sample(format_string, sample_string):

    keys = format_string.split(":")
    values = sample_string.split(":")

    result = {}

    for i, key in enumerate(keys):

        if i < len(values):
            result[key] = values[i]

        else:
            result[key] = "."

    return result


# =============================================================================
# ESCANEO VCF
# =============================================================================

def scan_vcf(annotation, fasta):

    print()
    print("=" * 80)
    print("ESCANEANDO VCF — BUB1B")
    print("=" * 80)

    variants = []

    records_scanned = 0

    sample_name = ""

    with open_text(VCF_PATH) as fh:

        for line in fh:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                header = line.rstrip("\n").split("\t")

                if len(header) >= 10:
                    sample_name = header[9]

                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            records_scanned += 1

            chrom = fields[0]
            pos = safe_int(fields[1])

            if pos is None:
                continue

            ref = fields[3]
            alt_string = fields[4]

            qual = fields[5]
            filt = fields[6]
            info_string = fields[7]

            # El VCF puede usar 15 y el GTF chr15 (o viceversa).
            # Se comparan con nombres normalizados para evitar perder variantes.
            if normalize_chrom(chrom) != normalize_chrom(BUB1B_CHROM_GTF):
                continue

            # Solo BUB1B.
            if (
                pos < BUB1B_START
                or pos > BUB1B_END
            ):
                continue

            info = parse_info(
                info_string
            )

            format_string = (
                fields[8]
                if len(fields) >= 9
                else ""
            )

            sample_string = (
                fields[9]
                if len(fields) >= 10
                else ""
            )

            sample = parse_sample(
                format_string,
                sample_string
            )

            gt = sample.get(
                "GT",
                "."
            )

            pgt = sample.get(
                "PGT",
                "."
            )

            pid = sample.get(
                "PID",
                "."
            )

            dp = sample.get(
                "DP",
                "."
            )

            gq = sample.get(
                "GQ",
                "."
            )

            ad = sample.get(
                "AD",
                "."
            )

            af = info.get(
                "AF",
                "."
            )

            ac = info.get(
                "AC",
                "."
            )

            an = info.get(
                "AN",
                "."
            )

            # Multialélico
            alt_alleles = alt_string.split(",")

            for alt_index, alt in enumerate(
                alt_alleles,
                start=1
            ):

                is_snv = (
                    len(ref) == 1
                    and len(alt) == 1
                    and ref.upper() in "ACGT"
                    and alt.upper() in "ACGT"
                )

                if is_snv:

                    trit_ref = nucleotide_to_trit(
                        ref
                    )

                    trit_alt = nucleotide_to_trit(
                        alt
                    )

                    trit_string = (
                        f"{trit_ref:+d}"
                        f"→"
                        f"{trit_alt:+d}"
                    )

                else:

                    trit_string = "NA"

                observed = genotype_observed(
                    gt
                )

                variant = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "type": (
                        "SNV"
                        if is_snv
                        else "INDEL_OTHER"
                    ),
                    "gt": gt,
                    "pgt": pgt,
                    "pid": pid,
                    "ad": ad,
                    "dp": dp,
                    "gq": gq,
                    "qual": qual,
                    "filter": filt,
                    "af": af,
                    "ac": ac,
                    "an": an,
                    "trit": trit_string,
                    "observed": observed,
                }

                variant = annotate_variant(
                    variant,
                    annotation,
                    fasta
                )

                variants.append(
                    variant
                )

    print()
    print("Muestra:", sample_name)
    print("Registros VCF recorridos:", f"{records_scanned:,}")
    print("Variantes BUB1B:", len(variants))

    return variants, sample_name


# =============================================================================
# MÉTRICAS TRINOVA
# =============================================================================

def calculate_metrics(variants):

    snvs = [
        v
        for v in variants
        if v["type"] == "SNV"
    ]

    trits = []

    observed_alt_trits = []

    for v in snvs:

        ref_trit = nucleotide_to_trit(
            v["ref"]
        )

        alt_trit = nucleotide_to_trit(
            v["alt"]
        )

        trits.append(ref_trit)
        trits.append(alt_trit)

        if v["observed"]:

            observed_alt_trits.append(
                alt_trit
            )

    entropy_ref_alt = shannon_entropy(
        trits
    )

    entropy_observed = shannon_entropy(
        observed_alt_trits
    )

    return {
        "snv_count": len(snvs),
        "trit_count": len(trits),
        "observed_alt_count": len(
            observed_alt_trits
        ),
        "entropy_ref_alt": entropy_ref_alt,
        "entropy_observed": entropy_observed,
        "trit_distribution": Counter(
            trits
        ),
        "observed_distribution": Counter(
            observed_alt_trits
        ),
    }


# =============================================================================
# CLÍNICO
# =============================================================================

def load_clinical():

    print()
    print("=" * 80)
    print("CLINICO.JSON")
    print("=" * 80)

    try:

        with open(
            CLINICAL_PATH,
            "r",
            encoding="utf-8"
        ) as fh:

            data = json.load(fh)

        print()
        print("[OK] JSON cargado correctamente.")

        return data

    except Exception as e:

        print()
        print("[ERROR] No se pudo cargar clinico.json")
        print(e)

        return {}


def extract_clinical_features(data):

    features = []

    root = data.get(
        "rare_disease_genomic_challenge",
        data
    )

    if not isinstance(root, dict):
        return features

    clinical = root.get(
        "clinical_features",
        []
    )

    if isinstance(clinical, list):

        for item in clinical:

            if not isinstance(item, dict):
                continue

            features.append({
                "clinical_feature": item.get(
                    "clinical_feature",
                    ""
                ),
                "hpo_term": item.get(
                    "hpo_term",
                    ""
                ),
                "hpo_id": item.get(
                    "hpo_id",
                    ""
                ),
            })

    return features


# =============================================================================
# TABLA
# =============================================================================

OUTPUT_FIELDS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "type",
    "gt",
    "ad",
    "dp",
    "gq",
    "qual",
    "filter",
    "af",
    "ac",
    "an",
    "observed",
    "trit",
    "region",
    "transcripts",
    "fasta_ref",
    "ref_match_fasta",
    "splice_distance",
    "near_splice",
    "splice_boundary",
    "splice_transcript",
    "cds_transcript",
    "codon_position",
    "ref_codon",
    "alt_codon",
    "ref_aa",
    "alt_aa",
    "protein_change",
    "consequence",
]


def variant_to_output_dict(v):

    return {
        "chrom": v["chrom"],
        "pos": v["pos"],
        "ref": v["ref"],
        "alt": v["alt"],
        "type": v["type"],
        "gt": v["gt"],
        "ad": v["ad"],
        "dp": v["dp"],
        "gq": v["gq"],
        "qual": v["qual"],
        "filter": v["filter"],
        "af": v["af"],
        "ac": v["ac"],
        "an": v["an"],
        "observed": v["observed"],
        "trit": v["trit"],
        "region": v["REGION"],
        "transcripts": v["TRANSCRIPTS"],
        "fasta_ref": v["FASTA_REF"],
        "ref_match_fasta": v["REF_MATCH_FASTA"],
        "splice_distance": v[
            "SPLICE_DISTANCE"
        ],
        "near_splice": v[
            "NEAR_SPLICE"
        ],
        "splice_boundary": v[
            "SPLICE_BOUNDARY"
        ],
        "splice_transcript": v[
            "SPLICE_TRANSCRIPT"
        ],
        "cds_transcript": v[
            "CDS_TRANSCRIPT"
        ],
        "codon_position": v[
            "CODON_POSITION"
        ],
        "ref_codon": v[
            "REF_CODON"
        ],
        "alt_codon": v[
            "ALT_CODON"
        ],
        "ref_aa": v[
            "REF_AA"
        ],
        "alt_aa": v[
            "ALT_AA"
        ],
        "protein_change": v[
            "PROTEIN_CHANGE"
        ],
        "consequence": v[
            "CONSEQUENCE"
        ],
    }


def write_tsv(path, variants):

    if not variants:
        return

    rows = [
        variant_to_output_dict(v)
        for v in variants
    ]

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        fh.write(
            "\t".join(OUTPUT_FIELDS)
            + "\n"
        )

        for row in rows:

            values = []

            for field in OUTPUT_FIELDS:

                value = row.get(
                    field,
                    ""
                )

                if value is None:
                    value = ""

                values.append(
                    str(value)
                )

            fh.write(
                "\t".join(values)
                + "\n"
            )


# =============================================================================
# IMPRESIÓN
# =============================================================================

def print_variants(variants):

    print()
    print("=" * 150)
    print("VARIANTES FUNCIONALES BUB1B")
    print("=" * 150)

    header = (
        f"{'POS':>10} "
        f"{'REF':>4} "
        f"{'ALT':>4} "
        f"{'TYPE':<12} "
        f"{'GT':<5} "
        f"{'DP':>4} "
        f"{'GQ':>4} "
        f"{'REGION':<18} "
        f"{'SPLICE':>6} "
        f"{'CODON':<8} "
        f"{'AA':<8} "
        f"{'CONSEQUENCE':<15}"
    )

    print(header)
    print("-" * len(header))

    for v in variants:

        print(
            f"{v['pos']:>10} "
            f"{v['ref']:>4} "
            f"{v['alt']:>4} "
            f"{v['type']:<12} "
            f"{v['gt']:<5} "
            f"{str(v['dp']):>4} "
            f"{str(v['gq']):>4} "
            f"{v['REGION']:<18} "
            f"{str(v['SPLICE_DISTANCE']):>6} "
            f"{(v['REF_CODON'] + '>' + v['ALT_CODON']) if v['REF_CODON'] else '-':<8} "
            f"{(v['REF_AA'] + '>' + v['ALT_AA']) if v['REF_AA'] else '-':<8} "
            f"{v['CONSEQUENCE'] or '-':<15}"
        )


# =============================================================================
# RESUMEN
# =============================================================================

def write_summary(
    path,
    variants,
    metrics,
    sample_name
):

    snvs = [
        v for v in variants
        if v["type"] == "SNV"
    ]

    indels = [
        v for v in variants
        if v["type"] != "SNV"
    ]

    cds = [
        v for v in variants
        if v["REGION"] == "CDS"
    ]

    exons = [
        v for v in variants
        if v["REGION"] in (
            "CDS",
            "EXON_NON_CDS"
        )
    ]

    splice = [
        v for v in variants
        if v["NEAR_SPLICE"] == "YES"
    ]

    consequence_counts = Counter(
        v["CONSEQUENCE"]
        for v in variants
        if v["CONSEQUENCE"]
    )

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "TRINOVA.OS\n"
        )

        fh.write(
            "ZOOM FUNCIONAL BUB1B\n"
        )

        fh.write(
            "=" * 80
            + "\n\n"
        )

        fh.write(
            f"Muestra: {sample_name}\n"
        )

        fh.write(
            f"BUB1B: {BUB1B_CHROM_GTF}:"
            f"{BUB1B_START}-{BUB1B_END}\n\n"
        )

        fh.write(
            "VARIANTES\n"
        )

        fh.write(
            f"Total: {len(variants)}\n"
        )

        fh.write(
            f"SNV: {len(snvs)}\n"
        )

        fh.write(
            f"INDEL/otros: {len(indels)}\n"
        )

        fh.write(
            f"EXON/CDS: {len(exons)}\n"
        )

        fh.write(
            f"CDS: {len(cds)}\n"
        )

        fh.write(
            f"Cerca de splice: {len(splice)}\n\n"
        )

        fh.write(
            "TRINOVA\n"
        )

        fh.write(
            f"Shannon REF/ALT: "
            f"{metrics['entropy_ref_alt']:.9f}\n"
        )

        fh.write(
            f"Shannon observado: "
            f"{metrics['entropy_observed']:.9f}\n"
        )

        fh.write(
            f"Valores ternarios: "
            f"{metrics['trit_count']}\n"
        )

        fh.write("\nDistribución ternaria:\n")

        for trit in (-1, 0, 1):

            count = metrics[
                "trit_distribution"
            ].get(
                trit,
                0
            )

            fh.write(
                f"  {trit:+d}: {count}\n"
            )

        fh.write(
            "\nCONSECUENCIAS CDS\n"
        )

        if consequence_counts:

            for consequence, count in (
                consequence_counts.most_common()
            ):

                fh.write(
                    f"  {consequence}: "
                    f"{count}\n"
                )

        else:

            fh.write(
                "  No calculadas.\n"
            )

        fh.write(
            "\nVARIANTES CDS\n"
        )

        for v in cds:

            fh.write(
                f"  {v['pos']} "
                f"{v['ref']}>{v['alt']} "
                f"GT={v['gt']} "
                f"DP={v['dp']} "
                f"GQ={v['gq']} "
                f"{v['CONSEQUENCE'] or 'NA'} "
                f"{v['PROTEIN_CHANGE'] or 'NA'} "
                f"transcript="
                f"{v['CDS_TRANSCRIPT'] or 'NA'}\n"
            )

        fh.write(
            "\nVARIANTES CERCANAS A SPLICE\n"
        )

        for v in splice:

            fh.write(
                f"  {v['pos']} "
                f"{v['ref']}>{v['alt']} "
                f"dist={v['SPLICE_DISTANCE']} "
                f"GT={v['gt']} "
                f"transcript="
                f"{v['SPLICE_TRANSCRIPT']}\n"
            )


# =============================================================================
# CLINICAL TXT
# =============================================================================

def write_clinical(path, clinical_data):

    features = extract_clinical_features(
        clinical_data
    )

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "TRINOVA.OS\n"
        )

        fh.write(
            "CLINICO — BUB1B\n"
        )

        fh.write(
            "=" * 80
            + "\n\n"
        )

        fh.write(
            "Fenotipo recuperado desde clinico.json\n\n"
        )

        for item in features:

            fh.write(
                f"{item['clinical_feature']}\t"
                f"{item['hpo_term']}\t"
                f"{item['hpo_id']}\n"
            )


# =============================================================================
# ARCHIVOS
# =============================================================================

def output_path(filename):

    # Crear la carpeta "audit 1" junto al script.
    output_dir = os.path.join(
        os.path.dirname(
            os.path.abspath(__file__)
        ),
        "audit 1"
    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    return os.path.join(
        output_dir,
        filename
    )


# =============================================================================
# MAIN
# =============================================================================


def audit_bub1b_track1(variants):

    """
    Auditoría final para preparar la revisión Track 1.

    No asigna patogenicidad ni causalidad. Solo revisa:
      - identidad de las variantes BUB1B recuperadas
      - GT/PGT/PID/DP/GQ/AF
      - concordancia REF vs FASTA ya calculada por annotate_variant()
      - las dos variantes CDS
      - evidencia de fase disponible en PGT/PID
      - candidatos potenciales para compound-het, sin llamarlo causal
    """

    print()
    print("=" * 80)
    print("AUDITORÍA BUB1B — TRACK 1 / FASE Y CANDIDATOS CDS")
    print("=" * 80)

    audit_rows = []

    for v in sorted(
        variants,
        key=lambda x: int(x.get("pos", 0))
    ):

        pos = int(v.get("pos", 0))
        gt = v.get("gt", ".")
        pgt = v.get("pgt", ".")
        pid = v.get("pid", ".")
        region = v.get("REGION", "-")
        consequence = v.get("CONSEQUENCE", "-")
        ref_match = v.get("REF_MATCH_FASTA", "-")

        audit_rows.append({
            "pos": pos,
            "ref": v.get("ref", "."),
            "alt": v.get("alt", "."),
            "type": v.get("type", "."),
            "gt": gt,
            "pgt": pgt,
            "pid": pid,
            "dp": v.get("dp", "."),
            "gq": v.get("gq", "."),
            "af": v.get("af", "."),
            "region": region,
            "consequence": consequence,
            "ref_match_fasta": ref_match,
        })

    print()
    print(
        f"{'POS':>10} {'REF':>6} {'ALT':>6} "
        f"{'GT':>6} {'PGT':>6} {'PID':>12} "
        f"{'DP':>5} {'GQ':>5} {'REGION':>20} "
        f"{'FASTA':>8} {'CONSEQUENCE':>18}"
    )
    print("-" * 145)

    for row in audit_rows:

        print(
            f"{row['pos']:>10} "
            f"{row['ref']:>6} "
            f"{row['alt']:>6} "
            f"{str(row['gt']):>6} "
            f"{str(row['pgt']):>6} "
            f"{str(row['pid']):>12} "
            f"{str(row['dp']):>5} "
            f"{str(row['gq']):>5} "
            f"{str(row['region']):>20} "
            f"{str(row['ref_match_fasta']):>8} "
            f"{str(row['consequence']):>18}"
        )

    mismatches = [
        v for v in variants
        if v.get("REF_MATCH_FASTA") == "NO"
    ]

    cds_variants = [
        v for v in variants
        if v.get("REGION") == "CDS"
    ]

    print()
    print("=" * 80)
    print("1. REF vs FASTA")
    print("=" * 80)

    if not mismatches:

        print("[OK] Todas las variantes coinciden con FASTA.")

    else:

        print(
            f"[ATENCIÓN] {len(mismatches)} variante(s) "
            "no coinciden literalmente con FASTA."
        )

        for v in mismatches:

            print(
                f"  - {v.get('chrom')}:{v.get('pos')} "
                f"{v.get('ref')}>{v.get('alt')} | "
                f"FASTA_REF={v.get('FASTA_REF', '-') } | "
                f"TYPE={v.get('type', '-') }"
            )

    print()
    print("=" * 80)
    print("2. VARIANTES CDS")
    print("=" * 80)

    if not cds_variants:

        print("[ATENCIÓN] No hay variantes anotadas como CDS.")

    else:

        for v in cds_variants:

            print()
            print(
                f"POS {v.get('pos')}: "
                f"{v.get('ref')}>{v.get('alt')}"
            )
            print(
                f"  GT={v.get('gt', '.')} "
                f"PGT={v.get('pgt', '.')} "
                f"PID={v.get('pid', '.')}"
            )
            print(
                f"  DP={v.get('dp', '.')} "
                f"GQ={v.get('gq', '.')} "
                f"AF={v.get('af', '.')}"
            )
            print(
                f"  CODON={v.get('CODON', '-')} "
                f"AA={v.get('AA', '-')} "
                f"CONSEQUENCE={v.get('CONSEQUENCE', '-')}"
            )
            print(
                f"  REF_MATCH_FASTA="
                f"{v.get('REF_MATCH_FASTA', '-')}"
            )

    target_positions = {
        40209701,
        40220612,
    }

    targets = [
        v for v in variants
        if int(v.get("pos", 0)) in target_positions
    ]

    print()
    print("=" * 80)
    print("3. RELACIÓN ENTRE 40209701 Y 40220612")
    print("=" * 80)

    if len(targets) != 2:

        print(
            "[ATENCIÓN] No se recuperaron ambas posiciones objetivo "
            "en el conjunto BUB1B actual."
        )

    else:

        targets = sorted(
            targets,
            key=lambda x: int(x.get("pos", 0))
        )

        for v in targets:

            print(
                f"  {v.get('pos')}: "
                f"{v.get('ref')}>{v.get('alt')} | "
                f"GT={v.get('gt', '.')} | "
                f"PGT={v.get('pgt', '.')} | "
                f"PID={v.get('pid', '.')}"
            )

        heterozygous = all(
            str(v.get("gt", "")).replace("/", "|") in (
                "0|1",
                "1|0",
            )
            for v in targets
        )

        print()

        if heterozygous:
            print(
                "[OK] Ambas variantes son heterocigotas en GT."
            )
        else:
            print(
                "[ATENCIÓN] Las dos variantes no aparecen ambas "
                "como heterocigotas en GT."
            )

        pgt_values = [
            str(v.get("pgt", "."))
            .replace("/", "|")
            for v in targets
        ]

        pid_values = [
            str(v.get("pid", "."))
            for v in targets
        ]

        pgt_informative = all(
            p in ("0|1", "1|0")
            for p in pgt_values
        )

        same_pid = (
            pid_values[0] not in (".", "", "None")
            and pid_values[0] == pid_values[1]
        )

        print()
        print(
            f"PGT informativo en ambas: "
            f"{'SÍ' if pgt_informative else 'NO'}"
        )
        print(
            f"PID compartido: "
            f"{'SÍ' if same_pid else 'NO'}"
        )

        if pgt_informative and same_pid:

            if pgt_values[0] != pgt_values[1]:

                phase_relation = "TRANS_EVIDENCE"

            else:

                phase_relation = "CIS_EVIDENCE"

            print()
            print(
                f"[FASE] {phase_relation}"
            )
            print(
                "  Esta etiqueta describe la información PGT/PID del VCF; "
                "no constituye por sí sola una afirmación de causalidad."
            )

        else:

            print()
            print(
                "[FASE] INDETERMINADA con los campos PGT/PID disponibles."
            )
            print(
                "  No se debe llamar compound-het únicamente por tener "
                "GT 0/1 en ambas variantes."
            )

    print()
    print("=" * 80)
    print("4. INTERPRETACIÓN PARA TRACK 1")
    print("=" * 80)
    print(
        "Este bloque NO asigna patogenicidad, causalidad ni ACMG."
    )
    print(
        "Su objetivo es dejar auditada la evidencia necesaria antes "
        "de construir el CSV final: variante, genotipo, fase y efecto."
    )

    return audit_rows


def write_track1_audit(path, variants):

    rows = audit_bub1b_track1(variants)

    fields = [
        "pos",
        "ref",
        "alt",
        "type",
        "gt",
        "pgt",
        "pid",
        "dp",
        "gq",
        "af",
        "region",
        "consequence",
        "ref_match_fasta",
    ]

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        handle.write(
            "\t".join(fields) + "\n"
        )

        for row in rows:

            handle.write(
                "\t".join(
                    str(row.get(field, "."))
                    for field in fields
                )
                + "\n"
            )

    return path


def main():

    print()
    print("=" * 80)
    print("TRINOVA.OS")
    print("ZOOM FUNCIONAL BUB1B")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # ARCHIVOS
    # -------------------------------------------------------------------------

    print()
    print("COMPROBACIÓN DE ARCHIVOS")
    print()

    paths = [
        ("VCF", VCF_PATH),
        ("VCF INDEX", VCF_INDEX),
        ("GTF", GTF_PATH),
        ("FASTA", FASTA_PATH),
        ("CLINICO", CLINICAL_PATH),
    ]

    for label, path in paths:

        if os.path.exists(path):

            print(
                f"[OK] {label}"
            )

        else:

            print(
                f"[ERROR] {label}: "
                f"{path}"
            )

    missing = [
        label
        for label, path in paths
        if not os.path.exists(path)
    ]

    if missing:

        print()
        print(
            "Faltan archivos:",
            ", ".join(missing)
        )

        input(
            "\nPulsa ENTER para cerrar..."
        )

        return

    # -------------------------------------------------------------------------
    # GENCODE
    # -------------------------------------------------------------------------

    annotation = load_bub1b_annotation()

    # -------------------------------------------------------------------------
    # FASTA
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("PREPARANDO FASTA GRCh38")
    print("=" * 80)

    fasta = FastaReader(
        FASTA_PATH
    )

    # -------------------------------------------------------------------------
    # VCF
    # -------------------------------------------------------------------------

    variants, sample_name = scan_vcf(
        annotation,
        fasta
    )

    if not variants:

        print()
        print(
            "[ERROR] No se encontraron "
            "variantes dentro de BUB1B."
        )

        input(
            "\nPulsa ENTER para cerrar..."
        )

        return

    # -------------------------------------------------------------------------
    # MÉTRICAS
    # -------------------------------------------------------------------------

    metrics = calculate_metrics(
        variants
    )

    # -------------------------------------------------------------------------
    # CLÍNICO
    # -------------------------------------------------------------------------

    clinical_data = load_clinical()

    clinical_features = (
        extract_clinical_features(
            clinical_data
        )
    )

    # -------------------------------------------------------------------------
    # IMPRESIÓN
    # -------------------------------------------------------------------------

    print_variants(
        variants
    )

    # -------------------------------------------------------------------------
    # RESUMEN
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("RESULTADO FUNCIONAL")
    print("=" * 80)

    snvs = [
        v for v in variants
        if v["type"] == "SNV"
    ]

    indels = [
        v for v in variants
        if v["type"] != "SNV"
    ]

    cds = [
        v for v in variants
        if v["REGION"] == "CDS"
    ]

    exons = [
        v for v in variants
        if v["REGION"] in (
            "CDS",
            "EXON_NON_CDS"
        )
    ]

    splice = [
        v for v in variants
        if v["NEAR_SPLICE"] == "YES"
    ]

    fasta_mismatch = [
        v for v in variants
        if v["REF_MATCH_FASTA"] == "NO"
    ]

    print()
    print(
        f"BUB1B: "
        f"{BUB1B_CHROM_GTF} "
        f"{BUB1B_START} - {BUB1B_END}"
    )

    print()
    print(
        f"Total variantes: {len(variants)}"
    )

    print(
        f"SNV: {len(snvs)}"
    )

    print(
        f"INDEL/otros: {len(indels)}"
    )

    print(
        f"EXON/CDS: {len(exons)}"
    )

    print(
        f"CDS: {len(cds)}"
    )

    print(
        f"Cerca de splice: {len(splice)}"
    )

    print(
        f"REF no coincide con FASTA: "
        f"{len(fasta_mismatch)}"
    )

    print()
    print(
        f"Shannon TriNova REF/ALT: "
        f"{metrics['entropy_ref_alt']:.9f}"
    )

    print(
        f"Shannon observado: "
        f"{metrics['entropy_observed']:.9f}"
    )

    print()
    print(
        "Fenotipos clínicos recuperados:",
        len(clinical_features)
    )

    # -------------------------------------------------------------------------
    # ARCHIVOS
    # -------------------------------------------------------------------------

    all_path = output_path(
        "bub1b_funcional_todas.tsv"
    )

    snv_path = output_path(
        "bub1b_funcional_snv.tsv"
    )

    cds_path = output_path(
        "bub1b_funcional_cds.tsv"
    )

    splice_path = output_path(
        "bub1b_funcional_splice.tsv"
    )

    summary_path = output_path(
        "bub1b_funcional_summary.txt"
    )

    clinical_path = output_path(
        "bub1b_funcional_clinical.txt"
    )

    audit_path = output_path(
        "bub1b_track1_audit.tsv"
    )

    write_tsv(
        all_path,
        variants
    )

    write_tsv(
        snv_path,
        snvs
    )

    write_tsv(
        cds_path,
        cds
    )

    write_tsv(
        splice_path,
        splice
    )

    write_summary(
        summary_path,
        variants,
        metrics,
        sample_name
    )

    write_clinical(
        clinical_path,
        clinical_data
    )

    write_track1_audit(
        audit_path,
        variants
    )

    print()
    print("=" * 80)
    print("ARCHIVOS GENERADOS")
    print("=" * 80)

    print()
    print(
        "[OK]",
        os.path.basename(all_path)
    )

    print(
        "[OK]",
        os.path.basename(snv_path)
    )

    print(
        "[OK]",
        os.path.basename(cds_path)
    )

    print(
        "[OK]",
        os.path.basename(splice_path)
    )

    print(
        "[OK]",
        os.path.basename(summary_path)
    )

    print(
        "[OK]",
        os.path.basename(clinical_path)
    )

    print(
        "[OK]",
        os.path.basename(audit_path)
    )

    print()
    print("=" * 80)
    print("FIN DEL ZOOM FUNCIONAL BUB1B")
    print("=" * 80)

    print()
    print(
        "IMPORTANTE: este análisis describe "
        "consecuencias moleculares y no asigna "
        "patogenicidad o causalidad."
    )

    print()

    input(
        "Pulsa ENTER para cerrar..."
    )


# =============================================================================
# EJECUCIÓN
# =============================================================================

if __name__ == "__main__":
    main()
