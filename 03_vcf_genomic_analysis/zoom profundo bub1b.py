# ================================================================
# TRINOVA.OS
# ZOOM PROFUNDO BUB1B
# ================================================================
#
# REFERENCIA:
#   GRCh38
#   GENCODE v50 / Ensembl 116
#
# GEN:
#   BUB1B
#   ENSG00000156970.15
#   chr15:40160984-40221137
#   strand +
#
# OBJETIVO
# --------
# 1. Mantener el motor ternario original.
# 2. Analizar TODO BUB1B.
# 3. Separar SNV de INDEL.
# 4. Recuperar GT / AD / DP / GQ / PL.
# 5. Recuperar QUAL / FILTER / INFO.
# 6. Determinar si cada variante cae en:
#       - gene
#       - transcript
#       - exon
#       - CDS
#       - UTR
#       - intron/intergenic dentro del span
# 7. Calcular Shannon TriNova sobre SNV.
# 8. Calcular distribución ternaria.
# 9. Calcular Shannon de los alelos ALT observados.
# 10. Generar archivos TSV.
# 11. Leer clinico.json y mostrar sus fenotipos.
#
# IMPORTANTE
# ----------
# La parte REF/ALT conserva la metodología original:
#
#       A  = +1
#       T/U = -1
#       C/G = 0
#
# Los INDEL no se fuerzan artificialmente a un estado ternario.
#
# ================================================================


import os
import gzip
import json
import math
import re

from collections import Counter, defaultdict


# ================================================================
# RUTAS
# ================================================================

VCF_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\WGS_proband01.vcf.gz"
)

VCF_INDEX_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\WGS_proband01.vcf.gz.tbi"
)

GTF_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\gencode.v50.annotation.gtf.gz"
)

FASTA_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\GRCh38.primary_assembly.genome.fa.gz"
)

CLINICO_PATH = (
    r"C:\Users\ericr\AppData\Local\Programs\Python\Python313"
    r"\MVA-DATA\mva_data\clinico.json"
)


# ================================================================
# BUB1B CONFIRMADO EN GENCODE v50
# ================================================================

BUB1B_CHROM = "chr15"

BUB1B_START = 40160984

BUB1B_END = 40221137

BUB1B_GENE_ID = "ENSG00000156970.15"

BUB1B_GENE_NAME = "BUB1B"


# ================================================================
# ARCHIVOS DE SALIDA
# ================================================================

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "zoom bub1b resultados"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_ALL = os.path.join(OUTPUT_DIR, "bub1b_todas_variantes.tsv")

OUTPUT_SNV = os.path.join(OUTPUT_DIR, "bub1b_snv.tsv")

OUTPUT_INDEL = os.path.join(OUTPUT_DIR, "bub1b_indel.tsv")

OUTPUT_OBSERVED = os.path.join(OUTPUT_DIR, "bub1b_observed.tsv")

OUTPUT_EXON = os.path.join(OUTPUT_DIR, "bub1b_exon.tsv")

OUTPUT_CDS = os.path.join(OUTPUT_DIR, "bub1b_cds.tsv")

OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, "bub1b_deep_summary.txt")

OUTPUT_CLINICAL = os.path.join(OUTPUT_DIR, "bub1b_clinical.txt")


# ================================================================
# TRINOVA.OS
# CODIFICACIÓN TERNARIA ORIGINAL
# ================================================================

def nucleotide_to_trit(base):

    if base is None:
        return None

    base = base.upper()

    if base == "A":
        return 1

    if base in ("T", "U"):
        return -1

    if base in ("C", "G"):
        return 0

    return None


# ================================================================
# REF / ALT
# ================================================================

def tcmp(ref, alt):

    return (
        nucleotide_to_trit(ref),
        nucleotide_to_trit(alt)
    )


# ================================================================
# SHANNON
# ================================================================

def shannon_entropy(values):

    values = [
        v for v in values
        if v is not None
    ]

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


# ================================================================
# DISTRIBUCIÓN
# ================================================================

def ternary_distribution(values):

    values = [
        v for v in values
        if v is not None
    ]

    total = len(values)

    counts = Counter(values)

    result = {}

    for state in (-1, 0, 1):

        count = counts.get(
            state,
            0
        )

        percentage = (
            100.0 * count / total
            if total > 0
            else 0.0
        )

        result[state] = (
            count,
            percentage
        )

    return result


# ================================================================
# NORMALIZAR CROMOSOMA
# ================================================================

def normalize_chrom(chrom):

    chrom = chrom.strip()

    if chrom.lower().startswith("chr"):

        return chrom[3:]

    return chrom


# ================================================================
# PARSER GTF
# ================================================================

def parse_gtf_attributes(text):

    attributes = {}

    pattern = re.compile(
        r'(\S+)\s+"([^"]*)"'
    )

    for match in pattern.finditer(text):

        key = match.group(1)

        value = match.group(2)

        attributes[key] = value

    return attributes


# ================================================================
# CARGAR ANOTACIÓN DE BUB1B
# ================================================================

def load_bub1b_annotation():

    print()
    print("=" * 80)
    print("CARGANDO ANOTACIÓN GENCODE v50 DE BUB1B")
    print("=" * 80)
    print()

    annotations = []

    exon_intervals = []

    cds_intervals = []

    utr_intervals = []

    transcript_intervals = []

    gene_intervals = []

    with gzip.open(
        GTF_PATH,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as gtf:

        for line in gtf:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chrom = fields[0]

            if normalize_chrom(chrom) != normalize_chrom(
                BUB1B_CHROM
            ):
                continue

            try:

                start = int(fields[3])

                end = int(fields[4])

            except ValueError:

                continue

            if end < BUB1B_START:
                continue

            if start > BUB1B_END:
                continue

            feature = fields[2]

            strand = fields[6]

            attributes = parse_gtf_attributes(
                fields[8]
            )

            gene_name = attributes.get(
                "gene_name",
                ""
            )

            if gene_name != BUB1B_GENE_NAME:
                continue

            item = {

                "chrom": chrom,

                "start": start,

                "end": end,

                "feature": feature,

                "strand": strand,

                "gene_id": attributes.get(
                    "gene_id",
                    "."
                ),

                "gene_name": gene_name,

                "gene_type": attributes.get(
                    "gene_type",
                    "."
                ),

                "transcript_id": attributes.get(
                    "transcript_id",
                    "."
                ),

                "transcript_type": attributes.get(
                    "transcript_type",
                    "."
                ),

                "transcript_name": attributes.get(
                    "transcript_name",
                    "."
                ),

                "exon_number": attributes.get(
                    "exon_number",
                    "."
                ),

                "exon_id": attributes.get(
                    "exon_id",
                    "."
                ),

                "protein_id": attributes.get(
                    "protein_id",
                    "."
                ),

                "ccdsid": attributes.get(
                    "ccdsid",
                    "."
                ),

                "tag": attributes.get(
                    "tag",
                    "."
                )
            }

            annotations.append(item)

            if feature == "gene":

                gene_intervals.append(item)

            elif feature == "transcript":

                transcript_intervals.append(item)

            elif feature == "exon":

                exon_intervals.append(item)

            elif feature == "CDS":

                cds_intervals.append(item)

            elif feature == "UTR":

                utr_intervals.append(item)

    print(
        "Anotaciones BUB1B:",
        len(annotations)
    )

    print(
        "Genes:",
        len(gene_intervals)
    )

    print(
        "Transcripts:",
        len(transcript_intervals)
    )

    print(
        "Exones:",
        len(exon_intervals)
    )

    print(
        "CDS:",
        len(cds_intervals)
    )

    print(
        "UTR:",
        len(utr_intervals)
    )

    print()

    return {

        "all": annotations,

        "gene": gene_intervals,

        "transcript": transcript_intervals,

        "exon": exon_intervals,

        "cds": cds_intervals,

        "utr": utr_intervals
    }


# ================================================================
# BUSCAR ANOTACIÓN PARA UNA POSICIÓN
# ================================================================

def position_in_interval(pos, interval):

    return (
        interval["start"]
        <= pos
        <= interval["end"]
    )


def annotate_position(pos, annotation):

    exon_hits = []

    cds_hits = []

    utr_hits = []

    transcript_hits = []

    for item in annotation["exon"]:

        if position_in_interval(
            pos,
            item
        ):

            exon_hits.append(item)

    for item in annotation["cds"]:

        if position_in_interval(
            pos,
            item
        ):

            cds_hits.append(item)

    for item in annotation["utr"]:

        if position_in_interval(
            pos,
            item
        ):

            utr_hits.append(item)

    for item in annotation["transcript"]:

        if position_in_interval(
            pos,
            item
        ):

            transcript_hits.append(item)

    if exon_hits:

        region = "exon"

    elif cds_hits:

        region = "CDS"

    elif utr_hits:

        region = "UTR"

    elif transcript_hits:

        region = "transcript_intron"

    else:

        region = "gene_span_other"

    exon_numbers = sorted(
        set(
            x["exon_number"]
            for x in exon_hits
            if x["exon_number"] != "."
        )
    )

    transcript_ids = sorted(
        set(
            x["transcript_id"]
            for x in (
                exon_hits
                + cds_hits
                + utr_hits
                + transcript_hits
            )
            if x["transcript_id"] != "."
        )
    )

    protein_ids = sorted(
        set(
            x["protein_id"]
            for x in cds_hits
            if x["protein_id"] != "."
        )
    )

    return {

        "region": region,

        "exon_numbers": ",".join(
            exon_numbers
        ),

        "transcript_ids": ",".join(
            transcript_ids
        ),

        "protein_ids": ",".join(
            protein_ids
        ),

        "exon_count": len(
            exon_hits
        ),

        "cds_count": len(
            cds_hits
        ),

        "utr_count": len(
            utr_hits
        )
    }


# ================================================================
# PARSEAR INFO VCF
# ================================================================

def parse_info(text):

    result = {}

    if text == ".":
        return result

    for item in text.split(";"):

        if "=" in item:

            key, value = item.split(
                "=",
                1
            )

            result[key] = value

        else:

            result[item] = True

    return result


# ================================================================
# PARSEAR SAMPLE
# ================================================================

def parse_sample(
    format_string,
    sample_string
):

    result = {}

    if (
        format_string == "."
        or sample_string == "."
    ):

        return result

    keys = format_string.split(":")

    values = sample_string.split(":")

    for i, key in enumerate(keys):

        if i < len(values):

            result[key] = values[i]

        else:

            result[key] = "."

    return result


# ================================================================
# PARSEAR GT
# ================================================================

def parse_gt(gt):

    if not gt or gt == ".":

        return []

    gt = gt.replace(
        "|",
        "/"
    )

    result = []

    for item in gt.split("/"):

        if item == ".":

            result.append(None)

        else:

            try:

                result.append(
                    int(item)
                )

            except ValueError:

                result.append(None)

    return result


# ================================================================
# ¿ALT OBSERVADO?
# ================================================================

def alt_observed(
    gt,
    alt_index
):

    alleles = parse_gt(gt)

    return (
        alt_index in alleles
    )


# ================================================================
# ESCANEAR VCF
# ================================================================

def scan_vcf(annotation):

    print()
    print("=" * 80)
    print("ESCANEANDO VCF — BUB1B")
    print("=" * 80)
    print()

    variants = []

    sample_name = "."

    total_records = 0

    interval_records = 0

    snv_count = 0

    indel_count = 0

    multiallelic_count = 0

    with gzip.open(
        VCF_PATH,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as vcf:

        for line in vcf:

            if line.startswith("##"):

                continue

            if line.startswith("#CHROM"):

                header = line.rstrip(
                    "\n"
                ).split("\t")

                if len(header) >= 10:

                    sample_name = header[9]

                continue

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) < 8:

                continue

            total_records += 1

            chrom = fields[0]

            if normalize_chrom(chrom) != normalize_chrom(
                BUB1B_CHROM
            ):

                continue

            try:

                pos = int(fields[1])

            except ValueError:

                continue

            if (
                pos < BUB1B_START
                or pos > BUB1B_END
            ):

                continue

            interval_records += 1

            ref = fields[3]

            alt_string = fields[4]

            qual = fields[5]

            filt = fields[6]

            info_string = fields[7]

            info = parse_info(
                info_string
            )

            format_string = (
                fields[8]
                if len(fields) >= 9
                else "."
            )

            sample_string = (
                fields[9]
                if len(fields) >= 10
                else "."
            )

            sample = parse_sample(
                format_string,
                sample_string
            )

            gt = sample.get(
                "GT",
                "."
            )

            ad = sample.get(
                "AD",
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

            pl = sample.get(
                "PL",
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

            alts = alt_string.split(",")

            if len(alts) > 1:

                multiallelic_count += 1

            for alt_index, alt in enumerate(
                alts,
                start=1
            ):

                is_snv = (
                    len(ref) == 1
                    and len(alt) == 1
                )

                if is_snv:

                    snv_count += 1

                    ref_trit, alt_trit = tcmp(
                        ref,
                        alt
                    )

                    trit = (
                        f"{ref_trit}->{alt_trit}"
                    )

                else:

                    indel_count += 1

                    ref_trit = None

                    alt_trit = None

                    trit = "NA"

                annotation_info = annotate_position(
                    pos,
                    annotation
                )

                observed = alt_observed(
                    gt,
                    alt_index
                )

                variants.append({

                    "CHROM": chrom,

                    "POS": pos,

                    "REF": ref,

                    "ALT": alt,

                    "ALT_INDEX": alt_index,

                    "TYPE": (
                        "SNV"
                        if is_snv
                        else "INDEL_OTHER"
                    ),

                    "QUAL": qual,

                    "FILTER": filt,

                    "GT": gt,

                    "AD": ad,

                    "DP": dp,

                    "GQ": gq,

                    "PL": pl,

                    "PGT": pgt,

                    "PID": pid,

                    "AF": info.get(
                        "AF",
                        "."
                    ),

                    "AC": info.get(
                        "AC",
                        "."
                    ),

                    "AN": info.get(
                        "AN",
                        "."
                    ),

                    "INFO_DP": info.get(
                        "DP",
                        "."
                    ),

                    "MQ": info.get(
                        "MQ",
                        "."
                    ),

                    "QD": info.get(
                        "QD",
                        "."
                    ),

                    "FS": info.get(
                        "FS",
                        "."
                    ),

                    "SOR": info.get(
                        "SOR",
                        "."
                    ),

                    "REF_TRIT": ref_trit,

                    "ALT_TRIT": alt_trit,

                    "TRIT": trit,

                    "OBSERVED": observed,

                    "REGION": annotation_info[
                        "region"
                    ],

                    "EXON_NUMBERS": annotation_info[
                        "exon_numbers"
                    ],

                    "TRANSCRIPT_IDS": annotation_info[
                        "transcript_ids"
                    ],

                    "PROTEIN_IDS": annotation_info[
                        "protein_ids"
                    ],

                    "EXON_COUNT": annotation_info[
                        "exon_count"
                    ],

                    "CDS_COUNT": annotation_info[
                        "cds_count"
                    ],

                    "UTR_COUNT": annotation_info[
                        "utr_count"
                    ]
                })

    print(
        "Muestra:",
        sample_name
    )

    print()

    print(
        "Registros VCF recorridos:",
        total_records
    )

    print(
        "Dentro de BUB1B:",
        interval_records
    )

    print(
        "SNV:",
        snv_count
    )

    print(
        "INDEL/otros:",
        indel_count
    )

    print(
        "Multialélicos:",
        multiallelic_count
    )

    print()

    return (
        variants,
        sample_name
    )


# ================================================================
# ANALISIS TRINARIO
# ================================================================

def analyze_ternary(variants):

    values = []

    transitions = Counter()

    snvs = []

    for v in variants:

        if v["TYPE"] != "SNV":
            continue

        if (
            v["REF_TRIT"] is None
            or v["ALT_TRIT"] is None
        ):
            continue

        values.append(
            v["REF_TRIT"]
        )

        values.append(
            v["ALT_TRIT"]
        )

        transitions[
            v["TRIT"]
        ] += 1

        snvs.append(v)

    entropy = shannon_entropy(
        values
    )

    distribution = ternary_distribution(
        values
    )

    return {

        "values": values,

        "entropy": entropy,

        "distribution": distribution,

        "transitions": transitions,

        "snvs": snvs
    }


# ================================================================
# ANALISIS OBSERVADO
# ================================================================

def analyze_observed(variants):

    values = []

    observed = []

    for v in variants:

        if not v["OBSERVED"]:
            continue

        if v["TYPE"] != "SNV":
            continue

        if v["ALT_TRIT"] is None:
            continue

        values.append(
            v["ALT_TRIT"]
        )

        observed.append(v)

    entropy = shannon_entropy(
        values
    )

    distribution = ternary_distribution(
        values
    )

    return {

        "values": values,

        "entropy": entropy,

        "distribution": distribution,

        "variants": observed
    }


# ================================================================
# ESCRIBIR TSV
# ================================================================

COLUMNS = [

    "CHROM",
    "POS",
    "REF",
    "ALT",
    "ALT_INDEX",
    "TYPE",

    "QUAL",
    "FILTER",

    "GT",
    "AD",
    "DP",
    "GQ",
    "PL",
    "PGT",
    "PID",

    "AF",
    "AC",
    "AN",

    "INFO_DP",
    "MQ",
    "QD",
    "FS",
    "SOR",

    "REF_TRIT",
    "ALT_TRIT",
    "TRIT",

    "OBSERVED",

    "REGION",

    "EXON_NUMBERS",

    "TRANSCRIPT_IDS",

    "PROTEIN_IDS",

    "EXON_COUNT",

    "CDS_COUNT",

    "UTR_COUNT"
]


def write_tsv(
    path,
    variants
):

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as out:

        out.write(
            "\t".join(COLUMNS)
            + "\n"
        )

        for v in variants:

            row = []

            for column in COLUMNS:

                value = v.get(
                    column,
                    "."
                )

                if value is None:

                    value = "."

                row.append(
                    str(value)
                )

            out.write(
                "\t".join(row)
                + "\n"
            )


# ================================================================
# CLÍNICO
# ================================================================

def extract_phenotypes(obj):

    results = []

    if isinstance(obj, dict):

        for key, value in obj.items():

            key_lower = str(
                key
            ).lower()

            if (
                "phenotype" in key_lower
                or "hpo" in key_lower
                or "clinical" in key_lower
                or "manifest" in key_lower
            ):

                results.append(
                    (
                        key,
                        value
                    )
                )

            results.extend(
                extract_phenotypes(value)
            )

    elif isinstance(obj, list):

        for item in obj:

            results.extend(
                extract_phenotypes(item)
            )

    return results


def inspect_clinical():

    print()
    print("=" * 80)
    print("CLINICO.JSON")
    print("=" * 80)
    print()

    if not os.path.exists(
        CLINICO_PATH
    ):

        print(
            "No se encontró clinico.json."
        )

        return {}

    try:

        with open(
            CLINICO_PATH,
            "r",
            encoding="utf-8"
        ) as f:

            data = json.load(f)

    except Exception as e:

        print(
            "ERROR leyendo clinico.json:"
        )

        print(e)

        return {}

    print(
        "JSON cargado correctamente."
    )

    print()

    if isinstance(
        data,
        dict
    ):

        print(
            "Claves principales:"
        )

        for key in data.keys():

            print(
                "  -",
                key
            )

    print()

    phenotypes = extract_phenotypes(
        data
    )

    print(
        "Campos clínicos detectados:",
        len(phenotypes)
    )

    for key, value in phenotypes:

        print()

        print(
            "Campo:",
            key
        )

        print(
            "Valor:",
            value
        )

    # Guardar copia legible.
    with open(
        OUTPUT_CLINICAL,
        "w",
        encoding="utf-8"
    ) as out:

        json.dump(
            data,
            out,
            ensure_ascii=False,
            indent=2
        )

    return data


# ================================================================
# TABLA DE RESULTADOS
# ================================================================

def print_variants(
    variants
):

    print()
    print("=" * 150)
    print("VARIANTES EN BUB1B")
    print("=" * 150)
    print()

    header = (
        f"{'POS':>10} "
        f"{'REF':>5} "
        f"{'ALT':>5} "
        f"{'TYPE':>11} "
        f"{'GT':>7} "
        f"{'DP':>5} "
        f"{'GQ':>5} "
        f"{'AF':>8} "
        f"{'REGION':>18} "
        f"{'TRIT':>8} "
        f"{'OBS':>5}"
    )

    print(header)

    print("-" * len(header))

    for v in variants:

        print(
            f"{v['POS']:>10} "
            f"{v['REF']:>5} "
            f"{v['ALT']:>5} "
            f"{v['TYPE']:>11} "
            f"{v['GT']:>7} "
            f"{v['DP']:>5} "
            f"{v['GQ']:>5} "
            f"{str(v['AF']):>8} "
            f"{v['REGION']:>18} "
            f"{v['TRIT']:>8} "
            f"{str(v['OBSERVED']):>5}"
        )


# ================================================================
# RESUMEN
# ================================================================

def write_summary(
    path,
    variants,
    sample_name,
    ternary,
    observed
):

    max_entropy = math.log2(3)

    coherence = (
        1.0
        - ternary["entropy"]
        / max_entropy
    )

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as out:

        out.write(
            "TRINOVA.OS — DEEP ZOOM BUB1B\n"
        )

        out.write(
            "=" * 80
            + "\n\n"
        )

        out.write(
            "GEN\n"
        )

        out.write(
            f"Gene: {BUB1B_GENE_NAME}\n"
        )

        out.write(
            f"Gene ID: {BUB1B_GENE_ID}\n"
        )

        out.write(
            f"Region: {BUB1B_CHROM}:"
            f"{BUB1B_START}-"
            f"{BUB1B_END}\n"
        )

        out.write(
            "Reference: GRCh38\n"
        )

        out.write(
            "Annotation: GENCODE v50\n\n"
        )

        out.write(
            "SAMPLE\n"
        )

        out.write(
            f"{sample_name}\n\n"
        )

        out.write(
            "VARIANTS\n"
        )

        out.write(
            f"Total: {len(variants)}\n"
        )

        out.write(
            "SNV: "
            + str(
                sum(
                    v["TYPE"] == "SNV"
                    for v in variants
                )
            )
            + "\n"
        )

        out.write(
            "INDEL/other: "
            + str(
                sum(
                    v["TYPE"] != "SNV"
                    for v in variants
                )
            )
            + "\n\n"
        )

        out.write(
            "TRINOVA — REF/ALT\n"
        )

        out.write(
            f"Shannon: "
            f"{ternary['entropy']:.9f}\n"
        )

        out.write(
            f"Maximum theoretical: "
            f"{max_entropy:.9f}\n"
        )

        out.write(
            f"Relative coherence: "
            f"{coherence:.9f}\n\n"
        )

        out.write(
            "TERNARY DISTRIBUTION\n"
        )

        for state in (-1, 0, 1):

            count, percentage = (
                ternary["distribution"][state]
            )

            out.write(
                f"{state:+d}: "
                f"{count} "
                f"({percentage:.6f}%)\n"
            )

        out.write("\n")

        out.write(
            "TRANSITIONS\n"
        )

        for transition, count in sorted(
            ternary["transitions"].items()
        ):

            out.write(
                f"{transition}: {count}\n"
            )

        out.write("\n")

        out.write(
            "OBSERVED SAMPLE — ALT TRITS\n"
        )

        out.write(
            f"Observed SNVs: "
            f"{len(observed['variants'])}\n"
        )

        out.write(
            f"Observed Shannon: "
            f"{observed['entropy']:.9f}\n\n"
        )

        out.write(
            "OBSERVED DISTRIBUTION\n"
        )

        for state in (-1, 0, 1):

            count, percentage = (
                observed["distribution"][state]
            )

            out.write(
                f"{state:+d}: "
                f"{count} "
                f"({percentage:.6f}%)\n"
            )

        out.write("\n")

        out.write(
            "REGIONS\n"
        )

        region_counts = Counter(
            v["REGION"]
            for v in variants
        )

        for region, count in (
            region_counts.most_common()
        ):

            out.write(
                f"{region}: {count}\n"
            )


# ================================================================
# MAIN
# ================================================================

def main():

    print()
    print("=" * 80)
    print("TRINOVA.OS")
    print("DEEP ZOOM BUB1B")
    print("=" * 80)
    print()

    # ------------------------------------------------------------
    # ARCHIVOS
    # ------------------------------------------------------------

    print(
        "COMPROBACIÓN DE ARCHIVOS"
    )

    print()

    files = [

        (
            "VCF",
            VCF_PATH
        ),

        (
            "VCF INDEX",
            VCF_INDEX_PATH
        ),

        (
            "GTF",
            GTF_PATH
        ),

        (
            "FASTA",
            FASTA_PATH
        ),

        (
            "CLINICO",
            CLINICO_PATH
        )
    ]

    for label, path in files:

        if os.path.exists(path):

            print(
                "[OK]",
                label
            )

        else:

            print(
                "[NO ENCONTRADO]",
                label
            )

    print()

    if not os.path.exists(
        VCF_PATH
    ):

        print(
            "ERROR: falta el VCF."
        )

        return

    if not os.path.exists(
        GTF_PATH
    ):

        print(
            "ERROR: falta el GTF."
        )

        return

    # ------------------------------------------------------------
    # ANOTACIÓN
    # ------------------------------------------------------------

    annotation = (
        load_bub1b_annotation()
    )

    # ------------------------------------------------------------
    # VCF
    # ------------------------------------------------------------

    variants, sample_name = (
        scan_vcf(
            annotation
        )
    )

    # ------------------------------------------------------------
    # TRINARIO
    # ------------------------------------------------------------

    ternary = analyze_ternary(
        variants
    )

    # ------------------------------------------------------------
    # OBSERVADO
    # ------------------------------------------------------------

    observed = analyze_observed(
        variants
    )

    # ------------------------------------------------------------
    # CLÍNICO
    # ------------------------------------------------------------

    clinical = inspect_clinical()

    # ------------------------------------------------------------
    # TABLA
    # ------------------------------------------------------------

    print_variants(
        variants
    )

    # ------------------------------------------------------------
    # ARCHIVOS
    # ------------------------------------------------------------

    snv = [
        v for v in variants
        if v["TYPE"] == "SNV"
    ]

    indel = [
        v for v in variants
        if v["TYPE"] != "SNV"
    ]

    observed_all = [
        v for v in variants
        if v["OBSERVED"]
    ]

    exon = [
        v for v in variants
        if v["EXON_COUNT"] > 0
    ]

    cds = [
        v for v in variants
        if v["CDS_COUNT"] > 0
    ]

    write_tsv(
        OUTPUT_ALL,
        variants
    )

    write_tsv(
        OUTPUT_SNV,
        snv
    )

    write_tsv(
        OUTPUT_INDEL,
        indel
    )

    write_tsv(
        OUTPUT_OBSERVED,
        observed_all
    )

    write_tsv(
        OUTPUT_EXON,
        exon
    )

    write_tsv(
        OUTPUT_CDS,
        cds
    )

    write_summary(
        OUTPUT_SUMMARY,
        variants,
        sample_name,
        ternary,
        observed
    )

    # ------------------------------------------------------------
    # RESULTADO FINAL
    # ------------------------------------------------------------

    print()
    print("=" * 80)
    print("RESULTADO FINAL")
    print("=" * 80)
    print()

    print(
        "BUB1B:",
        BUB1B_CHROM,
        BUB1B_START,
        "-",
        BUB1B_END
    )

    print()

    print(
        "Total variantes:",
        len(variants)
    )

    print(
        "SNV:",
        len(snv)
    )

    print(
        "INDEL/otros:",
        len(indel)
    )

    print(
        "Observadas por GT:",
        len(observed_all)
    )

    print()

    print(
        "Shannon TriNova REF/ALT:",
        f"{ternary['entropy']:.9f}"
    )

    print(
        "Shannon observado:",
        f"{observed['entropy']:.9f}"
    )

    print()

    print(
        "Variantes en EXON:",
        len(exon)
    )

    print(
        "Variantes en CDS:",
        len(cds)
    )

    print()

    print(
        "ARCHIVOS GENERADOS"
    )

    print()

    for path in [

        OUTPUT_ALL,
        OUTPUT_SNV,
        OUTPUT_INDEL,
        OUTPUT_OBSERVED,
        OUTPUT_EXON,
        OUTPUT_CDS,
        OUTPUT_SUMMARY,
        OUTPUT_CLINICAL

    ]:

        print(
            "[OK]",
            path
        )

    print()
    print("=" * 80)
    print("FIN DEL DEEP ZOOM BUB1B")
    print("=" * 80)
    print()


# ================================================================
# EJECUTAR
# ================================================================

if __name__ == "__main__":

    main()
