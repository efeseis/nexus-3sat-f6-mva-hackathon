import os
import re
import gzip
import json
import math
import csv
from collections import defaultdict

# ============================================================
# TRINOVA.OS R7.1
# CHR15 VARIANT ARCHITECTURE + FUNCTIONAL TRIAGE
#
# Objetivos:
#   1. Separar allele balance de population AF.
#   2. Anotar localmente función usando GENCODE v50.
#   3. No asumir que dos hets = compound-het.
#   4. Reducir el dominio de pares candidatos.
#   
#
# IMPORTANTE:
#   - phase = UNKNOWN
#   - population AF puede ser UNKNOWN
#   - esto NO demuestra causalidad clínica
# ============================================================

VCF_PATH = r"C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\proband01.vcf.gz"
GTF_PATH = r"C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\gencode.v50.annotation.gtf.gz"
CLINICAL_PATH = r"C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\clinico.json"


OUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "variantes_arch_data"
)

TARGET_CHROM = "15"

# ============================================================
# CONFIG
# ============================================================

MAX_GENE_PAIRS_PER_GENE = 100
TOP_GENES = 100
TOP_PAIRS = 200

# Ventana aproximada para splice region
SPLICE_WINDOW = 8

# Promotor aproximado para priorización local
PROMOTER_WINDOW = 2000

# ============================================================
# UTILIDADES
# ============================================================

def safe_float(x):
    try:
        if x is None or x == "." or x == "":
            return None
        return float(x)
    except Exception:
        return None


def safe_int(x):
    try:
        if x is None or x == "." or x == "":
            return None
        return int(x)
    except Exception:
        return None


def normalize_chrom(chrom):
    chrom = chrom.strip()
    chrom = chrom.replace("chr", "").replace("CHR_", "")
    return chrom


def parse_info(info):
    out = {}

    if info == ".":
        return out

    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = True

    return out


def parse_gt(sample_value, format_keys):
    data = dict(zip(format_keys, sample_value.split(":")))
    gt = data.get("GT", ".")

    alleles = re.split(r"[\/|]", gt)

    if len(alleles) != 2:
        return data, None

    try:
        a = int(alleles[0])
        b = int(alleles[1])
    except Exception:
        return data, None

    return data, (a, b)


def genotype_class(gt):
    if gt is None:
        return "UNKNOWN"

    a, b = gt

    if a == 0 and b == 0:
        return "HOM_REF"

    if (a == 0 and b > 0) or (a > 0 and b == 0):
        return "HET"

    if a > 0 and b > 0:
        return "HOM_ALT"

    return "UNKNOWN"


def extract_ad(ad_string):
    if ad_string in (None, "", "."):
        return None

    try:
        return [int(x) for x in ad_string.split(",")]
    except Exception:
        return None


def allele_balance_from_ad(ad, alt_index=1):
    """
    Allele balance del ALT dentro del individuo.

    Esto NO es population AF.
    """
    if not ad:
        return None

    if len(ad) <= alt_index:
        return None

    ref_depth = ad[0]
    alt_depth = ad[alt_index]

    total = ref_depth + alt_depth

    if total <= 0:
        return None

    return alt_depth / total


def population_af_from_info(info_dict, alt_index=1):
    """
    Solo acepta AF explícitamente presente en INFO.

    NO usa AD/DP como sustituto.
    """
    raw = info_dict.get("AF")

    if raw is None:
        return None

    values = raw.split(",")

    if alt_index - 1 >= len(values):
        return None

    return safe_float(values[alt_index - 1])


def parse_gtf_attributes(attr_string):
    attrs = {}

    for match in re.finditer(
        r'(\S+)\s+"([^"]*)"',
        attr_string
    ):
        key = match.group(1)
        value = match.group(2)
        attrs[key] = value

    return attrs


# ============================================================
# CLÍNICA
# ============================================================

def load_clinical(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    hpos = []

    def recursive_find(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, str):
                    if re.fullmatch(r"HP:\d{7}", v):
                        hpos.append(v)
                recursive_find(v)

        elif isinstance(obj, list):
            for item in obj:
                recursive_find(item)

    recursive_find(data)

    return sorted(set(hpos))


# ============================================================
# GENCODE
# ============================================================

def load_chr15_gencode(path):
    """
    Carga genes, transcripts, exons y CDS de chr15.

    No intenta convertir esto en una predicción clínica.
    Solo proporciona contexto funcional/local.
    """

    genes = {}
    transcripts = {}
    exons = defaultdict(list)
    cds = defaultdict(list)

    print("    Leyendo GENCODE v50 chr15...")

    with gzip.open(path, "rt", encoding="utf-8") as f:

        for line in f:

            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            chrom = normalize_chrom(parts[0])

            if chrom != TARGET_CHROM:
                continue

            feature = parts[2]

            try:
                start = int(parts[3])
                end = int(parts[4])
            except Exception:
                continue

            strand = parts[6]
            attrs = parse_gtf_attributes(parts[8])

            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)

            transcript_id = attrs.get(
                "transcript_id",
                ""
            ).split(".")[0]

            if not gene_id:
                continue

            if feature == "gene":

                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "start": start,
                    "end": end,
                    "strand": strand
                }

            elif feature == "transcript":

                if transcript_id:

                    transcripts[transcript_id] = {
                        "transcript_id": transcript_id,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "start": start,
                        "end": end,
                        "strand": strand
                    }

            elif feature == "exon":

                if transcript_id:
                    exons[transcript_id].append(
                        (start, end)
                    )

            elif feature == "CDS":

                if transcript_id:
                    cds[transcript_id].append(
                        (start, end)
                    )

    return genes, transcripts, exons, cds


# ============================================================
# INDEX GENÓMICO
# ============================================================

def build_gene_index(genes, bin_size=1_000_000):

    index = defaultdict(list)

    for gene_id, g in genes.items():

        b1 = g["start"] // bin_size
        b2 = g["end"] // bin_size

        for b in range(b1, b2 + 1):
            index[b].append(gene_id)

    return index


def interval_overlap(pos, ref, start, end):

    var_start = pos
    var_end = pos + max(len(ref), 1) - 1

    return not (
        var_end < start or
        var_start > end
    )


def annotate_gene_region(
    pos,
    ref,
    gene,
    transcripts,
    exons,
    cds
):

    if not interval_overlap(
        pos,
        ref,
        gene["start"],
        gene["end"]
    ):
        return None

    region_types = set()
    transcript_hits = set()
    cds_hits = set()

    # Promotor aproximado
    if gene["strand"] == "+":

        if gene["start"] - PROMOTER_WINDOW <= pos < gene["start"]:
            region_types.add("PROMOTER")

    else:

        if gene["end"] < pos <= gene["end"] + PROMOTER_WINDOW:
            region_types.add("PROMOTER")

    for tid, t in transcripts.items():

        if t["gene_id"] != gene["gene_id"]:
            continue

        if not interval_overlap(
            pos,
            ref,
            t["start"],
            t["end"]
        ):
            continue

        transcript_hits.add(tid)

        exon_hit = False
        cds_hit = False
        splice_hit = False

        for ex_start, ex_end in exons.get(tid, []):

            if interval_overlap(
                pos,
                ref,
                ex_start,
                ex_end
            ):

                exon_hit = True

                if (
                    abs(pos - ex_start) <= SPLICE_WINDOW or
                    abs(pos - ex_end) <= SPLICE_WINDOW
                ):
                    splice_hit = True

        for cds_start, cds_end in cds.get(tid, []):

            if interval_overlap(
                pos,
                ref,
                cds_start,
                cds_end
            ):
                cds_hit = True
                cds_hits.add(tid)

        if splice_hit:
            region_types.add("SPLICE_REGION")

        elif cds_hit:
            region_types.add("CDS")

        elif exon_hit:
            region_types.add("EXON")

        else:
            region_types.add("INTRON")

    if not region_types:
        return None

    if "SPLICE_REGION" in region_types:
        priority = 5
        label = "SPLICE_REGION"

    elif "CDS" in region_types:
        priority = 4
        label = "CDS"

    elif "EXON" in region_types:
        priority = 3
        label = "EXON"

    elif "PROMOTER" in region_types:
        priority = 2
        label = "PROMOTER"

    elif "INTRON" in region_types:
        priority = 1
        label = "INTRON"

    else:
        priority = 0
        label = "OTHER"

    return {
        "region": label,
        "functional_priority": priority,
        "transcripts": sorted(transcript_hits),
        "cds_transcripts": sorted(cds_hits)
    }


def find_gene_annotations(
    pos,
    ref,
    genes,
    index
):

    bin_size = 1_000_000
    b = pos // bin_size

    candidates = []

    for gene_id in index.get(b, []):

        g = genes[gene_id]

        ann = annotate_gene_region(
            pos,
            ref,
            g,
            GLOBAL_TRANSCRIPTS,
            GLOBAL_EXONS,
            GLOBAL_CDS
        )

        if ann:

            row = {
                "gene_id": gene_id,
                "gene_name": g["gene_name"],
                "gene_start": g["start"],
                "gene_end": g["end"],
                "gene_strand": g["strand"],
                **ann
            }

            candidates.append(row)

    return candidates


# ============================================================
# VARIANT SCORE
# ============================================================

def variant_score(v):

    score = 0.0

    # Calidad
    qual = v.get("qual")

    if qual is not None:

        score += min(
            max(qual / 100.0, 0.0),
            2.0
        )

    # GQ
    gq = v.get("gq")

    if gq is not None:

        score += min(
            max(gq / 50.0, 0.0),
            2.0
        )

    # Functional local annotation
    score += v.get(
        "functional_priority",
        0
    ) * 1.25

    # Allele balance:
    # no es AF poblacional
    ab = v.get("allele_balance")

    if ab is not None:

        if 0.25 <= ab <= 0.75:
            score += 0.75

        elif 0.15 <= ab <= 0.85:
            score += 0.25

    # Population AF:
    # solo si realmente está presente
    paf = v.get("population_af")

    if paf is not None:

        if paf < 0.001:
            score += 2.0

        elif paf < 0.01:
            score += 1.0

        elif paf < 0.05:
            score += 0.25

    return score


# ============================================================
# VCF CHR15
# ============================================================

def read_chr15_vcf(
    path,
    genes,
    gene_index
):

    variants = []

    sample_name = None

    print("    Procesando VCF chr15...")

    with gzip.open(
        path,
        "rt",
        encoding="utf-8"
    ) as f:

        for line in f:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                header = line.rstrip("\n").split("\t")

                if len(header) >= 10:
                    sample_name = header[9]

                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 10:
                continue

            chrom = normalize_chrom(parts[0])

            if chrom != TARGET_CHROM:
                continue

            try:
                pos = int(parts[1])
            except Exception:
                continue

            ref = parts[3]
            alt_field = parts[4]

            qual = safe_float(parts[5])

            info = parse_info(parts[7])

            format_keys = parts[8].split(":")

            sample_data = parts[9]

            sample_values, gt = parse_gt(
                sample_data,
                format_keys
            )

            gt_class = genotype_class(gt)

            if gt_class not in (
                "HET",
                "HOM_ALT"
            ):
                continue

            ad = extract_ad(
                sample_values.get("AD")
            )

            dp = safe_int(
                sample_values.get("DP")
            )

            gq = safe_float(
                sample_values.get("GQ")
            )

            alts = alt_field.split(",")

            for alt_index, alt in enumerate(alts, start=1):

                # Si genotype no contiene este ALT
                if gt is not None:

                    if alt_index not in gt:
                        continue

                allele_balance = allele_balance_from_ad(
                    ad,
                    alt_index
                )

                population_af = population_af_from_info(
                    info,
                    alt_index
                )

                annotations = find_gene_annotations(
                    pos,
                    ref,
                    genes,
                    gene_index
                )

                # Si no cae dentro de gen, conservar
                # como intergenic.
                if not annotations:

                    annotations = [{
                        "gene_id": "",
                        "gene_name": "INTERGENIC",
                        "gene_start": None,
                        "gene_end": None,
                        "gene_strand": ".",
                        "region": "INTERGENIC",
                        "functional_priority": 0,
                        "transcripts": [],
                        "cds_transcripts": []
                    }]

                for ann in annotations:

                    v = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,

                        "gt": gt,
                        "gt_class": gt_class,

                        "dp": dp,
                        "gq": gq,
                        "qual": qual,

                        # IMPORTANTE:
                        "allele_balance": allele_balance,
                        "population_af": population_af,

                        "gene_id": ann["gene_id"],
                        "gene_name": ann["gene_name"],

                        "region": ann["region"],
                        "functional_priority":
                            ann["functional_priority"],

                        "transcripts":
                            ",".join(
                                ann["transcripts"]
                            ),

                        "cds_transcripts":
                            ",".join(
                                ann["cds_transcripts"]
                            )
                    }

                    v["score"] = variant_score(v)

                    variants.append(v)

    return variants, sample_name


# ============================================================
# ARQUITECTURA POR GEN
# ============================================================

def aggregate_genes(variants):

    by_gene = defaultdict(list)

    for v in variants:

        gene_id = v["gene_id"]

        if not gene_id:
            continue

        by_gene[gene_id].append(v)

    rows = []

    for gene_id, vals in by_gene.items():

        hets = [
            v for v in vals
            if v["gt_class"] == "HET"
        ]

        homs = [
            v for v in vals
            if v["gt_class"] == "HOM_ALT"
        ]

        functional = [
            v for v in vals
            if v["functional_priority"] >= 3
        ]

        splice = [
            v for v in vals
            if v["functional_priority"] >= 5
        ]

        pop_af_known = [
            v["population_af"]
            for v in vals
            if v["population_af"] is not None
        ]

        rare = [
            v for v in vals
            if (
                v["population_af"] is not None and
                v["population_af"] < 0.01
            )
        ]

        # Cuando no hay population AF,
        # NO llamamos al variante "rare".
        af_unknown = [
            v for v in vals
            if v["population_af"] is None
        ]

        max_score = max(
            (v["score"] for v in vals),
            default=0.0
        )

        functional_scores = sorted(
            (
                v["score"]
                for v in functional
            ),
            reverse=True
        )

        # Evita que cantidad bruta de variantes
        # domine el ranking.
        architecture = max_score

        if len(hets) >= 2:
            architecture += 1.0

        if len(functional) >= 1:
            architecture += 1.5

        if len(functional) >= 2:
            architecture += 1.5

        if len(splice) >= 1:
            architecture += 1.0

        if len(rare) >= 1:
            architecture += 1.0

        # Bonus solo cuando realmente hay
        # dos variantes funcionales het.
        functional_hets = [
            v for v in hets
            if v["functional_priority"] >= 3
        ]

        if len(functional_hets) >= 2:
            architecture += 2.0

        gene_name = vals[0]["gene_name"]

        rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,

            "variants": len(vals),
            "hets": len(hets),
            "hom_alt": len(homs),

            "functional": len(functional),
            "splice": len(splice),

            "rare_known": len(rare),
            "af_unknown": len(af_unknown),

            "max_variant_score": round(
                max_score,
                6
            ),

            "architecture_score": round(
                architecture,
                6
            )
        })

    rows.sort(
        key=lambda x: (
            x["architecture_score"],
            x["functional"],
            x["hets"]
        ),
        reverse=True
    )

    for i, row in enumerate(rows, start=1):
        row["rank"] = i

    return rows


# ============================================================
# COMPOUND-HET TRIAGE
# ============================================================

def pair_level(a, b):

    # Level 1:
    # dos heterocigotas
    level = 1

    # Level 2:
    # ambas con population AF conocida y <1%
    if (
        a["population_af"] is not None and
        b["population_af"] is not None and
        a["population_af"] < 0.01 and
        b["population_af"] < 0.01
    ):
        level = 2

    # Level 3:
    # además ambas funcionales
    if (
        level >= 2 and
        a["functional_priority"] >= 3 and
        b["functional_priority"] >= 3
    ):
        level = 3

    return level


def pair_score(a, b):

    level = pair_level(a, b)

    score = (
        a["score"] +
        b["score"]
    )

    score += level * 1.5

    # Bonus si ambas tienen CDS/splice
    if (
        a["functional_priority"] >= 4 and
        b["functional_priority"] >= 4
    ):
        score += 2.0

    # Balance razonable
    if (
        a["allele_balance"] is not None and
        b["allele_balance"] is not None
    ):

        if (
            0.20 <= a["allele_balance"] <= 0.80 and
            0.20 <= b["allele_balance"] <= 0.80
        ):
            score += 0.5

    return score


def build_compound_pairs(variants):

    by_gene = defaultdict(list)

    for v in variants:

        if v["gt_class"] == "HET":
            if v["gene_id"]:
                by_gene[v["gene_id"]].append(v)

    pairs = []

    for gene_id, vals in by_gene.items():

        # Solo genes con al menos 2 hets
        if len(vals) < 2:
            continue

        # Priorizamos primero funcionales.
        functional = [
            v for v in vals
            if v["functional_priority"] >= 3
        ]

        if len(functional) >= 2:
            pool = functional
        else:
            pool = vals

        n = len(pool)

        for i in range(n):

            for j in range(i + 1, n):

                a = pool[i]
                b = pool[j]

                score = pair_score(a, b)

                level = pair_level(a, b)

                pairs.append({
                    "gene_id": gene_id,
                    "gene_name": a["gene_name"],

                    "level": level,

                    "score": round(
                        score,
                        6
                    ),

                    "chrom_1": a["chrom"],
                    "pos_1": a["pos"],
                    "ref_1": a["ref"],
                    "alt_1": a["alt"],

                    "chrom_2": b["chrom"],
                    "pos_2": b["pos"],
                    "ref_2": b["ref"],
                    "alt_2": b["alt"],

                    "region_1": a["region"],
                    "region_2": b["region"],

                    "population_af_1":
                        a["population_af"],

                    "population_af_2":
                        b["population_af"],

                    "allele_balance_1":
                        a["allele_balance"],

                    "allele_balance_2":
                        b["allele_balance"],

                    "phase": "UNKNOWN",

                    "interpretation":
                        "CANDIDATE_COMPOUND_HET"
                })

    pairs.sort(
        key=lambda x: (
            x["level"],
            x["score"]
        ),
        reverse=True
    )

    return pairs


# ============================================================
# TRACK 1
# ============================================================

def build_track1(
    variants,
    pairs,
    limit=1000
):

    rows = []

    # Primero pares
    for p in pairs:

        rows.append({
            "proband_id": "WGS_EX2312012",

            "chrom_1": p["chrom_1"],
            "pos_1": p["pos_1"],
            "ref_1": p["ref_1"],
            "alt_1": p["alt_1"],

            "chrom_2": p["chrom_2"],
            "pos_2": p["pos_2"],
            "ref_2": p["ref_2"],
            "alt_2": p["alt_2"],

            "epcr": "",

            "finding_type":
                p["interpretation"],

            "notes":
                f'{p["gene_name"]}; '
                f'level={p["level"]}; '
                f'phase=UNKNOWN; '
                f'region1={p["region_1"]}; '
                f'region2={p["region_2"]}; '
                f'popAF1={p["population_af_1"]}; '
                f'popAF2={p["population_af_2"]}; '
                f'AB1={p["allele_balance_1"]}; '
                f'AB2={p["allele_balance_2"]}'
        })

        if len(rows) >= limit:
            return rows

    # Después variantes individuales
    sorted_vars = sorted(
        variants,
        key=lambda x: x["score"],
        reverse=True
    )

    for v in sorted_vars:

        rows.append({
            "proband_id": "WGS_EX2312012",

            "chrom_1": v["chrom"],
            "pos_1": v["pos"],
            "ref_1": v["ref"],
            "alt_1": v["alt"],

            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",

            "epcr": "",

            "finding_type":
                "SINGLE_VARIANT",

            "notes":
                f'{v["gene_name"]}; '
                f'region={v["region"]}; '
                f'GT={v["gt_class"]}; '
                f'AB={v["allele_balance"]}; '
                f'popAF={v["population_af"]}; '
                f'QUAL={v["qual"]}; '
                f'GQ={v["gq"]}; '
                f'score={v["score"]:.4f}'
        })

        if len(rows) >= limit:
            break

    return rows


# ============================================================
# TSV
# ============================================================

def write_tsv(path, rows):

    if not rows:
        return

    os.makedirs(
        os.path.dirname(path),
        exist_ok=True
    )

    fields = list(rows[0].keys())

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# ============================================================
# INSPECCIÓN
# ============================================================

def inspect_gene(
    gene_name,
    variants,
    pairs
):

    vals = [
        v for v in variants
        if v["gene_name"] == gene_name
    ]

    if not vals:
        print(
            f"\nCONTROL: {gene_name} no encontrado."
        )
        return

    print("\n")
    print("=" * 110)
    print(
        f"CONTROL POST-RANKING — {gene_name}"
    )
    print("=" * 110)

    hets = [
        v for v in vals
        if v["gt_class"] == "HET"
    ]

    homs = [
        v for v in vals
        if v["gt_class"] == "HOM_ALT"
    ]

    functional = [
        v for v in vals
        if v["functional_priority"] >= 3
    ]

    rare_known = [
        v for v in vals
        if (
            v["population_af"] is not None and
            v["population_af"] < 0.01
        )
    ]

    print(
        f"Variants:       {len(vals)}"
    )
    print(
        f"Hets:           {len(hets)}"
    )
    print(
        f"Hom-alt:        {len(homs)}"
    )
    print(
        f"Functional:     {len(functional)}"
    )
    print(
        f"Rare known AF:  {len(rare_known)}"
    )

    print("\nVARIANTES:")

    for i, v in enumerate(
        sorted(
            vals,
            key=lambda x: x["score"],
            reverse=True
        ),
        start=1
    ):

        print(
            f"{i:02d}. "
            f"{v['chrom']}:{v['pos']} "
            f"{v['ref']}>{v['alt']} "
            f"GT={v['gt_class']} "
            f"AB={v['allele_balance']} "
            f"popAF={v['population_af']} "
            f"DP={v['dp']} "
            f"GQ={v['gq']} "
            f"QUAL={v['qual']} "
            f"REGION={v['region']} "
            f"SCORE={v['score']:.4f}"
        )

    gene_pairs = [
        p for p in pairs
        if p["gene_name"] == gene_name
    ]

    print("\nPARES:")

    if not gene_pairs:
        print("    ninguno")
    else:
        for p in gene_pairs[:20]:

            print(
                f"    L{p['level']} "
                f"{p['pos_1']}:{p['ref_1']}>{p['alt_1']} "
                f"+ "
                f"{p['pos_2']}:{p['ref_2']}>{p['alt_2']} "
                f"score={p['score']:.4f} "
                f"phase={p['phase']}"
            )


# ============================================================
# MAIN
# ============================================================

def main():

    global GLOBAL_TRANSCRIPTS
    global GLOBAL_EXONS
    global GLOBAL_CDS

    os.makedirs(
        OUT_DIR,
        exist_ok=True
    )

    print("=" * 80)
    print(
        "TRINOVA.OS R7.1 — CHR15 "
        "VARIANT ARCHITECTURE + FUNCTIONAL TRIAGE"
    )
    print("=" * 80)

    print("\nVCF:")
    print(VCF_PATH)

    print("\nGTF:")
    print(GTF_PATH)

    print("\nCLINICAL:")
    print(CLINICAL_PATH)

    print("\nOUTPUT:")
    print(OUT_DIR)

    # --------------------------------------------------------
    # 1
    # --------------------------------------------------------

    print("\n[1/8] Perfil clínico...")

    hpos = load_clinical(
        CLINICAL_PATH
    )

    print(
        f"    HPO encontrados: {len(hpos)}"
    )

    for h in hpos:
        print(
            f"      {h}"
        )

    # --------------------------------------------------------
    # 2
    # --------------------------------------------------------

    print("\n[2/8] GENCODE v50 chr15...")

    (
        genes,
        transcripts,
        exons,
        cds
    ) = load_chr15_gencode(
        GTF_PATH
    )

    GLOBAL_TRANSCRIPTS = transcripts
    GLOBAL_EXONS = exons
    GLOBAL_CDS = cds

    print(
        f"    Genes chr15: {len(genes):,}"
    )

    print(
        f"    Transcripts: {len(transcripts):,}"
    )

    print(
        f"    Exon blocks: "
        f"{sum(len(x) for x in exons.values()):,}"
    )

    print(
        f"    CDS blocks: "
        f"{sum(len(x) for x in cds.values()):,}"
    )

    # --------------------------------------------------------
    # 3
    # --------------------------------------------------------

    print("\n[3/8] Índice genómico...")

    gene_index = build_gene_index(
        genes
    )

    print(
        f"    Bins: {len(gene_index)}"
    )

    # --------------------------------------------------------
    # 4
    # --------------------------------------------------------

    print("\n[4/8] VCF real exclusivamente chr15...")

    variants, sample_name = read_chr15_vcf(
        VCF_PATH,
        genes,
        gene_index
    )

    print(
        f"    Sample: {sample_name}"
    )

    print(
        f"    Variantes procesadas: "
        f"{len(variants):,}"
    )

    print(
        "    Separación explícita:"
    )

    print(
        "      allele_balance = AD/DP del individuo"
    )

    print(
        "      population_af  = INFO/AF solamente"
    )

    # --------------------------------------------------------
    # 5
    # --------------------------------------------------------

    print("\n[5/8] Arquitectura por gen...")

    gene_rows = aggregate_genes(
        variants
    )

    print(
        f"    Genes rankeados: "
        f"{len(gene_rows):,}"
    )

    # --------------------------------------------------------
    # 6
    # --------------------------------------------------------

    print("\n[6/8] Compound-het triage...")

    pairs = build_compound_pairs(
        variants
    )

    print(
        f"    Pares candidatos: "
        f"{len(pairs):,}"
    )

    levels = defaultdict(int)

    for p in pairs:
        levels[p["level"]] += 1

    print(
        f"    Level 1: {levels[1]:,}"
    )

    print(
        f"    Level 2: {levels[2]:,}"
    )

    print(
        f"    Level 3: {levels[3]:,}"
    )

    # --------------------------------------------------------
    # 7
    # --------------------------------------------------------

    print("\n[7/8] Ranking...")

    print("\n")
    print("=" * 115)
    print(
        "TOP 30 GENES — CHR15"
    )
    print("=" * 115)

    print(
        f"{'RANK':>4} "
        f"{'GENE':<24} "
        f"{'VAR':>5} "
        f"{'HET':>5} "
        f"{'FUNC':>5} "
        f"{'SPL':>5} "
        f"{'RARE':>5} "
        f"{'MAX':>8} "
        f"{'ARCH':>8}"
    )

    print("-" * 115)

    for row in gene_rows[:30]:

        print(
            f"{row['rank']:4d} "
            f"{row['gene_name'][:24]:<24} "
            f"{row['variants']:5d} "
            f"{row['hets']:5d} "
            f"{row['functional']:5d} "
            f"{row['splice']:5d} "
            f"{row['rare_known']:5d} "
            f"{row['max_variant_score']:8.3f} "
            f"{row['architecture_score']:8.3f}"
        )

    # --------------------------------------------------------
    # TOP PAIRS
    # --------------------------------------------------------

    print("\n")
    print("=" * 120)
    print(
        "TOP 30 PARES COMPOUND-HET TRIAGE"
    )
    print("=" * 120)

    for i, p in enumerate(
        pairs[:30],
        start=1
    ):

        print(
            f"{i:02d}. "
            f"L{p['level']} "
            f"{p['gene_name'][:20]:<20} "
            f"{p['pos_1']}:{p['ref_1']}>{p['alt_1']} "
            f"+ "
            f"{p['pos_2']}:{p['ref_2']}>{p['alt_2']} "
            f"score={p['score']:.3f} "
            f"phase={p['phase']} "
            f"regions="
            f"{p['region_1']}/{p['region_2']}"
        )

    # --------------------------------------------------------
    # 8
    # --------------------------------------------------------

    print("\n[8/8] Escribiendo resultados...")

    top_genes = gene_rows[:TOP_GENES]

    top_pairs = pairs[:TOP_PAIRS]

    track1 = build_track1(
        variants,
        pairs
    )

    write_tsv(
        os.path.join(
            OUT_DIR,
            "chr15_gene_ranking.tsv"
        ),
        top_genes
    )

    write_tsv(
        os.path.join(
            OUT_DIR,
            "chr15_compound_het_triage.tsv"
        ),
        top_pairs
    )

    write_tsv(
        os.path.join(
            OUT_DIR,
            "track1_predictions.csv"
        ),
        track1
    )

    # --------------------------------------------------------
    # Variantes informativas
    # --------------------------------------------------------

    informative = sorted(
        variants,
        key=lambda x: x["score"],
        reverse=True
    )[:1000]

    write_tsv(
        os.path.join(
            OUT_DIR,
            "chr15_top_variants.tsv"
        ),
        informative
    )

    # --------------------------------------------------------
    # Resumen JSON
    # --------------------------------------------------------

    summary = {
        "pipeline": "TRINOVA.OS R7.1",
        "chromosome": TARGET_CHROM,
        "sample": sample_name,

        "hpo_count": len(hpos),

        "genes_chr15": len(genes),
        "transcripts_chr15": len(transcripts),

        "variants_processed": len(variants),

        "genes_ranked": len(gene_rows),

        "compound_pairs": len(pairs),

        "compound_levels": {
            "level_1": levels[1],
            "level_2": levels[2],
            "level_3": levels[3]
        },

        "population_af_policy":
            "INFO/AF only; never inferred from AD/DP",

        "allele_balance_policy":
            "AD/DP within proband",

        "phase_policy":
            "UNKNOWN unless independently established",

        "phenotype_policy":
            "HPO loaded but no HPO-gene association database applied",

        "clinical_interpretation":
            "exploratory prioritization only"
    }

    with open(
        os.path.join(
            OUT_DIR,
            "summary.json"
        ),
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            summary,
            f,
            indent=2
        )

    # --------------------------------------------------------
    # CONTROL POST-RANKING
    # --------------------------------------------------------
    #
    # Solo inspección posterior.
    # NO participa en el ranking.
    # --------------------------------------------------------

    inspect_gene(
        "BUB1B",
        variants,
        pairs
    )

    print("\n")
    print("=" * 80)
    print("R7.1 READY")
    print("=" * 80)

    print(
        f"Output: {OUT_DIR}"
    )

    print(
        "\nArchivos:"
    )

    print(
        "  chr15_gene_ranking.tsv"
    )

    print(
        "  chr15_compound_het_triage.tsv"
    )

    print(
        "  chr15_top_variants.tsv"
    )

    print(
        "  track1_predictions.csv"
    )

    print(
        "  summary.json"
    )

    print(
        "\nIMPORTANTE:"
    )

    print(
        "  population_af NO se deriva de AD/DP."
    )

    print(
        "  allele_balance y population_af son métricas distintas."
    )

    print(
        "  phase permanece UNKNOWN."
    )

    print(
        "  BUB1B solo se inspecciona después del ranking."
    )

    print(
        "  HPO todavía no tiene asociación gen->fenotipo local."
    )

    print(
        "  El resultado es priorización computacional exploratoria."
    )


if __name__ == "__main__":
    main()
