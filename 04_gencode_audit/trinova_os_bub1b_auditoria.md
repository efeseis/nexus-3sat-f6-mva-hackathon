# BUB1B — Variant, Reference, and Phase Audit

## 1. Objective

To record the results from `trinova_os_bub1b_auditoria_completa.py` and connect them with the Track 1 final summary table, the GENCODE v50 inspection, and the clinical hypothesis of the MVA case.

This document differentiates between the records observed in the VCF, the results calculated by the script, and the interpretations that require additional evidence.

## 2. Reference Context

The GENCODE v50 inspection identifies BUB1B via `ENSG00000156970.15`, at `chr15:40160984–40221137`, on GRCh38.

The analysis reports traversing 5,012,204 records from the sample's VCF (`WGS_PROBAND01`) and finding 14 variants within the analyzed interval.

## 3. Variants Matching the Final Table

### Variant 1: chr15:40209701 T>G

- Reported genotype: `0/1`.
- DP: 46.
- GQ: 99.
- Assigned region: CDS.
- REF–FASTA concordance: Yes.
- Reported consequence: `STOP_GAINED`.

### Variant 2: chr15:40220612 T>G

- Reported genotype: `0/1`.
- DP: 28.
- GQ: 99.
- Assigned region: CDS.
- REF–FASTA concordance: Yes.
- Reported consequence: `MISSENSE`.

Both variants match the positions and alleles of the primary row in the final summary table. The concordance of the REF alleles with the FASTA sequence supports this initial reference check.

Consequence labels are recorded as script outputs. The phase audit block presents the `CODON` and `AA` fields without values; therefore, codon and amino acid changes must be verified independently before considering them reproduced results.

## 4. REF–FASTA Check Results

The script detected four records that do not literally match the reference sequence:

| Position | VCF REF | VCF ALT | Base or sequence reported by FASTA |
|---|---|---|---|
| 40180642 | `CT` | `C` | `C` |
| 40181093 | `AT` | `A` | `A` |
| 40216568 | `TA` | `T` | `T` |
| 40216570 | `TA` | `T` | `T` |

All four records are indels. A literal discrepancy in an indel may depend on allele representation and normalization, but it should not be ignored: the complete sequence at the corresponding position must be verified to ensure that both the VCF and the FASTA utilize the same reference and coordinate convention.

These four discrepancies do not correspond to the two CDS variants in the final table, which the script marks as concordant with FASTA.

## 5. Phase Audit Results

Both CDS variants present a `0/1` genotype. However, the `PGT` and `PID` fields appear as `.` for both.

Consequently, the audit classifies the phase as undetermined with the available information.

It cannot be confirmed that the variants are in *trans*, nor can it be stated that they constitute a biallelic compound heterozygous variant pair solely from the `0/1` genotypes.

The phase requires additional informative evidence, such as reads that allow establishing the relationship between alleles or adequate phasing data.

## 6. Relationship to the MVA Hypothesis

BUB1B remains a biologically reasonable candidate for investigation due to BUBR1's function in controlling chromosome segregation.

The reported presence of two heterozygous variants in the CDS region makes studying their potential relevance pertinent. Nevertheless, the molecular consequence must be validated independently, and the phase remains unresolved.

The available data do not demonstrate pathogenicity, causality, or a biallelic mechanism responsible for the phenotype on their own.

## 7. Conclusion and Status

This audit strengthens the traceability between the final table and the VCF records for the two BUB1B positions. It also identifies four REF–FASTA discrepancies in other records and documents that the phase of the two candidate variants is undetermined.

**Status:** Reported positions and genotypes match; REF is concordant with FASTA for both CDS variants; protein consequences pending independent reproduction; phase undetermined; clinical causality not established.

The evaluated final table is preserved without modifications. This document constitutes complementary audit and reproducibility evidence.
