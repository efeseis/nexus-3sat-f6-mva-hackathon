= RESTART: C:/Users/ericr/AppData/Local/Programs/Python/Python313/ADN_ estudio/ff/validador independiente_fix.py
============================================================================
BioNexus — VALIDACIÓN TRACK 1
Reproducción histórica + auditoría independiente
============================================================================

============================================================================
1. VERIFICACIÓN DE ARCHIVOS
============================================================================
[OK] clinico.json: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\clinico.json
[OK] GENCODE v50 GTF: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\gencode.v50.annotation.gtf.gz
[OK] FASTA GRCh38: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\GRCh38.primary_assembly.genome.fa.gz
[OK] VCF WGS: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\WGS_proband01.vcf.gz
[OK] índice VCF TBI: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\proband01.vcf.gz.tbi

============================================================================
2. VALIDACIÓN DEL JSON CLÍNICO
============================================================================
JSON sintácticamente válido: sí
Tipo raíz: dict
Claves raíz: ['rare_disease_genomic_challenge']
Claves del desafío: ['document_type', 'confidentiality', 'clinical_features']

============================================================================
3. ANOTACIÓN GENCODE V50: BUB1B
============================================================================
Líneas GTF no comentadas: 11,242,714
Conteos específicos de BUB1B: {'gene': 1, 'transcript': 31, 'exon': 511, 'CDS': 366, 'start_codon': 23, 'stop_codon': 25, 'UTR': 166}
Gen: BUB1B | ID: ENSG00000156970.15
Intervalo: chr15:40160984-40221137
Hebbra: +
Longitud inclusiva: 60,154 bp

============================================================================
4. VCF WGS Y MÉTRICAS DEL MOTOR HISTÓRICO
============================================================================
Muestra: WGS_proband01
Registros VCF procesados: 5,012,204
SNV de una base elegibles: 3,973,868
Registros no elegibles/omitidos: 1,038,336
Variantes encontradas dentro de BUB1B: 14
Ningún cromosoma necesitó fallback sintético.

TOP 10 — COHERENCIA
 1. chr 7 | C=0.6408 | trits=300
 2. chr12 | C=0.6201 | trits=300
 3. chr16 | C=0.6090 | trits=300
 4. chr14 | C=0.6077 | trits=300
 5. chr11 | C=0.6063 | trits=300
 6. chr 1 | C=0.6035 | trits=300
 7. chr 8 | C=0.6021 | trits=300
 8. chr 4 | C=0.6008 | trits=300
 9. chr 5 | C=0.6008 | trits=300
10. chr18 | C=0.6008 | trits=300

TOP 10 — ENTROPÍA
 1. chr15 | S=1.5354 | trits=300
 2. chr21 | S=1.5224 | trits=300
 3. chr22 | S=1.5205 | trits=300
 4. chr 3 | S=1.5106 | trits=300
 5. chr X | S=1.5025 | trits=300
 6. chr 5 | S=1.5020 | trits=300
 7. chr13 | S=1.4969 | trits=300
 8. chr20 | S=1.4969 | trits=300
 9. chr 2 | S=1.4936 | trits=300
10. chr Y | S=1.4896 | trits=300

COMPARACIÓN CON REFERENCIAS HISTÓRICAS
coherence chr 7: histórico=0.6408 | actual=0.6408 | delta=0.0
coherence chr12: histórico=0.6201 | actual=0.6201 | delta=0.0
coherence chr16: histórico=0.6090 | actual=0.609 | delta=0.0
entropy   chr15: histórico=1.5354 | actual=1.5354 | delta=0.0
entropy   chr21: histórico=1.5224 | actual=1.5224 | delta=0.0
entropy   chr22: histórico=1.5205 | actual=1.5205 | delta=0.0

PAR PRINCIPAL BUB1B DEL CSV TRACK 1
chr15:40209701 T>G | ENCONTRADA | GT=0/1 DP=46 GQ=99 PS= PID= PGT=
chr15:40220612 T>G | ENCONTRADA | GT=0/1 DP=28 GQ=99 PS= PID= PGT=
Fase: Fase no establecida por los campos GT/PS/PID/PGT disponibles para el par.

============================================================================
5. VALIDACIÓN DE REF CONTRA FASTA GRCh38
============================================================================
chr15:40168264 G>A | FASTA=G | MATCH
chr15:40180642 CT>C | FASTA=CT | MATCH
chr15:40181093 AT>A | FASTA=AT | MATCH
chr15:40182809 A>C | FASTA=A | MATCH
chr15:40186692 G>A | FASTA=G | MATCH
chr15:40187404 A>G | FASTA=A | MATCH
chr15:40192892 C>T | FASTA=C | MATCH
chr15:40193435 T>C | FASTA=T | MATCH
chr15:40198239 G>GT | FASTA=G | MATCH
chr15:40209701 T>G | FASTA=T | MATCH
chr15:40216470 A>G | FASTA=A | MATCH
chr15:40216568 TA>T | FASTA=TA | MATCH
chr15:40216570 TA>T | FASTA=TA | MATCH
chr15:40220612 T>G | FASTA=T | MATCH
Líneas FASTA leídas: 38,977,868
Tiempo de lectura FASTA: 67.6 s

============================================================================
VALIDACIÓN FINALIZADA
Estado: COMPLETED_WITH_REVIEW_REQUIRED
Carpeta de resultados: C:\Users\ericr\AppData\Local\Programs\Python\Python313\MVA-DATA\mva_data\validacion_track1_bionexus
CSV original: no modificado.
============================================================================
