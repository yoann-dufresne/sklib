# Benchmark datasets

Two ways to feed a dataset to the harness:

1. **Catalogued** (single download, wired in `benchmark/scripts/lib.sh`): just name it in
   `DATASETS=...`. Available: `sarscov2`, `ecoli` (local fixtures); `yeast`, `celegans`,
   `chr21`, `chr20`, `chr1`, `chm13` (downloaded from UCSC on first use).

2. **Custom / multi-file** (pangenomes, metagenomes, read sets): build a single sanitized
   FASTA and drop it where the harness caches genomes — `prepare_genome` then reuses it
   verbatim, no code change:
   ```bash
   GEN=benchmark/data/genomes
   python3 benchmark/scripts/e2e_helpers.py sanitize my_input.fa > "$GEN/<name>.sanitized.fa"
   DATASETS="<name> ecoli" bash benchmark/scripts/construct.sh
   ```
   ("sanitize" = uppercase + split at non-ACGT runs, so sklib and KMC agree.)

The single-genome ladder below fits in 62 GB RAM. Pangenomes/deep metagenomes exceed it —
run those on a larger-RAM machine (the harness is portable; CSVs concatenate).

## Scaling ladder (catalogued, fits in RAM)
| name | organism | ~size | URL source |
|---|---|---|---|
| sarscov2 | SARS-CoV-2 | 30 kb | local fixture |
| ecoli | E. coli | 4.6 Mb | local fixture |
| yeast | S. cerevisiae sacCer3 | 12 Mb | UCSC |
| celegans | C. elegans ce11 | 100 Mb | UCSC |
| chr21 / chr1 | human chr | 40 / 230 Mb | UCSC hg38 |
| chm13 | human T2T-CHM13v2 | 3.1 Gb | UCSC hs1 |

## Bacterial pangenome (high redundancy → best-case compaction)
Install the NCBI `datasets` CLI (no sudo):
```bash
curl -sL 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' \
  -o ~/.local/bin/datasets && chmod +x ~/.local/bin/datasets
```
Download N complete E. coli assemblies, concatenate, register:
```bash
datasets download genome taxon "Escherichia coli" --assembly-level complete \
  --limit 200 --include genome --filename ecoli_pan.zip
unzip -p ecoli_pan.zip 'ncbi_dataset/data/*/*.fna' > ecoli_pan.fa
python3 benchmark/scripts/e2e_helpers.py sanitize ecoli_pan.fa > benchmark/data/genomes/ecoli_pan.sanitized.fa
```
Then `DATASETS="ecoli_pan"`. Swap taxon for `"Salmonella enterica"` for a second pangenome.

## Rice (large plant genome / pangenome)
Single reference IRGSP-1.0 Nipponbare via NCBI (RefSeq `GCF_001433935.1`, ~380 Mb):
```bash
datasets download genome accession GCF_001433935.1 --include genome --filename rice.zip
unzip -p rice.zip 'ncbi_dataset/data/*/*.fna' > rice.fa
python3 benchmark/scripts/e2e_helpers.py sanitize rice.fa > benchmark/data/genomes/rice.sanitized.fa
```
For a rice *pangenome*, fetch several accessions from the 3K-RGP / NCBI and concatenate as above.

## Metagenome (high diversity, many singletons → worst-case compaction)
Defined mock — ZymoBIOMICS, or any assembled metagenome. Example: assembled contigs of a
gut WGS sample. To use raw reads as the indexed set, sanitize the reads FASTA the same way.
For a controlled diversity sweep, mix M bacterial genomes from the pangenome step.

## Real reads (query workload, mixed status)
Instead of `simreads`, query a real read set of the indexed organism. ENA gives direct
FASTQ URLs (no SRA toolkit); convert to FASTA and point the query at it:
```bash
# e.g. an E. coli Illumina run (replace the accession)
curl -sL 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR.../<run>_1.fastq.gz' \
  | zcat | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' > reads.fa
# then query an existing index directly:
build/bin/sskm query -l <index>.sskm -i reads.fa -o /dev/null
```
(The harness `reads` workload uses simulated reads; real reads are an add-on per organism.)
