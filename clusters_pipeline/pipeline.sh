#!/bin/bash

# Configuration
DIR1="../../ASSEMBLIES/DSI_10/"
DIR2="../../ASSEMBLIES/DSII_10_L20/"
REFERENCE="../../references/giant_virus_1.fa"
THREADS=8
GENOME_LEN=106365
OUT_DIR="../../OUTPUTS/CLUSTERS/SIZE_10_L20"
TMP_DIR="./tmp"
CONDA_ENV="mmseq2"

# 1 & 2. Setup Environment
mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"
ALL_CONTIGS="$OUT_DIR/all_contigs.fasta"
MMSEQ_CLUSTER_PRE="$OUT_DIR/mmseq_clusters"
MMSEQ_SEARCH_OUT="$OUT_DIR/mmseq_search_results.tsv"
CLUSTERS_SUBDIR="$OUT_DIR/CLUSTERS"

# 3. Concatenate Fasta Files
echo "Concatenating fasta files..."
cat "$DIR1"/*.fa "$DIR1"/*.fasta "$DIR2"/*.fa "$DIR2"/*.fasta > "$ALL_CONTIGS" 2>/dev/null

# 4. Reference Search
echo "Searching against reference..."
conda run -n "$CONDA_ENV" mmseqs easy-search "$ALL_CONTIGS" "$REFERENCE" "$MMSEQ_SEARCH_OUT" "$TMP_DIR" \
    --threads "$THREADS" \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen" \
    --search-type 3

# 5. MMseqs2 Linclust
echo "Clustering sequences..."
conda run -n "$CONDA_ENV" mmseqs easy-linclust "$ALL_CONTIGS" "$MMSEQ_CLUSTER_PRE" "$TMP_DIR" \
    --threads "$THREADS" \
    -c 0.5 \
    --min-seq-id 0.95 \
    --cov-mode 1 \
    --cluster-mode 2

# 6. Extract Cluster Fastas
# Added explicit -t flag for the TSV input.
echo "Extracting cluster fastas..."
python extact_mmseq_clusters.py \
    -i "$ALL_CONTIGS" \
    -t "${MMSEQ_CLUSTER_PRE}_cluster.tsv" \
    -m 3 \
    -o "$CLUSTERS_SUBDIR"

# 7. Statistical Analysis
# Fixed --clusters_dir to point to the actual directory where fastas are located.
echo "Analyzing clusters..."
python analyse_clusters.py \
    --clusters_dir "$CLUSTERS_SUBDIR" \
    --before_dir "$DIR1/*.fa" \
    --after_dir "$DIR2/*.fa" \
    --mmseq_tsv "$MMSEQ_SEARCH_OUT" \
    --all_contigs "$ALL_CONTIGS" \
    --output_json "$OUT_DIR/stats.json" \
    --genome_len "$GENOME_LEN"

echo "Pipeline complete. Results in $OUT_DIR"
