import argparse
import os
import sys
import numpy as np
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description="Split FASTA by MMseqs2 clusters with stats.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-t", "--tsv", required=True, help="MMseqs2 cluster TSV file")
    parser.add_argument("-o", "--outdir", default="clusters_output", help="Output directory")
    parser.add_argument("-m", "--min_size", type=int, default=1, help="Minimum sequences per cluster")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Map Sequence IDs to Cluster IDs
    print(f"[*] Reading cluster mapping from {args.tsv}...")
    cluster_map = {}
    cluster_counts = defaultdict(int)

    with open(args.tsv, 'r') as f:
        for line in f:
            c_id, seq_id = line.strip().split('\t')
            cluster_map[seq_id] = c_id
            cluster_counts[c_id] += 1

    # 2. Parse FASTA and group sequences
    print(f"[*] Parsing FASTA and calculating lengths...")
    grouped_seqs = defaultdict(list)
    centroid_lengths = {}

    with open(args.input, 'r') as handle:
        for title, seq in tqdm(SimpleFastaParser(handle), desc="Processing Seqs"):
            sid = title.split()[0]
            cid = cluster_map.get(sid)

            if cid:
                # Store sequence
                grouped_seqs[cid].append(f">{title}\n{seq}")
                # If this sequence is the representative (centroid), store its length
                if sid == cid:
                    centroid_lengths[cid] = len(seq)

    # 3. Filter and Statistics
    final_clusters = {cid: seqs for cid, seqs in grouped_seqs.items() if len(seqs) >= args.min_size}
    sizes = [len(s) for s in final_clusters.values()]
    c_lens = [centroid_lengths[cid] for cid in final_clusters.keys() if cid in centroid_lengths]

    print("\n" + "="*30)
    print("CLUSTER SIZE STATISTICS (Quantiles)")
    print("="*30)
    if sizes:
        for q in range(0, 101, 10):
            print(f"{q}% quantile: {np.percentile(sizes, q):.1f} seqs")

    print("\n" + "="*30)
    print("CENTROID LENGTH STATISTICS")
    print("="*30)
    if c_lens:
        print(f"Mean Length:   {np.mean(c_lens):.2f} bp")
        print(f"Median Length: {np.median(c_lens):.2f} bp")
        print(f"Max Length:    {np.max(c_lens)} bp")
        print(f"Min Length:    {np.min(c_lens)} bp")
    print("="*30 + "\n")

    # 4. Writing files
    print(f"[*] Writing {len(final_clusters)} cluster files to {args.outdir}...")
    for cid, seq_list in tqdm(final_clusters.items(), desc="Writing Files"):
        count = len(seq_list)
        # Added count to filename
        out_path = os.path.join(args.outdir, f"{cid}_size_{count}.fa")

        with open(out_path, 'w') as f:
            f.write("\n".join(seq_list))

    print("[+] Done.")

if __name__ == "__main__":
    main()

