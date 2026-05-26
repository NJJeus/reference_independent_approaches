import os
import glob
import json
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from scipy.stats import fisher_exact

def clean_name(name):
    """Maintained from original script logic."""
    return name

def calculate_unique_coverage(df, start_col, end_col):
    """Calculates non-overlapping alignment coverage."""
    if df.empty: return 0
    intervals = []
    for _, row in df.iterrows():
        s, e = row[start_col], row[end_col]
        intervals.append((min(s, e), max(s, e)))
    intervals.sort()
    
    merged_length = 0
    if not intervals: return 0
    curr_start, curr_end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= curr_end:
            curr_end = max(curr_end, next_end)
        else:
            merged_length += (curr_end - curr_start + 1)
            curr_start, curr_end = next_start, next_end
    merged_length += (curr_end - curr_start + 1)
    return merged_length

def main():
    parser = argparse.ArgumentParser(description="Process viral clusters and calculate statistics.")
    parser.add_argument("--clusters_dir", default="./CLUSTERS_L_20/", help="Path to cluster fasta files")
    parser.add_argument("--before_dir", default="./DSI_100/*.fa", help="Glob pattern for 'Before' samples")
    parser.add_argument("--after_dir", default="./DSII_100_L20/*.fa", help="Glob pattern for 'After' samples")
    parser.add_argument("--mmseq_tsv", default="./mmseq2_search_l20.tsv", help="MMseq2 results tsv")
    parser.add_argument("--all_contigs", default="./all_contigs_l_20.fa", help="Path to all_contigs fasta")
    parser.add_argument("--output_json", default="stats.json", help="Output JSON file name")
    parser.add_argument("--genome_len", type=int, default=106365, help="Reference genome length for coverage")
    
    args = parser.parse_args()

    # 1. Read Cluster Data
    clusters_data = {}
    for filename in os.listdir(args.clusters_dir):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            cluster_id = filename.rsplit('.', 1)[0]
            sequences = []
            sample_names = set()
            total_length = []
            
            for record in SeqIO.parse(os.path.join(args.clusters_dir, filename), "fasta"):
                s_length = len(record.seq)
                if s_length < 1000:
                    continue
                total_length.append(s_length)
                
                parts = record.id.replace('VIROME5__M', 'M_').split('_')
                sample_name = f"{parts[0]}"
                
                sequences.append({
                    "seq_id": record.id,
                    "length": s_length,
                    "sample": sample_name
                })
                sample_names.add(sample_name)

            if len(sequences) < 4:
                continue
                
            clusters_data[cluster_id] = {
                "sequences": sequences,
                "samples": sample_names,
                "Mean_Length": sum(total_length) / len(sequences),
                "Max_Length": max(total_length)
            }

    # 2. Before/After Comparison & Fisher Test
    before = set([os.path.basename(i).split('.')[0] for i in glob.glob(args.before_dir)])
    after = set([os.path.basename(i).split('.')[0] for i in glob.glob(args.after_dir)])
    l_b, l_a = len(before), len(after)

    for k in clusters_data.keys():
        b = len(clusters_data[k]['samples'] & before)
        a = len(clusters_data[k]['samples'] & after)
        clusters_data[k]['Before_After'] = (b, a)
        fisher = fisher_exact([[b, a], [l_b-b, l_a-a]], alternative='less')
        clusters_data[k]['Fisher_Effect'] = fisher[0]
        clusters_data[k]['p_value'] = fisher[1]

    # 3. Filter Significant Clusters
    alpha = 0.05 / len(clusters_data) if clusters_data else 0.05
    good_clusters = {k: v for k, v in clusters_data.items() if v['p_value'] < alpha}
    
    for k, v in good_clusters.items():
        v['Total_Length'] = sum([i['length'] for i in v['sequences']])

    # 4. MMseq2 and Alignment Logic
    mmseq_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".split(' ')
    mmseq_results = pd.read_table(args.mmseq_tsv, header=None, names=mmseq_cols)
    positive_seq = sum([[seq['seq_id'] for seq in clus['sequences']] for clus in good_clusters.values()], [])

    all_seq = [[rec.id, len(rec.seq)] for rec in SeqIO.parse(args.all_contigs, 'fasta')]
    all_seq_df = pd.DataFrame(all_seq, columns=['SeqID', 'Length'])

    blast_df = mmseq_results.query('pident > 90 & length > 1000').drop_duplicates('qseqid').copy()
    blast_df['qseqid_clean'] = blast_df['qseqid'].apply(clean_name)

    true_positive_df = blast_df[blast_df['qseqid_clean'].isin(positive_seq)]
    
    # Coverage Calculation
    tp_alignment = true_positive_df.groupby('sseqid').apply(
        lambda x: calculate_unique_coverage(x, 'sstart', 'send')
    ).sum()

    # Length-based Stats
    tp_seqs = true_positive_df.qseqid.to_list()
    tp_total_len = all_seq_df.query("SeqID in @tp_seqs").Length.sum()
    
    fp_contigs = [i for i in positive_seq if i not in blast_df.qseqid.to_list()]
    fp_total_len = all_seq_df.query("SeqID in @fp_contigs").Length.sum()

    # 5. Compile Statistics
    stats = {
        "true_positive_alignment_abs": int(tp_alignment),
        "coverage_ratio": float(tp_alignment / args.genome_len),
        "precision_by_length": float(tp_total_len / (tp_total_len + fp_total_len)) if (tp_total_len + fp_total_len) > 0 else 0,
        "significant_clusters_count": len(good_clusters),
        "total_clusters_processed": len(clusters_data)
    }

    with open(args.output_json, 'w') as f:
        json.dump(stats, f, indent=4)
    
    print(f"Statistics written to {args.output_json}")

if __name__ == "__main__":
    main()
