import sys
import json
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

# --- Configuration ---
# Columns defined in the request:
# 0-qseqid, 1-sseqid, 2-pident, 3-length, 4-mismatch, 5-gapopen, 6-qstart, 7-qend, 8-sstart, 9-send, 10-evalue, 11-bitscore
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

def parse_fasta(filepath):
    """Reads a FASTA file using SeqIO and returns a dictionary of {sequence_id: length}."""
    lengths = {}
    try:
        # Use SeqIO.parse to read FASTA records efficiently
        for record in SeqIO.parse(filepath, "fasta"):
            # SeqIO record.id automatically handles the part before the first space
            lengths[record.id] = len(record.seq)
        return lengths
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file {filepath}: {e}", file=sys.stderr)
        print("Hint: If 'Bio' module is missing, install Biopython ('pip install biopython').", file=sys.stderr)
        sys.exit(1)

def parse_blastn(filepath):
    """Reads the BLASTn TSV file using pandas and returns a list of hit dictionaries."""
    try:
        # Load the TSV file using pandas
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=BLAST_COLUMNS,
            # Use 'str' to prevent pandas from inferring types for ID columns
            dtype={'qseqid': str, 'sseqid': str}
        )

        # Convert coordinates and bitscore to appropriate types
        df['sstart'] = df['sstart'].astype(int)
        df['send'] = df['send'].astype(int)
        # Use pd.to_numeric with errors='coerce' to turn non-numeric values into NaN
        df['bitscore'] = pd.to_numeric(df['bitscore'], errors='coerce')

        # Drop rows where bitscore conversion failed (e.g., headers or non-numeric scores/coordinates)
        df = df.dropna(subset=['bitscore'])

        # Convert DataFrame rows back to a list of dictionaries for compatibility
        hits = df.to_dict('records')
        return hits
    except FileNotFoundError:
        print(f"Error: BLAST file not found at {filepath}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Warning: BLAST file {filepath} is empty.", file=sys.stderr)
        return []
    except Exception as e:
        print(f"Error reading BLAST file {filepath}: {e}", file=sys.stderr)
        print("Hint: If 'pandas' module is missing, install it ('pip install pandas').", file=sys.stderr)
        sys.exit(1)

def intervals_overlap(a_start, a_end, b_start, b_end):
    """
    Checks if two intervals [a_start, a_end] and [b_start, b_end] overlap by more than 10%
    of the length of the shorter interval.
    """
    # 1. Normalize intervals (ensure start <= end)
    a_min, a_max = min(a_start, a_end), max(a_start, a_end)
    b_min, b_max = min(b_start, b_end), max(b_start, b_end)

    # 2. Calculate the intersection (overlap)
    overlap_start = max(a_min, b_min)
    overlap_end = min(a_max, b_max)

    overlap_length = overlap_end - overlap_start

    # If overlap_length is 0 or negative, there is no physical intersection
    if overlap_length <= 0:
        return False

    # 3. Calculate lengths of the two intervals
    len_A = a_max - a_min
    len_B = b_max - b_min

    # Should not happen with typical BLAST results, but ensures robustness against zero-length intervals
    if len_A <= 0 or len_B <= 0:
        return False

    # 4. Determine the minimum length and allowed tolerance (10%)
    min_len = min(len_A, len_B)
    allowed_overlap = min_len * 0.10

    # 5. Overlap is considered 'true' if the intersection is greater than the allowed tolerance
    return overlap_length > allowed_overlap

def filter_overlapping_hits(all_hits):
    """
    Greedy filtering: Keeps the hit with the best bitscore for subject regions,
    discarding all lower-scoring hits that overlap on the subject sequence by more
    than the allowed tolerance.
    """
    # 1. Sort all hits by bitscore (descending)
    sorted_hits = sorted(all_hits, key=lambda x: x['bitscore'], reverse=True)

    best_hits = []
    # Tracks the subject regions already covered by a higher-scoring hit
    covered_subject_regions = defaultdict(list)  # {sseqid: [(start, end), ...]}

    for hit in sorted_hits:
        sseqid = hit['sseqid']
        s_start = hit['sstart']
        s_end = hit['send']

        # Normalize subject coordinates (start <= end)
        norm_start = min(s_start, s_end)
        norm_end = max(s_start, s_end)

        is_overlapping = False
        # Check against already covered regions for the same subject sequence
        for covered_start, covered_end in covered_subject_regions[sseqid]:
            # Use the updated overlap check
            if intervals_overlap(norm_start, norm_end, covered_start, covered_end):
                is_overlapping = True
                break

        if not is_overlapping:
            # This hit has the highest bitscore in this region, so we keep it
            best_hits.append(hit)
            # Mark this region as covered
            covered_subject_regions[sseqid].append((norm_start, norm_end))

    return best_hits

def merge_intervals(intervals):
    """
    Merges a list of [start, end] intervals into a minimal set of non-overlapping intervals
    and calculates the total length of the covered regions.
    """
    if not intervals:
        return 0

    # 1. Sort by start coordinate
    intervals.sort(key=lambda x: x[0])

    merged = []
    current_start, current_end = intervals[0]

    for next_start, next_end in intervals[1:]:
        # If the next interval starts before or at the current interval's end, they overlap
        if next_start <= current_end:
            # Extend the current interval to the maximum of the two ends
            current_end = max(current_end, next_end)
        else:
            # No overlap, save the current merged interval and start a new one
            merged.append([current_start, current_end])
            current_start, current_end = next_start, next_end

    # Add the last merged interval
    merged.append([current_start, current_end])

    # Calculate total length (sum of (end - start))
    total_length = sum(end - start for start, end in merged)
    return total_length

def calculate_stats(query_lengths, subject_lengths, all_hits, best_hits):
    """
    Calculates the required statistics based on the total query and subject lengths,
    the total raw hits, and the list of non-overlapping best hits.
    """
    # --- Step 4 & 5: Query Contig Analysis ---
    aligned_qseqids = set(hit['qseqid'] for hit in best_hits)
    all_qseqids = set(query_lengths.keys())
    unaligned_qseqids = all_qseqids - aligned_qseqids

    # 4. Count of aligned/unaligned query contigs
    num_aligned_contigs = len(aligned_qseqids)
    num_unaligned_contigs = len(unaligned_qseqids)

    # 5. Total length of aligned/unaligned contigs
    total_aligned_length = sum(query_lengths.get(qid, 0) for qid in aligned_qseqids)
    total_unaligned_length = sum(query_lengths.get(qid, 0) for qid in unaligned_qseqids)
    total_query_length = sum(query_lengths.values())

    # --- Subject Coverage Calculation ---
    subject_coverage_intervals = defaultdict(list)
    for hit in best_hits:
        s_start = hit['sstart']
        s_end = hit['send']
        # Normalize and add the interval [min, max]
        subject_coverage_intervals[hit['sseqid']].append([min(s_start, s_end), max(s_start, s_end)])

    total_aligned_subject_length = 0
    for sseqid, intervals in subject_coverage_intervals.items():
        # Merge intervals for each subject sequence to get non-redundant coverage
        total_aligned_subject_length += merge_intervals(intervals)

    # --- Step 5: Subject Analysis ---
    total_subject_length = sum(subject_lengths.values())
    
    percent_aligned_subject = 0
    if total_subject_length > 0:
        percent_aligned_subject = (total_aligned_subject_length / total_subject_length) * 100

    # --- Step 3: Best Hit Stats ---
    max_bitscore = max(hit['bitscore'] for hit in best_hits) if best_hits else 0
    num_non_overlapping_best_hits = len(best_hits)

    results = {
        "query_contig_statistics": {
            "total_contigs": len(all_qseqids),
            "num_aligned_contigs": num_aligned_contigs,
            "num_unaligned_contigs": num_unaligned_contigs,
            "total_aligned_length": total_aligned_length,
            "total_unaligned_length": total_unaligned_length,
            "total_query_length": total_query_length
        },
        "subject_statistics": {
            "total_subject_length": total_subject_length,
            "total_aligned_subject_length": total_aligned_subject_length,
            "percent_subject_aligned": round(percent_aligned_subject, 2)
        },
        "blast_hit_statistics": {
            "total_raw_hits": len(all_hits), # Total hits before filtering
            "num_non_overlapping_best_hits": num_non_overlapping_best_hits, # Hits after filtering
            "max_best_hit_bitscore": max_bitscore
        }
    }
    return results

def write_json(filepath, data):
    """Writes the results dictionary to a JSON file."""
    try:
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"\nSuccessfully wrote analysis results to {filepath}")
    except Exception as e:
        print(f"Error writing JSON file {filepath}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) != 5:
        print("Usage: python blast_analysis.py <query_contigs.fasta> <subject.fasta> <blastn_results.tsv> <output.json>", file=sys.stderr)
        sys.exit(1)

    query_file = sys.argv[1]
    subject_file = sys.argv[2]
    blast_file = sys.argv[3]
    output_file = sys.argv[4]

    print(f"Parsing query contigs from: {query_file} (using Bio.SeqIO)...")
    query_lengths = parse_fasta(query_file)

    print(f"Parsing subject sequences from: {subject_file} (using Bio.SeqIO)...")
    subject_lengths = parse_fasta(subject_file)

    print(f"Parsing BLASTn results from: {blast_file} (using pandas)...")
    all_hits = parse_blastn(blast_file)

    if not all_hits:
        print("Warning: No valid BLAST hits found. Generating statistics based on zero hits.")

    print("Filtering overlapping subject hits (keeping best bitscore) with 10% overlap tolerance...")
    best_hits = filter_overlapping_hits(all_hits)

    print("Calculating final statistics...")
    results = calculate_stats(query_lengths, subject_lengths, all_hits, best_hits)

    write_json(output_file, results)
    print("Analysis complete.")


if __name__ == "__main__":
    main()
