import argparse
import os
import random
from typing import List

# Import Biopython modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Fixed parameters for the simulation
READ_LENGTH = 150
MIN_FRAGMENT_LENGTH = 400
MAX_FRAGMENT_LENGTH = 600
# Fixed high-quality score (Phred+33 score of 41, represented by 'J')
QUALITY_SCORES = 'J' * READ_LENGTH

def parse_fasta(filepath: str) -> List[SeqRecord]:
    """
    Uses Bio.SeqIO to parse a FASTA file and filters sequences that are too short 
    to produce a fragment of MAX_FRAGMENT_LENGTH + READ_LENGTH.
    """
    sequences = []
    # Minimum required length for a sequence to produce a fragment
    min_len = MAX_FRAGMENT_LENGTH + READ_LENGTH

    try:
        # SeqIO.parse returns an iterator of SeqRecord objects
        for record in SeqIO.parse(filepath, 'fasta'):
            # Ensure the sequence is in uppercase for consistency
            record.seq = record.seq.upper() 
            
            if len(record.seq) >= min_len:
                sequences.append(record)
            else:
                print(f"Warning: Ignored sequence '{record.id}' in {os.path.basename(filepath)} because it was shorter than the required {min_len}bp.")

    except FileNotFoundError:
        print(f"Error: Input file not found: {filepath}")
        return []
    except Exception as e:
        print(f"An error occurred while parsing {filepath}: {e}")
        return []

    return sequences

def main():
    """Main function to parse arguments and run the Biopython-based simulation."""
    parser = argparse.ArgumentParser(
        description="Simulate paired-end FASTQ reads from multiple FASTA files using Biopython."
    )
    parser.add_argument(
        '--fasta_files', 
        nargs='+', 
        help="One or more input FASTA files."
    )
    parser.add_argument(
        '--N', 
        type=int, 
        required=True, 
        help="N: Number of read pairs to generate per sequence in the FASTA file."
    )
    parser.add_argument(
        '--X', 
        type=int, 
        required=True, 
        help="X: Number of independent FASTQ file pairs to generate for each input FASTA file."
    )
    parser.add_argument(
        '-O', 
        '--output_dir', 
        type=str, 
        required=True, 
        help="The output directory where FASTQ files will be written."
    )

    args = parser.parse_args()

    # 1. Setup output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory set to: {os.path.abspath(output_dir)}")
    print("-" * 30)

    # 2. Process each FASTA file
    for fasta_file in args.fasta_files:
        print(f"Processing FASTA file: {os.path.basename(fasta_file)}")
        
        # Parse FASTA data using Biopython
        fasta_records = parse_fasta(fasta_file)
        
        if not fasta_records:
            print(f"Skipping {os.path.basename(fasta_file)}.")
            continue
            
        print(f"Found {len(fasta_records)} sequences suitable for simulation.")

        # Determine the base name for the output files
        input_file_name_base = os.path.splitext(os.path.basename(fasta_file))[0]
        out_dir = os.path.join(output_dir, input_file_name_base)
        try:
            os.mkdir(out_dir)
        except Exception:
            None
        output_base_path = os.path.join(out_dir, input_file_name_base)

        # 3. Generate X pairs of FASTQ files
        for i in range(args.X):
            output_r1 = f"{output_base_path}_pair_{i+1}_R1.fastq"
            output_r2 = f"{output_base_path}_pair_{i+1}_R2.fastq"
            
            print(f"  Generating pair {i + 1} of {args.X}...")

            read_count = 0
            
            try:
                # Open the paired FASTQ files for writing
                with open(output_r1, 'w') as fq1, open(output_r2, 'w') as fq2:
                    
                    for record in fasta_records:
                        header = record.id
                        seq = record.seq # Bio.Seq object
                        seq_len = len(seq)
                        
                        # Generate N reads per sequence
                        for j in range(args.N):
                            # A. Choose a random fragment length L (e.g., 400-600bp)
                            frag_len = random.randint(MIN_FRAGMENT_LENGTH, MAX_FRAGMENT_LENGTH)
                            
                            # B. Choose a random start position P
                            max_p = seq_len - frag_len # Ensure fragment ends within sequence bounds
                            start_pos = random.randint(0, max_p)
                            
                            # --- Read 1 (Forward strand, 5' end of fragment) ---
                            r1_seq_obj = seq[start_pos : start_pos + READ_LENGTH]
                            
                            # --- Read 2 (Reverse Complement of the 3' end of fragment) ---
                            # The start position of the R2 sequence in the forward strand
                            r2_start_pos = start_pos + frag_len - READ_LENGTH
                            r2_fragment_end_seq_obj = seq[r2_start_pos : r2_start_pos + READ_LENGTH]
                            
                            # Use Biopython's built-in reverse complement function (best practice!)
                            r2_seq_obj = r2_fragment_end_seq_obj.reverse_complement()
                            
                            # Create Read ID
                            # Format: @{FASTA_ID}:{N_index}_{X_index}_pos_{start}_{frag_len}/read_number
                            read_id_base = f"{header}:{j}_{i+1}_pos_{start_pos}_{frag_len}"

                            # Write Read 1 (R1)
                            fq1.write(f"@{read_id_base}/1\n")
                            fq1.write(f"{str(r1_seq_obj)}\n")
                            fq1.write("+\n")
                            fq1.write(f"{QUALITY_SCORES}\n")
                            
                            # Write Read 2 (R2)
                            fq2.write(f"@{read_id_base}/2\n")
                            fq2.write(f"{str(r2_seq_obj)}\n")
                            fq2.write("+\n")
                            fq2.write(f"{QUALITY_SCORES}\n")
                            
                            read_count += 1
                        
                print(f"    Wrote {read_count * 2} reads to {os.path.basename(output_r1)} and {os.path.basename(output_r2)}")

            except Exception as e:
                print(f"An error occurred while writing FASTQ pair {i+1}: {e}")
                
        print("-" * 30)

if __name__ == '__main__':
    print("NOTE: This script requires the Biopython library. Install with: pip install biopython")
    main()

