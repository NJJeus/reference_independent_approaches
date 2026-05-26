import pandas as pd
import numpy as np
import argparse
import sys

# Define the translation table globally as it is a constant
utranslation_table = str.maketrans('0123', 'ATGC')

def int_to_k_mer(i, k_size):
    """
    Converts an integer representation of a k-mer (base 4) into its
    DNA sequence string (A, T, G, C) and pads it with 'A's to the k_size.
    """
    global utranslation_table
    
    # Convert integer 'i' to base 4 string
    seq = np.base_repr(i, 4).translate(utranslation_table)
    
    # Pad with 'A's to the required k_size
    return ''.join(['A'] * (k_size - len(seq))) + seq

def main():
    # 1. Setup argparse
    parser = argparse.ArgumentParser(
        description="Convert integer-encoded K-mers from a CSV file to DNA sequences and write them to a new file."
    )
    
    parser.add_argument(
        '--input', 
        type=str, 
        required=True, 
        help='Input CSV file containing K-mer integer IDs (e.g., p_values3_no_bonferoni.csv).'
    )
    
    parser.add_argument(
        '--output', 
        type=str, 
        required=True, 
        help='Output file to write the K-mer sequences in FASTA-like format (e.g., k_mers_no_bonferoni).'
    )
    
    parser.add_argument(
        '--k_size', 
        type=int, 
        default=22, 
        help='The length of the K-mer (used for padding). Default is 22.'
    )
    
    args = parser.parse_args()

    # 2. Validation (Optional but recommended)
    if args.k_size <= 0:
        print("Error: k_size must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    # 3. Read the input CSV
    try:
        # Assuming the K-mer integer IDs are in the first column (index 0)
        filtered_k_mers = pd.read_csv(
            args.input, 
            header=None, 
            usecols=[0]
        )[0].to_list()
    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Process and write the output
    print(f"Processing {len(filtered_k_mers)} K-mers (size {args.k_size})...")
    
    try:
        # Open the output file in write mode ('w') - or append ('a') if that's preferred
        with open(args.output, 'w') as handle:
            for k_int in filtered_k_mers:
                k_seq = int_to_k_mer(k_int, args.k_size)
                # FASTA-like format: >ID\nSEQUENCE\n
                to_write = f">{k_int}\n{k_seq}\n"
                handle.write(to_write)
        print(f"Successfully wrote results to '{args.output}'")
        
    except Exception as e:
        print(f"Error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
