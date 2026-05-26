from functools import lru_cache
import pandas as pd
import tqdm
import numpy as np
from scipy.stats import fisher_exact
import argparse
import sys # Import sys for checking if we are using an interactive environment

@lru_cache(maxsize=None)
def perform_fisher(before: int, after: int, TOTAL_SAMPLES1: int, TOTAL_SAMPLES2: int):
	"""
	Performs Fisher's exact test (one-sided, 'greater') on the abundance data.
	
	The 2x2 table is structured as:
		| After Group | Before Group |
	----------------------------------
	Kmer Present | after       | before       |
	Kmer Absent  | kafter      | kbefore      |
	
	kafter = TOTAL_SAMPLES - after
	kbefore = TOTAL_SAMPLES - before
	"""
	kbefore = TOTAL_SAMPLES1 - before
	kafter = TOTAL_SAMPLES2 - after
	
				
	res = fisher_exact([[after, before], [kafter, kbefore]], alternative='greater')
	return np.array([res.statistic, res.pvalue])


def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(
        description='Perform Fisher\'s exact test on k-mer abundance data.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to the input file (e.g., /mnt/SAN96/RUNS/gurov_sa/for_vera_22_abund/22_itog.txt).\n'
             'Expected columns (space-separated, no header by default):\n'
             '1. kmer (int64)\n'
             '2. Before (uint64)\n'
             '3. After (uint64)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='output.csv',
        help='Path to the output CSV file (e.g., TEST_OUTPUT/test4.csv).\n'
             'Output columns: kmer, statistic, p_value. Appended if file exists.'
    )
    
    parser.add_argument(
        '-N1', '--samples1',
        type=int,
        required=True,
        help='Total number of samples (N) for the Fisher\'s exact test (e.g., 367).'
    )
    parser.add_argument(
        '-N2', '--samples2',
        type=int,
        required=False,
        default=False,
        help='Total number of samples (N) for the Fisher\'s exact test (e.g., 367).'
    )

    args = parser.parse_args()

    INPUT_FILE = args.input
    OUTPUT_FILE = args.output
    TOTAL_SAMPLES1 = args.samples1
    TOTAL_SAMPLES2 = args.samples2 if args.samples2 else args.samples1
    
    # Clear output file if it exists, as the loop is set to append ('a')
    # If you want to append across multiple runs, remove this part.
    try:
        with open(OUTPUT_FILE, 'w') as f:
            f.write('')
    except IOError as e:
        print(f"Warning: Could not clear or open output file {OUTPUT_FILE} for writing: {e}")
        return

    
  
    CHUNKSIZE = 200000 
    
    print(f"Processing input file: {INPUT_FILE}")
    print(f"Total samples (N): {TOTAL_SAMPLES1+TOTAL_SAMPLES2}")
    print(f"Output file: {OUTPUT_FILE}\n")

    pbar = tqdm.tqdm(desc="Processing Chunks", unit="chunk")
    
    try:
        for chunk in pd.read_table(
            INPUT_FILE,
            sep=' ', 
            chunksize=CHUNKSIZE, 
            header=None, 
            names=['kmer', 'Before', 'After'],
            dtype={'kmer':'int64', 'Before': 'uint64', 'After':'uint64'}
        ):
            # Apply the Fisher's exact test function to each row
            results = np.vstack(chunk.apply(
                lambda row: perform_fisher(row.Before, row.After, TOTAL_SAMPLES1, TOTAL_SAMPLES2), 
                axis=1
            ))
            
            # Assign results to new columns
            chunk['statistic'] = results[:, 0]
            chunk['p_value'] = results[:, 1]
            
            # Write only the required columns to the output file, appending results
            chunk[['kmer', 'statistic', 'p_value']].to_csv(
                OUTPUT_FILE, 
                mode='a', 
                header=False, 
                index=False
            )
            
            pbar.update(1)
            
            # Original script had 'break' here - removing it to process all chunks
            # break 
            
    except FileNotFoundError:
        print(f"\nError: Input file not found at {INPUT_FILE}")
    except Exception as e:
        print(f"\nAn error occurred during processing: {e}")
        
    pbar.close()
    print("\nProcessing complete.")


if __name__ == "__main__":
    main()
