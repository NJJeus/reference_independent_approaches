import pandas as pd
import numpy as np
import tqdm
import pickle as pkl
import json # Import json for JSON output
import argparse
import sys # Import sys for printing

# --- Helper function for JSON serialization ---

# Function to convert NumPy arrays to lists for JSON serialization
def convert_to_json_friendly(obj):
    """Recursively converts NumPy arrays to Python lists."""
    if isinstance(obj, np.ndarray):
        # Convert NumPy array elements to list, ensuring they are Python native types
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_json_friendly(v) for k, v in obj.items()}
    return obj

# --- Core Logic Functions (Kept as is) ---

dreverse = str.maketrans('0123', '1032')
def get_reverse_compl(n, k):
    """Calculates the reverse complement of a k-mer (represented as a base-4 integer)."""
    # Convert integer to base-4, pad with leading zeros, reverse, and translate (0<->1, 2<->3)
    quart = np.base_repr(n, 4).rjust(k, '0')[::-1].translate(dreverse)
    return int(quart, base=4)

# Vectorize the function for array operations
get_reverse_compl = np.vectorize(get_reverse_compl)

def get_full_seq(arr, k):
    """Returns the k-mer set including both forward and reverse complement strands."""
    return np.concatenate([arr, get_reverse_compl(arr, k=k)])

# --- Main Script Logic ---

def main(all_kmers_path, positive_kmers_path, wuhan_pkl_path, output_json_path, output_json_path_count, k): # Changed pkl to json
    """
    Main function to load k-mer lists, perform classification against reference
    sequences (all keys in the pickle file), and save the results.
    """

    # 1. Load Common K-mer Sets

    print('Read all k-mers file')
    try:
        data_all_k_mers = pd.read_csv(all_kmers_path, usecols=[0], header=None).to_numpy().ravel()
    except FileNotFoundError:
        print(f"Error: All k-mers file not found at {all_kmers_path}", file=sys.stderr)
        return

    print('Read positive k-mers')
    try:
        k_mers_positive = pd.read_csv(positive_kmers_path, usecols=[0], header=None).to_numpy().ravel()
    except FileNotFoundError:
        print(f"Error: Positive k-mers file not found at {positive_kmers_path}", file=sys.stderr)
        return

    print('Find negative k-mers (All K-mers - Positive K-mers)')
    k_mers_negative = np.setdiff1d(data_all_k_mers, k_mers_positive)

    print('Read reference k-mers dictionary')
    try:
        # Load the dictionary of reference k-mer sets
        wuhan_data_dict = pkl.load(open(wuhan_pkl_path, 'rb'))
    except FileNotFoundError:
        print(f"Error: Reference k-mers pickle file not found at {wuhan_pkl_path}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error loading pickle file: {e}", file=sys.stderr)
        return

    # 2. Iterate and Analyze for All Keys

    all_results = {}
    all_results_count = {}
    total_keys = len(wuhan_data_dict)
    print(f'Starting analysis for {total_keys} sequence keys found in the pickle file...')

    # Iterate over all keys in the loaded dictionary
    for key in wuhan_data_dict:
        print(f'\nProcessing key: {key}')

        # a. Prepare the current reference k-mer set (including reverse complements)
        wuhan_k_mers_raw = np.unique(wuhan_data_dict[key].ravel())
        wuhan_k_mers = get_full_seq(wuhan_k_mers_raw, k=k)

        # b. Calculate statistics (TP, FP, TN, FN)
        pbar = tqdm.tqdm(total=4, desc=f"Calculating stats for {key}")

        # True Positive (Positive in Test AND in Reference)
        true_positive = np.intersect1d(k_mers_positive, wuhan_k_mers); pbar.update(1)

        # False Positive (Positive in Test BUT NOT in Reference)
        false_positive = np.setdiff1d(k_mers_positive, wuhan_k_mers); pbar.update(1)

        # True Negative (Negative in Test AND NOT in Reference) - Note: This is within the scope of ALL k-mers
        true_negative = np.setdiff1d(k_mers_negative, wuhan_k_mers); pbar.update(1)

        # False Negative (Negative in Test BUT IN Reference)
        false_negative= np.intersect1d(k_mers_negative, wuhan_k_mers); pbar.update(1)

        pbar.close()

        # c. Find k-mers in reference that were not in the initial ALL k-mers set
        lost_wuhan = np.setdiff1d(wuhan_k_mers, data_all_k_mers)

        # d. Store results for the current key
        all_results[key] = {
            'TP': true_positive,
            'TN': true_negative,
            'FP': false_positive,
            'FN': false_negative,
            'Total_K_MERS': wuhan_k_mers
        }
        all_results_count[key] = {k:len(v) for k, v in all_results[key].items()}

    # 3. Save Final Results

    print('\nAnalysis complete. Dumping results.')
    
    # Convert NumPy arrays to lists for JSON serialization
    json_friendly_results = convert_to_json_friendly(all_results)
    json_friendly_results_count = convert_to_json_friendly(all_results_count)

    
    try:
        # SAVE AS JSON
        with open(output_json_path, 'w') as f:
            json.dump(json_friendly_results, f, indent=4)
        with open(output_json_path_count, 'w') as f:
            json.dump(json_friendly_results_count, f, indent=4)

        print(f'Results successfully dumped to {output_json_path}')
    except Exception as e:
        print(f"Error saving output JSON file: {e}", file=sys.stderr)


if __name__ == '__main__':
    # 4. Argparse Setup
    parser = argparse.ArgumentParser(description='Analyze k-mer statistics by comparing positive/negative sets against multiple reference sequences.')

    parser.add_argument('--all_kmers',
                        required=True,
                        help='Path to CSV file containing all analyzed k-mers (e.g., ../KMER_OUTPUT/k22_fisher_p_values.csv)')

    parser.add_argument('--positive_kmers',
                        required=True,
                        help='Path to CSV file containing positive/significant k-mers (e.g., ../KMER_OUTPUT/k22_fisher_p_values_bonferoni.csv)')

    parser.add_argument('--wuhan_pkl',
                        required=True,
                        help='Path to the pickle file containing reference k-mers (a dict where keys are sequence names and values are k-mer arrays, e.g., ../KMER_OUTPUT/wuhan_hu_1_k_mers.pkl)')

    # Changed argument name and help text
    parser.add_argument('--output_json',
                        required=True,
                        help='Path where the output statistics JSON file will be saved (e.g., ../RESULTS/k22_mers_statistics_fisher_bonferoni.json)')
    parser.add_argument('--output_json_count',
                        required=True,
                        help='Path where the output statistics JSON file will be saved (e.g., ../RESULTS/k22_mers_statistics_fisher_bonferoni.json)')
    parser.add_argument('--k',
                        required=True,
                        type=int,
                        help='k size')


    args = parser.parse_args()

    # Changed pkl argument name to json
    main(args.all_kmers, args.positive_kmers, args.wuhan_pkl, args.output_json, args.output_json_count, args.k)
