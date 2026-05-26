import argparse
from glob import glob
from Bio import SeqIO
import numpy as np
import pickle
import os
import tqdm

def lugat(kmer):
    translation_table = str.maketrans('ATGC', '0123')
    return int(kmer.translate(translation_table), base=4)

def process_fasta_files(input, output, k_size):
    kmer_dict = {}
    
    # Get all fasta files in input directory
    fasta_file = input
    

    file_name = fasta_file.split('/')[-1].split('.')[0]
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            sequence = str(record.seq).upper()
            kmers = []
            
            # Generate 22-mers with step 1
            for i in range(len(sequence) - (k_size-1)):
                kmer = sequence[i:i+k_size]
                if all(c in 'ATGC' for c in kmer):  # Only process valid kmers
                    kmers.append(lugat(kmer))
            
            # Convert to numpy array and store in dict
            kmer_dict[record.id] = np.array(kmers, dtype=np.int64)
            # Save the dictionary to pkl file
            with open(output, 'wb') as f:
                pickle.dump(kmer_dict, f)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process FASTA files to extract and convert 22-mers.')
    parser.add_argument('--input', help='Input directory containing FASTA files')
    parser.add_argument('--output', help='Output pickle file to save the kmer dictionary')
    parser.add_argument('-k', '--k_size', type=int)

    args = parser.parse_args()
    
    if os.path.isfile(args.input):
        input_var = [args.input]
        output_var = [args.output]
    elif os.path.isdir(args.input):
        input_var = glob(f"{args.input}/*.fa") + glob(f"{args.input}/*.fasta")
        print(input_var)
        if not os.path.isdir(args.output):
            os.mkdir(args.output)
        output_var = [f"{args.output}/{i.split('/')[-1].split('.')[0]}.pkl" for i in input_var]    
    
    for input_file, output_file in zip(input_var, output_var):
        process_fasta_files(input_file, output_file, args.k_size)
