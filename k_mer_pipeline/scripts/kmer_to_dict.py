from tqdm import tqdm
from pathlib import Path
import sys
import math
import argparse
import pickle

def merge_kmer_files(file_list, output_file):

    translation_table = str.maketrans('ATGC', '0123')
    kmer_counts = {}
    processed_files = 0

    # Обработка каждого файла
    for file_path in tqdm(file_list, desc="Processing files"):
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                kmer = int(parts[0].translate(translation_table), base=4)

                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1


        processed_files += 1


    for key in tqdm(list(kmer_counts.keys())):
        if kmer_counts[key] <= 5:  
            del kmer_counts[key]

  
    pickle.dump(kmer_counts, open(output_file.replace('.txt', '.pkl'), 'wb'))

    k_mers = list(kmer_counts.keys())
    k_mers.sort()

    # Сортировка и запись результатов
    with open(output_file, 'w') as out_f:
        # Сортируем k-mer по алфавиту
        for kmer in tqdm(k_mers, desc="Writing output"):
            out_f.write(f"{kmer} {kmer_counts[kmer]}\n")

    print(f"\nProcessed {processed_files} files")
    print(f"Total unique k-mers: {len(kmer_counts)}")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(
        description='Merge k-mer count files, converting DNA sequences to base-4 integers.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        type=str,
        nargs='+',  # <--- KEY CHANGE: Expects one or more file paths
        required=True,
        help='A list of one or more k-mer count file paths (e.g., file1.txt file2.txt ...).'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='merged_kmer_dict.txt',
        help='Name of the output file (e.g., after_kmer_dict.txt). Default: %(default)s'
    )

    args = parser.parse_args()
    file_list= args.input
    output_file = args.output

    # Сортируем файлы для предсказуемости
    file_list.sort()

    # Запускаем обработку
    merge_kmer_files(file_list, args.output)

