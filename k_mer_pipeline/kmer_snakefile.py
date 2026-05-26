# --- Configuration Section ---
"""
snakemake -s --use-conda -j 15 --cores 15 --config input_dir_1="" input_dir_2="" output_work_dir="" final_output_dir="" K_SIZE= references='' suffix1='_R1.fastq' suffix2='_R2.fastq'
"""

suffix_1 = config['suffix1']
suffix_2 = config['suffix2']

INPUT_PATHS = {
    "dataset1": config["input_dir_1"],
    "dataset2": config["input_dir_2"]
}


DIR_IDS = ["dataset1", "dataset2"]
K_SIZE = config['k_size']
REFERENCE_PATH = config['references']
# --- Sample and Target Definition ---

# Extract sample names and their source directory identifier (DIR_ID)
SAMPLES = []
for dir_id, input_path in INPUT_PATHS.items():
    samples_in_dir, = glob_wildcards(input_path + "{sample}" + suffix_1)
    for sample in samples_in_dir:
        SAMPLES.append((sample, dir_id))

SAMPLES_names = [s[0] for s in SAMPLES]
SAMPLES_dirs = [s[1] for s in SAMPLES]

REFERENCES, = glob_wildcards(REFERENCE_PATH + '{reference}' + '.fa')
print(REFERENCES)

# Print the extracted sample names and their associated directories for debugging
print("Samples (name, directory ID):", SAMPLES)

# 1. Targets for the per-dataset merged files (Intermediate step: kmer_to_dict.py)
# This uses the specific DIR_IDS list defined earlier
want_merged_datasets = expand(
    config["output_work_dir"] + '/{dir_id}/merged_ds.txt',
    dir_id=DIR_IDS
)

# 2. Target for the final merged output (Final step: merge_dict.py)
want_final_merge = expand(config["final_output_dir"] + '{reference}' + '_count.json', reference=REFERENCES)
want_final_merge2 = expand(config["final_output_dir"] + 'contigs_{reference}' + '.json', reference=REFERENCES)


rule all:
    input:
        want_final_merge, want_final_merge2# Set the final output as the ultimate target

# --- Workflow Rules ---

rule dsk_paired:
    input:
        read1=lambda wildcards: INPUT_PATHS[wildcards.dir_id] + wildcards.sample + suffix_1,
        read2=lambda wildcards: INPUT_PATHS[wildcards.dir_id] + wildcards.sample + suffix_2
    output:
        config["output_work_dir"] + '/{dir_id}/dsk_output/' + "{sample}_paired.h5"
    threads: 5
    params:
        k_size = K_SIZE
    conda:
        "envs/dsk.yaml"

    shell:
        """
        dsk -file {input.read1},{input.read2} -out {output} -kmer-size {params.k_size}
        """

rule dsk2ascii_paired:
    input:
        config["output_work_dir"] + '/{dir_id}/dsk_output/' + "{sample}_paired.h5"
    output:
        # Note: This is now the input for the new 'merge_dataset_kmers' rule
        config["output_work_dir"] + '/{dir_id}/dsk_ascii_output/' + "{sample}.dsk_paired.txt"
    conda:
        "envs/dsk.yaml"

    shell:
        "dsk2ascii -file {input} -out {output}"

# --- New Rules for Merging ---

rule merge_dataset_kmers:
    # This rule takes ALL the dsk2ascii outputs for a given {dir_id}
    input:
        # Filter the combined SAMPLES list to only include files for the current dir_id
        # Then, construct the input paths using 'expand'
        lambda wildcards: expand(
            config["output_work_dir"] + '/{dir_id}/dsk_ascii_output/' + "{sample}.dsk_paired.txt",
            dir_id=wildcards.dir_id,
            sample=[s[0] for s in SAMPLES if s[1] == wildcards.dir_id]
        )
    output:
        config["output_work_dir"] + '/{dir_id}/merged_ds.txt'
    params:
        # Create a string of all input files separated by spaces for the script
        files_list=lambda wildcards, input: " ".join(input)
    # Define a generic environment or ensure 'kmer_to_dict.py' is available
    # conda:
    #     "envs/python_scripts.yaml"

    shell:
        """
        python scripts/kmer_to_dict.py \
            -i {params.files_list} \
            -o {output}
        """

rule final_merge_comparison:
    # This rule takes the two merged files from the previous step
    input:
        before=config["output_work_dir"] + '/dataset1/merged_ds.txt',
        after=config["output_work_dir"] + '/dataset2/merged_ds.txt'
    output:
        config["output_work_dir"] + 'final_merged_comparison.txt'
    # conda:
    #     "envs/python_scripts.yaml"

    shell:
        """
        python scripts/merge_dict.py  \
            -b {input.before} \
            -a {input.after} \
            -o {output}
        """

rule run_fisher:
    input:
        rules.final_merge_comparison.output
    output:
        config["output_work_dir"] + 'fisher_exact_raw.txt'
    params:
        n_samples_before = sum(1 for s in SAMPLES if s[1] == "dataset1"),
        n_samples_after = sum(1 for s in SAMPLES if s[1] == "dataset2")
    shell:
        """
        python scripts/k_mer_fisher.py --input {input} --output {output} --samples1 {params.n_samples_before} --samples2 {params.n_samples_after}
        """

rule filter_fisher:
    input:
        rules.run_fisher.output
    output:
        config["output_work_dir"] + 'fisher_exact_filter.txt'
    shell:
        """
        n_k_mers=$(wc -l {input} | awk '{{print $1}}')
		p_value=0.05
		threshold=$(echo "scale=17; $p_value / $n_k_mers" | bc -l)
		awk -F ',' -v T="$threshold" '$3 < T' {input} > {output}
        """



rule reference_to_k_mer:
    input:
         REFERENCE_PATH + '{reference}' + '.fa'
    output:
        config["output_work_dir"] + '{reference}' + '.pkl'
    params:
        k_size = K_SIZE
    shell:
        """
        python scripts/contigs_to_kmers.py --input {input} --output {output}  -k {params.k_size}
        """

rule measure_statistics:
    input:
        raw_k_mers = rules.run_fisher.output,
        filtered_k_mers = rules.filter_fisher.output,
        reference = rules.reference_to_k_mer.output
    output:
        kmers = config["output_work_dir"] + '{reference}' + '_statistics.json',
        statistics = config["final_output_dir"] + '{reference}' + '_count.json'
    params:
        k_size = K_SIZE
    shell:
        """
        python scripts/measure_22_mers_statistics.py --all_kmers {input.raw_k_mers} --positive_kmers {input.filtered_k_mers} --wuhan_pkl {input.reference} --output_json {output.kmers} --output_json_count {output.statistics} --k {params.k_size}
        """

rule k_mer_int_to_seq:
    input:
        rules.filter_fisher.output
    output:
        config["output_work_dir"] + 'k_mers_sequences.fasta'
    shell:
        """
        python scripts/int_to_kmer.py --input {input} --output {output} --k_size 22 
        """

rule megahit_k_mers:
    input:
        rules.k_mer_int_to_seq.output
    output:
        config["output_work_dir"] + 'MEGAHIT/final.contigs.fa'
    conda:
        'envs/megahit.yaml'
    threads:
        5
    shell:
        f"""
        rm -r {config['output_work_dir'] + 'MEGAHIT'}
        megahit -r {{input}} -o {config['output_work_dir']}/MEGAHIT  --k-list 15,17,19,21 -t 5
        """

rule blastn:
    input:
        raw_k_mers = rules.megahit_k_mers.output,
        reference = REFERENCE_PATH + '{reference}' + '.fa'
    output:
        config["output_work_dir"] + 'BLASTN_RESULTS/{reference}' + '.tsv'
    conda:
        'envs/blast.yaml'
    shell:
        """
        blastn -query {input.raw_k_mers} -subject {input.reference} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > {output}
        """

rule statistics_blastn:
    input:
        k_mer_contigs = rules.megahit_k_mers.output,
        reference = REFERENCE_PATH + '{reference}' + '.fa',
        blastn_results = rules.blastn.output
    output:
        config["final_output_dir"] + 'contigs_{reference}' + '.json'
    shell:
        """
        python scripts/get_contigs_results.py {input.k_mer_contigs} {input.reference} {input.blastn_results} {output}
        """

