#!/bin/bash
# merge-illumina-runs.sh
# Copyright Jackson M. Tsuji, 2017
# Neufeld lab, University of Waterloo, Canada
# Created Feb. 5th
# Description: Bulk merges reads from the same sample but sequenced over multiple Illumina runs into the same file.

# Basic script stuff (from Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 397):
set -e
set -u
set -o pipefail

script_version=1.0.1

# If no input is provided, exit out and provide help
if [ $# == 0 ]
    then
    printf "\n$(basename $0): Merges reads from the same sample sequenced over multiple Illumina runs into the same file.\n"
    printf "Version: ${script_version}\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Usage: $(basename $0) merge_description_file.tsv number_of_runs output_folder_path | tee $(basename $0 .sh).log \n\n"
    printf "Usage details:\n"
    printf "* input files: MUST BE gzipped FastQ files with extension .fastq.gz\n"
    printf "* merge_description_file.tsv: tab-separated file with headers. See more detailed description below.\n"
	printf "* number_of_runs: number of sequencing runs to merge (e.g., 2).\n"
	printf "* output_folder_path: path to where you want the merged reads to be output.\n"
	printf "\n\n"
	printf "Column names for merge_description_file.tsv (MUST EXACTLY MATCH):\n"
	printf "* sampleID: the same of the sample. The merged file will be this name followed by _R1.[extension] or _R2.[extension]\n"
	printf "    **Importantly, files to be merged MUST BE gzipped FastQ files with extension .fastq.gz\n"
	printf "* filepath_run1_R1: filepath to the R1 reads of the first Illumina run for that sample\n"
	printf "* filepath_run1_R2: filepath to the R2 reads of the first Illumina run for that sample\n"
	printf "* filepath_run2_R1: filepath to the R1 reads of the second Illumina run for that sample\n"
	printf "* filepath_run2_R1: filepath to the R2 reads of the second Illumina run for that sample\n"
	printf "* ...repeat for as many runs as you have.\n\n"
	printf "Caveats to using this script:\n"
	printf "* Note that this script only works for non-interleaved paired-end reads.\n"
	printf "* Also note that you'll likely run into problems if you used different barcodes in the two Illumina runs for the same sample. If different barcodes were used, I'd recommend removing any primers/barcode tags upstream of merging the read files together to help prevent issues with tools such as PANDAseq.\n"
    printf "\n\n"
    exit 1
fi

# Set variables from user input:
merge_description_file=$1
number_of_runs=$2
output_folder_path=$3

function check_headers {
	# Description: checks that the table headers are okay and exits out if there is an error.
	# GLOBAL params: merge_description_file; number_of_runs
	# Input params: none
	# Output: none
	
	# Build an array of headers
	headers=($(head -n 1 ${merge_description_file}))
	
	# Count total number of columns
	total_cols=${#headers[@]}
	
	# Check sampleID
	header_name_check 0 "sampleID"
	
	# Check total columns matches expected number given number of runs
	expected_cols=$((${number_of_runs}*2+1))
	if [ $total_cols != $expected_cols ]; then
		echo "ERROR: expected to find ${expected_cols} columns in $(basename ${merge_description_file}), but found ${total_cols} instead. Exiting..."
		exit 1
	fi
	
	# Check run columns match expected names
	for i in $(seq 1 $((${total_cols}-1))); do
		# Intentionally skips column 0, because that was already checked above.
		
		expected_run_num=$(((${i}+1)/2))
		# For R1 vs. R2 reads
		expected_fwd_rev=$(((${i}+1)%2+1))
		
		expected_col_name="filepath_run${expected_run_num}_R${expected_fwd_rev}"
		
		header_name_check ${i} ${expected_col_name}
		
	done
	
	echo "Test passed: merge_description_file looks okay."
	echo ""
	
}

function header_name_check {
	# Description: checks if a column's header matches the expected name and exits if it does not match
	# GLOBAL params: merge_description_file; number_of_runs; headers
	# Input params: col_num (column number of header, starting numbering at 0); expected_name (exact name expected for header)
	# Output: none
	
	local col_num=$1
	local expected_name=$2
	
	# Get header name for column
	header_name=${headers[col_num]}
	
	if [ $header_name != $expected_name ]; then
		echo "ERROR: column $((${col_num}+1)) header '${header_name}' does not match expected name '${expected_name}'. Exiting..."
		exit 1
	fi

}

function merge_all_reads {
	# Description: iteratively merges reads for each sample as R1 and R2 files
	# GLOBAL params: merge_description_file; number_of_runs; output_folder_path
	# Input params: none
	# Output: writes files to disk
	
	echo "Merging reads from samples..."
	
	merge_R1_reads
	merge_R2_reads
	
	echo ""
	
}

function merge_R1_reads {
	# Description: iteratively merges the R1 reads from a given sample
	# GLOBAL params: merge_description_file; number_of_runs; output_folder_path
	# Input params: none
	# Output: writes files to disk
	
	# Temporarily change the internal fields separator (IFS) so that whitespaces do not create new entries. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
	OFS="$IFS"
	IFS=$'\t'$'\n'
	
	# Load basic info from merge_description_file
	sample_names=($(tail -n +2 ${merge_description_file} | cut -d $'\t' -f 1))
	number_of_samples=${#sample_names[@]}
	
	echo "R1 reads:"
	echo ""
	
	for i in $(seq 1 ${number_of_samples}); do
		# Reorder counter to start at zero
		j=$((i-1))
		
		# Determine name and path of output file
		sample_name=${sample_names[${j}]}
		output_path="${output_folder_path}/${sample_name}_R1.fastq.gz"
		
		# Make sure the output file doens't already exist
		if [ -f ${output_path} ]; then
			echo "ERROR: file ${output_path} already exists. Please delete before running this script. Exiting..."
			exit 1
		fi
		
		# Iteratively concatenate files from all runs together
		for run in $(seq 1 ${number_of_runs}); do
			# Get the sample names for that run
			run_col=$((${run}*2))
			run_filepaths=($(tail -n +2 ${merge_description_file} | cut -d $'\t' -f ${run_col}))
			input_filepath=${run_filepaths[${j}]}
			
			# Add that file's contents onto the output file
			echo "Run ${run}: ${input_filepath} --> ${output_path}"
			cat ${input_filepath} >> ${output_path}
		done
		
	done

	echo ""
	echo ""

	# Fix the IFS
	IFS="$OFS"

}

function merge_R2_reads {
	# Description: iteratively merges the R2 reads from a given sample
	# GLOBAL params: merge_description_file; number_of_runs; output_folder_path
	# Input params: none
	# Output: writes files to disk
	
	# Temporarily change the internal fields separator (IFS) so that whitespaces do not create new entries. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
	OFS="$IFS"
	IFS=$'\t'$'\n'
	
	# Load basic info from merge_description_file
	sample_names=($(tail -n +2 ${merge_description_file} | cut -d $'\t' -f 1))
	number_of_samples=${#sample_names[@]}
	
	echo "R2 reads:"
	echo ""
	
	for i in $(seq 1 ${number_of_samples}); do
		# Reorder counter to start at zero
		j=$((i-1))
		
		# Determine name and path of output file
		sample_name=${sample_names[${j}]}
		output_path="${output_folder_path}/${sample_name}_R2.fastq.gz"
		
		# Make sure the output file doens't already exist
		if [ -f ${output_path} ]; then
			echo "ERROR: file ${output_path} already exists. Please delete before running this script. Exiting..."
			exit 1
		fi
		
		# Iteratively concatenate files from all runs together
		for run in $(seq 1 ${number_of_runs}); do
			# Get the sample names for that run
			run_col=$(((${run}*2)+1))
			run_filepaths=($(tail -n +2 ${merge_description_file} | cut -d $'\t' -f ${run_col}))
			input_filepath=${run_filepaths[${j}]}
			
			# Add that file's contents onto the output file
			echo "Run ${run}: ${input_filepath} --> ${output_path}"
			cat ${input_filepath} >> ${output_path}
		done
		
	done
	
	echo ""
	echo ""

	# Fix the IFS
	IFS="$OFS"

}

function main {
	
	echo "Running $(basename $0), version $script_version at $(date)."
	echo ""
	echo "Merging files from ${number_of_runs} Illumina runs. Merging $(tail -n +2 ${merge_description_file} | wc -l) samples."
	echo ""

	# Test provided input file
	check_headers

	# Merge reads
	merge_all_reads
	
	echo ""
	echo "Read merging finished. Output saved in folder '$(realpath ${output_folder_path})'."
	echo ""

	echo "$(basename $0): finished at $(date)."
	echo ""

}

main
