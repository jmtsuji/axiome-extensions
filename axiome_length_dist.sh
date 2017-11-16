#!/bin/bash
# Created Feb 4, 2017, by Jackson Tsuji (Neufeld lab PhD student)
# Description: Generates length distribution plots for all PANDAseq output for a given AXIOME run
# **Run script from within the main AXIOME folder of interest.
# Output goes to misc_analysis/length_distributions
# Last updated: Feb 3, 2017

# Basic script stuff (from Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 397):
set -e
set -u
set -o pipefail

script_version=1.0.0
date_code=$(date '+%y%m%d')

# Settings:
work_dir=$(pwd) # Could be modified in the future if desired
combined_file_basename="all_samples"

echo "Running $(basename $0) version $script_version on ${date_code} (yymmdd)."
echo ""

cd $work_dir

# Test that PANDAseq directory exists, and exit if it does not
if [ ! -d source ]
# From http://stackoverflow.com/a/4906665, accessed Feb. 4, 2017
then
    print "Did not find 'source' directory for PANDAseq files. Job terminating."
    exit 1
fi

# Make output directory tree
mkdir -p misc_analysis/length_distributions/raw
mkdir -p misc_analysis/length_distributions/tables
mkdir -p misc_analysis/length_distributions/plots
mkdir -p misc_analysis/length_distributions/combined_samples

# Remove combined samples file if already exists (as regular or gzipped file)
rm -f misc_analysis/length_distributions/combined_samples/${combined_file_basename}_seqlengths.txt ${combined_file_basename}_seqlengths.txt.gz
# see http://stackoverflow.com/a/31318262, accessed Feb. 4, 2017

# Get list of PANDAseq output files (store as an array for looping)
# Note: name will be stored with full file path and extension for simplicity later on
pandaseq_files=($(find source -name "*.fasta.gz" -type f))
# Got help from http://stackoverflow.com/questions/2961673/find-missing-argument-to-exec, post by Marian on June 2, 2010; accessed May 13, 2016
echo "Found ${#pandaseq_files[@]} PANDAseq output files to generate length distributions from."
echo ""

# Making bakups
echo "Making length distribution plots using$(seq_length_plotter.R 2>&1 | head -n 1 | cut -d ',' -f 2) for the following samples..."
echo ""

for pandaseq_file in ${pandaseq_files[@]}
do
    # Make base name of file
    pandaseq_file_base=$(basename $pandaseq_file .fasta.gz)
    echo "${pandaseq_file_base}.fasta.gz"

    # Generate length distribution using awk
    zcat ${pandaseq_file} | grep -v '^>' | awk '{ print length($0) }' > misc_analysis/length_distributions/raw/${pandaseq_file_base}_seqlengths.txt

    # Add seq lengths to combined samples file
    cat misc_analysis/length_distributions/raw/${pandaseq_file_base}_seqlengths.txt >> misc_analysis/length_distributions/combined_samples/${combined_file_basename}_seqlengths.txt

    # Generate length table and plot using R script
    seq_length_plotter.R misc_analysis/length_distributions/raw/${pandaseq_file_base}_seqlengths.txt ${pandaseq_file_base}

    # Move table and PDF to correct folders -- something to be improved in the R script later, perhaps
    mv ${pandaseq_file_base}_lengths_table.csv misc_analysis/length_distributions/tables
    mv ${pandaseq_file_base}_lengths_plot.pdf misc_analysis/length_distributions/plots

    # Gzip sequence lengths file
    gzip --force misc_analysis/length_distributions/raw/${pandaseq_file_base}_seqlengths.txt
done

echo ""
echo "Done."

echo ""

echo "Generating length distribution for all samples in AXIOME run combined..."

# Generate length table and plot using R script
seq_length_plotter.R misc_analysis/length_distributions/combined_samples/${combined_file_basename}_seqlengths.txt ${combined_file_basename}

# Move table and PDF to correct folders -- something to be improved in the R script later, perhaps
mv ${combined_file_basename}_lengths_table.csv misc_analysis/length_distributions/combined_samples
mv ${combined_file_basename}_lengths_plot.pdf misc_analysis/length_distributions/combined_samples

# Gzip sequence lengths file
gzip --force misc_analysis/length_distributions/combined_samples/${combined_file_basename}_seqlengths.txt

echo "Done."
echo ""

echo "$(basename $0): finished. Output can be fount in misc_analysis/length_distributions"
echo ""
