# axiome-extensions
BETA personal scripts to extend or improve [AXIOME2](https://github.com/neufeld/AXIOME2) or [MESaS](https://github.com/neufeld/MESaS) functionality

Copyright Jackson M. Tsuji, 2018

Brief description of each script (more documentation coming!)
- `merge-illumina-reads.sh`: tool for upstream merging of read files (before running AXIOME) in case you ran the same samples over multiple Illumina runs and need to combine them. See help documentation within the script by running `./merge-illumina-reads.sh` in the terminal.
- `axiome_length_dist.sh`	Shows the length distribution of the reads merged by PANDAseq, to look for potential errors in the read merging process.
- `mesas-pcoa_JMT.R`	Improved PCoA plotting script, for use in RStudio. Cleaned up code and improved visuals for plotting.
- `mesas-pcoa_JMT-cl.R`	Improved PCoA plotting script (same as above), for use in command line. Effectively a replacement for the [original script](https://github.com/neufeld/MESaS/tree/master/scripts).
- `otutable_subset.R`		Makes a subset of the OTU table (e.g., to top 1%)
- `otutable_taxaplot.R`		Makes a taxa plot of the OTU table (like the [MetAnnotate barplot script](https://github.com/jmtsuji/metannotate-analysis))
