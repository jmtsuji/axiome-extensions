# Subset OTU table
# Adapted by Jackson M. Tsuji from BubblePlot script by Sara Victoria Coyotzi Alcaraz (Neufeld Lab)
# Description: Subsets OTUs from an OTU table with relatve abundance (in at least one sample) of specified value (e.g. >= 1%). This can be used as input for FAPROTAX
#               Works for rarefied and non-rarefied OTU tables.
# Version 1.1.1
# Last updated: May 5, 2017 (by Jackson Tsuji)

#####################################################
## Script settings: #################################
## Basic settings:
# Working directory
setwd("/Users/JTsuji/Downloads/")
# Functional table
otutable_name <- "otu_table.tab"
# Minimum percent abundance to subset
min_abundance <- 5
# Optional: print list of top OTUs per sample?
print_top_otu_list <- TRUE
# Optional: print list of total read counts for each sample before subset?
print_total_read_counts <- TRUE
# Optional: print taxa plot of top organisms in OTU table?
print_taxa_plot <- TRUE
#####################################################
#####################################################

# Import OTU table
OTU_table <- read.table(otutable_name , header = TRUE, sep="\t", comment.char = "", skip = 1)

# Standardize OTU ID header
colnames(OTU_table)[1] <- "OTU_ID"

# Check if representative sequence was included, and exclude if present, but keep a backup
if (is.na(match("ReprSequence", colnames(OTU_table))) == FALSE) {
  ReprSequence_exists <- TRUE
  print("Removing ReprSequence column")
  ReprSequence_col_num <- match("ReprSequence", colnames(OTU_table))
  OTU_table_ReprSequence <- OTU_table[,c(1,ReprSequence_col_num)]
  OTU_table <- OTU_table[,-(ReprSequence_col_num)]
} else {
  ReprSequence_exists <- FALSE
}

# Same for Consensus Lineage, but keep a backup
if (is.na(match("Consensus.Lineage", colnames(OTU_table))) == FALSE) {
  Consensus.Lineage_exists <- TRUE
  print("Removing Consensus.Lineage column")
  Consensus.Lineage_col_num <- match("Consensus.Lineage", colnames(OTU_table))
  OTU_table_lineage <- OTU_table[,c(1,Consensus.Lineage_col_num)]
  OTU_table <- OTU_table[,-(Consensus.Lineage_col_num)]
} else {
  Consensus.Lineage_exists <- FALSE
}

# Give the table informative row names
rownames(OTU_table)<- OTU_table$OTU_ID

# Remove redundant OTU_ID column
OTU_table <- OTU_table[,-1]

## Make secondary table with relative abundances
OTU_proportion <- OTU_table
# Divide each column by the column sum.
col_sums <- as.numeric(colSums(OTU_table))
for (i in c(1:length(col_sums))){
    OTU_proportion[,i] <- as.numeric(OTU_proportion[,i]/col_sums[i])
}
# Change proportions to percentages
OTU_proportion <- OTU_proportion * 100

# Subset those above or at 1% abundance for each sample and keep only those as a list of data frames (different length for each data frame)
OTU_select_prop <- lapply(1:ncol(OTU_proportion), function(x) { as.data.frame(subset(OTU_proportion, OTU_proportion[,x] >= min_abundance, select=x)) } )

### Remove sample if no OTUs are above threshold
# First, make empty vectors to fill into
removed_objects <- character()
removed_objects_num <- numeric()
# Then, run for loop to identify elements to remove
for (list_num in 1:ncol(OTU_proportion)) {
  # Check is list object is empty
  if (length(OTU_select_prop[[list_num]][,1]) == 0 ) {
    # Get name of element, for reporting
    removed_objects <- append(removed_objects, names(OTU_select_prop[[list_num]]))
    # Get number of list element for removal later
    removed_objects_num <- append(removed_objects_num, list_num)
  }
}
# Remove elements if anything was found
if (length(removed_objects) > 0 ){
  OTU_select_prop <- OTU_select_prop[-removed_objects_num]
}

# Report to user
if (length(removed_objects) > 0) {
  # Make name of output file
  output_removed_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_removed_samples.txt", sep = "")
  
  # Print warning (with two options for grammatical correctness)
  if (length(removed_objects) == 1) {
    print(paste("Removed ", length(removed_objects), " sample from OTU plot because it did not have any OTUs above the given threshold. Name of removed sample has been saved as ", output_removed_name, ". Note that the subset OTU table will still include the sample, however.", sep = ""))
  } else {
    print(paste("Removed ", length(removed_objects), " samples from OTU plot because they did not have any OTUs above the given threshold. Names of removed samples have been saved as ", output_removed_name, ". Note that the subset OTU table will still include the samples, however.", sep = ""))
  }
  
  # Write output
  write(removed_objects, file = output_removed_name, ncolumns = 1)
}
  

# Assign the sample names and OTU IDs as columns in the data frame
OTU_select_prop <- lapply(1:length(OTU_select_prop), function(x) { OTU_select_prop[[x]]$SampleName <- colnames(OTU_select_prop[[x]]); OTU_select_prop[[x]]$OTU_ID <- rownames(OTU_select_prop[[x]]); OTU_select_prop[[x]] } )

# Name the column for proportion as Percentage_abundance
OTU_select_prop <- lapply(1:length(OTU_select_prop), function(x) { colnames(OTU_select_prop[[x]])[1] = "Percentage_abundance"; OTU_select_prop[[x]] } )

# Give names to the list variables, which are of data frame type 
# If needed, remove from the colnames vector any samples that were removed from the data frame earlier
if (length(removed_objects) > 0) {
  select_names <- colnames(OTU_proportion)[-removed_objects_num]
  names(OTU_select_prop) <- select_names
} else {
  names(OTU_select_prop) <- colnames(OTU_proportion)
}

# From list of data frames to a single data frame
OTU_select_propDF <- do.call(rbind.data.frame, OTU_select_prop)

# Eliminate old row names
rownames(OTU_select_propDF) <- NULL

# Re-order columns for clarity
OTU_select_propDF <- OTU_select_propDF[,c(3,2,1)]

# Make OTU ID numeric for later
OTU_select_propDF$OTU_ID <- as.numeric(OTU_select_propDF$OTU_ID)

# Add consensus lineage and representative sequences back on
library(plyr)
library(dplyr)

if (Consensus.Lineage_exists == TRUE) {
  OTU_select_propDF <- dplyr::left_join(OTU_select_propDF, OTU_table_lineage, by = "OTU_ID")
}
if (ReprSequence_exists == TRUE) {
  OTU_select_propDF <- dplyr::left_join(OTU_select_propDF, OTU_table_ReprSequence, by = "OTU_ID")
}

# Print above table for reference for the user, if desired
if (print_top_otu_list == TRUE) {
  # Modify OTU table name as the output name for the table
  top_otu_list_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_top_otu_names.tsv", sep = "")
  
  # Write the table
  write.table(OTU_select_propDF, file = top_otu_list_name, sep = "\t", row.names = FALSE, col.names = TRUE)
}

## Make a taxa plot from subset, if desired
if (print_taxa_plot == TRUE) {
  # Load plotting packages
  library(ggplot2)   
  library(grid)
  
  # Check Consensus.Lineage exists
  if (Consensus.Lineage_exists == FALSE) {
    print("Cannot make a taxaplot because no Consensus.Lineage was provided in input OTU table")
  } else { 
    # Drop old leftover levels in data frame
    OTU_select_propDF$Consensus.Lineage <- droplevels(OTU_select_propDF$Consensus.Lineage)
    
    # Get length of legend
    legend_length <- length(levels(OTU_select_propDF$Consensus.Lineage))
    
    # Decide on number of legend columns as a ratio of legend length (30 to a column)
    legend_ncol <- ceiling(legend_length / 30)
    
    # Plot the Data (old version, then new version specific for paper)
    taxa_plot <- ggplot(OTU_select_propDF, aes(SampleName)) +
      geom_bar(aes(weight = Percentage_abundance, fill = Consensus.Lineage), position = "stack") +
      #scale_fill_manual() +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 9), 
            panel.border = element_rect(colour = "transparent"),
            axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
            legend.text = element_text(size = 5, face = "italic"), legend.title = element_blank(),
            legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", legend.position = "bottom", plot.margin = unit(c(2,2,2,2), "mm")) +
      guides(fill = guide_legend(ncol = legend_ncol)) +
      xlab("Sample") +
      ylab("Cumulative percentage abundance")
      
    print(taxa_plot)
    
    # Modify OTU table name as the output name for the table
    output_PDF_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_taxa_plot.pdf", sep = "")
    
    # Save the plot
    ggsave(output_PDF_name, width = 400, height = 300, units = c("mm"))
  }
}

## Make final subset OTU table
# Get top OTU IDs
top_OTUs <- unique(OTU_select_propDF$OTU_ID)

# Make a data frame
OTU_table_top <- data.frame("OTU_ID" = top_OTUs)

# Match OTU IDs to original OTU table
library(plyr)
library(dplyr)

# Add OTU_ID column back to original OTU table and delete row names
OTU_table$OTU_ID <- as.numeric(rownames(OTU_table))
rownames(OTU_table) <- NULL

OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table, by = "OTU_ID")

# Add consensus lineage and representative sequences back on
if (Consensus.Lineage_exists == TRUE) {
  OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_lineage, by = "OTU_ID")
}
if (ReprSequence_exists == TRUE) {
  OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_ReprSequence, by = "OTU_ID")
}

# Modify OTU table name as the output name for the table
OTU_table_top_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_top_", min_abundance, "percent.tsv", sep = "")

# Write the table
write.table(OTU_table_top, file = OTU_table_top_name, sep = "\t", row.names = FALSE, col.names = TRUE)

## Get some final output stats for the user
original_OTUs <- nrow(OTU_table)
top_OTUs <- length(top_OTUs)

num_samples <- length(col_sums)
hits_kept_overall <- as.numeric(colSums(OTU_table_top[,c(2:(num_samples + 1))]))
percent_kept_overall <- round(hits_kept_overall / col_sums * 100, 1)
percent_kept_overall_min <- min(percent_kept_overall)
percent_kept_overall_max <- max(percent_kept_overall)
percent_kept_overall_mean <- round(mean(percent_kept_overall), 1)

log_text <- paste("Output stats:\n\nNumber of samples: ", num_samples, 
            "\nNumber of OTUs in original OTU table: ", original_OTUs,
            "\nNumber of OTUs in subset OTU table (>= ", min_abundance, " percent): ", top_OTUs,
            "\n\nOTUs kept in subset OTU table represent between ", percent_kept_overall_min,
            " and ", percent_kept_overall_max, " of total reads for each sample (", 
            percent_kept_overall_mean, " percent on average)\n",
            sep = "")

# Modify OTU table name as the output name for the log
log_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_subset_log.txt", sep = "")

# Write the table
write(log_text, file = log_name, ncolumns = 1)

# Make a data frame of sample read counts and print for the user, if desired
if (print_total_read_counts == TRUE) {
  # Get percent hits kept for each sample with given abundance criterion
  percentage_kept <- unlist(lapply(1:length(OTU_select_prop), function(x) { sum(OTU_select_prop[[x]]$Percentage_abundance) } ))
  # Remove unused samples if needed
  if (length(removed_objects) > 0) {
    col_sums_forLog <- col_sums[-removed_objects_num]
  } else {
    col_sums_forLog <- col_sums
  }
  hits_kept <- percentage_kept / 100 * col_sums_forLog

  # Get percent hits kept for each sample in the actual OTU table (where abundant OTUs from one sample are also included for the others)
  num_samples <- length(col_sums_forLog)
  hits_kept_overall <- as.numeric(colSums(OTU_table_top[,c(2:(num_samples + 1))]))
  
  # Make data frame
  total_read_counts <- data.frame("Sample" = names(OTU_select_prop), "Read count before subset" = col_sums_forLog, 
                                  "Hits with abundance over threshold in sample" = hits_kept,
                                  "Hits kept in final OTU table" = hits_kept_overall)
  
  # Modify OTU table name as the output name for the table
  total_read_counts_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_initial_read_counts.tsv", sep = "")
  
  # Write the table
  write.table(total_read_counts, file = total_read_counts_name, sep = "\t", row.names = FALSE, col.names = TRUE)
}

