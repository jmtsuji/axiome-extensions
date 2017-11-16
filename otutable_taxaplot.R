# Taxa plot from OTU table
# Adapted by Jackson M. Tsuji from BubblePlot script by Sara Victoria Coyotzi Alcaraz (Neufeld Lab)
# Description: Subsets OTUs from an OTU table with relatve abundance (in at least one sample) of specified value (e.g. >= 1%, or top 10).
#               Works for rarefied and non-rarefied OTU tables. Then plots nicely.
# Version 3.0.2
# Last updated: Apr. 21, 2017 (by Jackson Tsuji)
# Required packages: plyr, dplyr, ggplot2, grid. If not installed, run install.packages("plyr"), and so on.

#####################################################
## Script settings: #################################
## Basic settings:
# Working directory
setwd("/Users/JTsuji/Documents/Research_General/PhD/04e_Biogeography/06_iTag_analysis/02_Sep2017_ELA_CyaMig_AXIOME/CyaMig_01/misc_analysis/mesas_representative_sequences0/")
# OTU table
otutable_name <- "otu_table_with_seq.tab"
# Amount to subset (top 'x' if greater than 1, or proportion if between 0-1)
subset_value <- 0.01
# Print taxa plot of top organisms in OTU table?
print_taxa_plot <- TRUE
# Taxonomic rank to display in legend?
plotting_rank <- "Family"
# Classification method ("greengenes" or "silva")
classification_method <- "greengenes"
# Optional: print list of top OTUs per sample?
print_top_otu_list <- TRUE
# Optional: print list of total read counts for each sample before subset?
print_total_read_counts <- FALSE
# Optional: print subset OTU table with only top OTUs?
print_subset_otu_table <- FALSE
# Optional: print log of subsetting?
print_subset_log <- FALSE
#####################################################
#####################################################

#####################################################
## Advanced tools for making designer plots #########
#####################################################
# Print list of taxon ranks with functional group and colour template in the ggplot for manual sorting?
print_plot_ranks <- FALSE

# Import taxa order with functional group names and colours (to sort ggplot rows)? Can use the taxon rank list that can be printed using the function above as a template.
import_taxa_order <- FALSE
plotting_info_file_name <- "plotting_info_top10_family.tsv"

# Make custom plot using imported metadata file?
# Need to custom-modify the plot code to meet your needs. Example is provided. See section 7 of code (line 421 onward)
# Need to also specify plotting_info_file_name above for this to work.
custom_plot_with_metadata <- FALSE
metadata_file_name <- "metadata_taxaplot.tsv"
#####################################################
#####################################################


#####################################################
################# Begin script ######################
#####################################################

#####################################################
#### 1. Import and pre-processing data ##############
#####################################################

# Import OTU table
OTU_table <- read.table(otutable_name, header = TRUE, sep="\t", comment.char = "", skip = 1, stringsAsFactors = FALSE)

# Standardize OTU ID header
colnames(OTU_table)[1] <- "OTU_ID"

# Make OTU ID character format, for consistency later
OTU_table$OTU_ID <- as.character(OTU_table$OTU_ID)

# Check if representative sequence was included, and exclude if present, but keep a backup
if (is.na(match("ReprSequence", colnames(OTU_table))) == FALSE) {
  ReprSequence_exists <- TRUE
  print("Found ReprSequence column. Will preserve in final tables.")
  ReprSequence_col_num <- match("ReprSequence", colnames(OTU_table))
  OTU_table_ReprSequence <- OTU_table[,c(1,ReprSequence_col_num)]
  OTU_table <- OTU_table[,-(ReprSequence_col_num)]
} else {
  ReprSequence_exists <- FALSE
}

# Same for Consensus Lineage, but keep a backup
if (is.na(match("Consensus.Lineage", colnames(OTU_table))) == FALSE) {
  Consensus.Lineage_exists <- TRUE
  print("Found Consensus.Lineage column. Will preserve in final tables.")
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
read_count <- as.numeric(colSums(OTU_table))
for (i in c(1:length(read_count))){
    OTU_proportion[,i] <- as.numeric(OTU_proportion[,i]/read_count[i])
}
# Change proportions to percentages
OTU_proportion <- OTU_proportion * 100


#####################################################
#### 2. Subset for top OTUs #########################
#####################################################

# Determine if subsetting above proportion threshold or rank
if (subset_value > 1) {
  print(paste("Provided rank value. Will subset top ", subset_value, ".", sep = ""))
  subset_type <- "rank"
} else if (subset_value < 1 & subset_value > 0) {
  print(paste("Provided propotion value. Will subset OTUs above ", subset_value * 100, "% abundance.", sep = ""))
  subset_type <- "proportion"
  # Convert subset value to percentage for ease later on
  subset_value <- subset_value * 100
} else {
  stop("Error: subset_type provided is not a proportion or rank. Did you input a negative number?")
}

# Subset OTUs for each sample and keep only those as a list of data frames (different length for each data frame)
if (subset_type == "rank") {
  # Make initial list separating each column
  OTU_subset <- lapply(1:ncol(OTU_proportion), function(x) { as.data.frame(OTU_proportion[,x, drop = FALSE]) } )
  # Drop function idea from http://stackoverflow.com/a/29325692 (accessed Apr. 6, 2017)
  
  # Make function for sorting top ranks
  rank_sort <- function(sample_position) {
    OTU_subset[[sample_position]] <- OTU_subset[[sample_position]][order(OTU_subset[[sample_position]], decreasing = TRUE),, drop = FALSE]
    OTU_subset[[sample_position]] <- OTU_subset[[sample_position]][c(1:subset_value),, drop = FALSE]
  }
  
  # Apply the function via lapply
  OTU_subset <- lapply(1:length(OTU_subset), rank_sort)
  
} else if (subset_type == "proportion") {
  OTU_subset <- lapply(1:ncol(OTU_proportion), function(x) { as.data.frame(subset(OTU_proportion, OTU_proportion[,x] >= subset_value, select=x)) } )
}

# Assign the sample names and OTU IDs as columns in the data frame
OTU_subset <- lapply(1:length(OTU_subset), function(x) { OTU_subset[[x]]$SampleName <- colnames(OTU_subset[[x]]); OTU_subset[[x]]$OTU_ID <- rownames(OTU_subset[[x]]); OTU_subset[[x]] } )

# Name the column for proportion as Percentage_abundance
OTU_subset <- lapply(1:length(OTU_subset), function(x) { colnames(OTU_subset[[x]])[1] = "Percentage_abundance"; OTU_subset[[x]] } )

# Give names to the list variables, which are of data frame type 
names(OTU_subset) <- colnames(OTU_proportion)

# From list of data frames to a single data frame
OTU_subsetDF <- do.call(rbind.data.frame, OTU_subset)

# Eliminate old row names
rownames(OTU_subsetDF) <- NULL

# Re-order columns for clarity
OTU_subsetDF <- OTU_subsetDF[,c(3,2,1)]

# Make OTU ID character for later
OTU_subsetDF$OTU_ID <- as.character(OTU_subsetDF$OTU_ID)

# Add consensus lineage and representative sequences back on
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(dplyr)))

if (Consensus.Lineage_exists == TRUE) {
  OTU_subsetDF <- dplyr::left_join(OTU_subsetDF, OTU_table_lineage, by = "OTU_ID")
}
if (ReprSequence_exists == TRUE) {
  OTU_subsetDF <- dplyr::left_join(OTU_subsetDF, OTU_table_ReprSequence, by = "OTU_ID")
}


#####################################################
#### 3. Parse consensus lineage values ##############
#####################################################

if (Consensus.Lineage_exists == FALSE) {
  print("Cannot make a taxaplot because no Consensus.Lineage was provided in input OTU table")
} else {
  # Drop old leftover levels in data frame, if a factor (shouldn't be anymore as of version 3.0.2)
  if (is.factor(OTU_subsetDF$Consensus.Lineage)) {
    OTU_subsetDF$Consensus.Lineage <- droplevels(OTU_subsetDF$Consensus.Lineage)
  }
  
  # Convert Consensus.Lineage to character (should already be as of version 3.0.2)
  OTU_subsetDF$Consensus.Lineage <- as.character(OTU_subsetDF$Consensus.Lineage)
  
  # Parse Consensus.Lineage
  Consensus.Lineage_list <- strsplit(OTU_subsetDF$Consensus.Lineage, "; ")
  
  Consensus.Lineage_list <- lapply(1:length(Consensus.Lineage_list), function(x) { as.data.frame(Consensus.Lineage_list[[x]], stringsAsFactors = FALSE) } )
  names(Consensus.Lineage_list) <- as.character(OTU_subsetDF$OTU_ID)
  
  # Loop to complete final steps of list manipulation (lapply was not cooperating for some reason)
  for (list_number in 1:length(Consensus.Lineage_list)) {
    # Change column name to OTU_ID for record-keeping
    colnames(Consensus.Lineage_list[[list_number]]) <- as.character(OTU_subsetDF$OTU_ID[list_number])
    
    # Parse taxonomic names to remove prefix based on classification method
    if (classification_method == "greengenes") {
      for (i in 1:nrow(Consensus.Lineage_list[[list_number]])) {
        # Make greengenes row vector for later
        greengenes_blank_rows <- numeric()
        
        # Only proceed if prefix appears to exist for that rank
        if (is.na(match(substr(Consensus.Lineage_list[[list_number]][i,1], 1, 3), c("k__","p__","c__","o__","f__","g__","s__"))) == FALSE) {
          
          # Get length of rank name
          ranklength <- nchar(Consensus.Lineage_list[[list_number]][i,1])
          
          # Remove first three characters of rank name (hard-coded number for now)
          Consensus.Lineage_list[[list_number]][i,1] <- substr(Consensus.Lineage_list[[list_number]][i,1], 4, ranklength)
          
        }
      }
      
      # Delete any blank greengenes rows
      Consensus.Lineage_list[[list_number]] <- Consensus.Lineage_list[[list_number]][Consensus.Lineage_list[[list_number]][,1] != "",, drop = FALSE]
      # drop = FALSE idea from http://stackoverflow.com/a/21025639, accessed 170421
      # row criterion idea from http://stackoverflow.com/a/6437778, accessed 170421
      
    } else if (classification_method == "silva") {
      for (i in 1:nrow(Consensus.Lineage_list[[list_number]])) {
        # Only proceed if prefix appears to exist for that rank
        if (substr(Consensus.Lineage_list[[list_number]][i,1], 1, 2) == "D_") {
          # Get length of rank name
          ranklength <- nchar(Consensus.Lineage_list[[list_number]][i,1])
          # Remove first five characters of rank name (hard-coded number for now)
          Consensus.Lineage_list[[list_number]][i,1] <- substr(Consensus.Lineage_list[[list_number]][i,1], 6, ranklength)
        }
      }
    } else {
      stop("Improper classification method: should be 'greengenes' or 'silva'.")
    }
    
    ## Add info from preceeding taxonomic rank to ones with ambiguous classification
    # General function
    # Info about function: requires that Consensus.Lineage_list and list_number exist. Works for just one object in a list as specified by list_number.
    #                     Output is a revised object to replace the one in the original list.
    #                     If there were no search hits, the output object is identical to the original object.
    # To improve: currently, will only find the ambiguous term of highest taxonomic rank and replace. Any of lower taxonomic rank will be left as-is.
    #             The function currently needs to be run multiple times if you have an instance of repeated ambiguous terms.
    fill_in_info <- function(search_term) {
      # Check if any entries match the search term
      if (is.na(match(search_term, Consensus.Lineage_list[[list_number]][,1])) == FALSE) {
        # If so, record their row number
        unc_row_num <- match(search_term, Consensus.Lineage_list[[list_number]][,1])
        
        # In case multiple entries match, find the one of highest taxonomic rank
        unc_row_num_min <- min(unc_row_num) # Actually, this no longer works, because 'match' function is used. Will only output the first hit anyway.
        
        # Append the searched term with info from the preceeding taxonomic rank
        new_names <- paste(Consensus.Lineage_list[[list_number]][unc_row_num,1], "_", Consensus.Lineage_list[[list_number]][(unc_row_num_min - 1),1], sep = "")
        
        # Create a copy of original list component updated with new names
        new_names_full <- Consensus.Lineage_list[[list_number]]
        new_names_full[unc_row_num,1] <- new_names
        
        # Return new list component to user
        return(new_names_full)
      } else {
        # If no matches, then return unmodified original list
        return(Consensus.Lineage_list[[list_number]])
      }
    }
    
    ## Run function for search terms of interest, then add new list component to the original list:
    Consensus.Lineage_list[[list_number]] <- fill_in_info("uncultured bacterium")
    Consensus.Lineage_list[[list_number]] <- fill_in_info("uncultured")
    Consensus.Lineage_list[[list_number]] <- fill_in_info("Ambiguous_taxa")
    
    # Store number of rows to see if less than 7 (full taxonomic rank)
    row_num <- nrow(Consensus.Lineage_list[[list_number]])
    if (row_num < 7) {
      # Fill out data frames shorter than 7 rows. Append "unclassified" to last taxonomic rank available.
      Consensus.Lineage_list[[list_number]][((row_num + 1):(7)),1] <- paste("Unclassified_", Consensus.Lineage_list[[list_number]][row_num,1], sep = "")
    }
    
    # Transpose data frames in preparation for recombining
    Consensus.Lineage_list[[list_number]] <- as.data.frame(t(Consensus.Lineage_list[[list_number]]), stringsAsFactors = FALSE)
    
    # Assign taxonomic names to rank columns
    colnames(Consensus.Lineage_list[[list_number]]) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  
  # Convert classifications to data frame
  Consensus.Lineage_df <- dplyr::bind_rows(Consensus.Lineage_list)
  # Consensus.Lineage_df$OTU_ID <- OTU_subsetDF$OTU_ID
  # OTU_subsetDF <- dplyr::left_join(OTU_subsetDF, Consensus.Lineage_df, by = "OTU_ID")
  
  # Bind classifications to original dataframe
  OTU_subsetDF <- cbind(OTU_subsetDF, Consensus.Lineage_df)
  
  # Re-order for clarity
  OTU_subsetDF <- OTU_subsetDF[,c(1:3,6:12,4:5)]
}

#####################################################
#### 4. Print plotting table ########################
#####################################################
# Print above table for reference for the user, if desired
if (print_top_otu_list == TRUE) {
  # Modify OTU table name as the output name for the table
  top_otu_list_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_top_otu_names.tsv", sep = "")
  
  # Write the table
  write.table(OTU_subsetDF, file = top_otu_list_name, sep = "\t", row.names = FALSE, col.names = TRUE)
}

#####################################################
#### 5. Print basic taxa plot #######################
#####################################################
## Make a taxa plot from subset, if desired
if (print_taxa_plot == TRUE) {
  # Load plotting packages
  suppressWarnings(suppressMessages(library(ggplot2)))
  library(grid)
  
  # Check Consensus.Lineage exists
  if (Consensus.Lineage_exists == FALSE) {
    print("Cannot make a taxaplot because no Consensus.Lineage was provided in input OTU table")
  } else { 
    # Choose rank name to plot by
    #plotting_rank <- "Family" # Moved to start variables
    # Find rank in data frame
    plotting_rank_col <- match(plotting_rank, colnames(OTU_subsetDF))
    
    # Get length of legend
    legend_length <- length(unique(OTU_subsetDF[,plotting_rank_col]))
    
    # Decide on number of legend columns as a ratio of legend length (30 to a column)
    legend_ncol <- ceiling(legend_length / 30)
    
    # Plot the Data (old version, then new version specific for paper)
    taxa_plot <- ggplot(OTU_subsetDF, aes(SampleName)) +
      geom_bar(aes_string(weight = "Percentage_abundance", fill = plotting_rank), position = "stack") +
      # geom_bar(aes(weight = Percentage_abundance, fill = Family), position = "stack") +
      #scale_fill_manual() +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 7), 
            panel.border = element_rect(colour = "transparent"),
            axis.text = element_text(size = 5, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
            legend.text = element_text(size = 5, face = "italic"), legend.title = element_blank(),
            legend.key = element_rect(colour = "transparent"), legend.key.size = unit(3, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", legend.position = "right", plot.margin = unit(c(2,2,2,2), "mm")) +
      guides(fill = guide_legend(ncol = legend_ncol)) +
      xlab("Sample") +
      ylab("Cumulative percentage abundance")
      
    print(taxa_plot)
    
    output_PDF_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_taxa_plot_auto.pdf", sep = "")
    
    # Save the plot
    ggsave(output_PDF_name, width = 150, height = 100, units = c("mm"))
  }
}

##############################################################################
#### 6. Print advanced plot (using advanced functions) #######################
##############################################################################
# Output template file for making plot with user customization
if (print_plot_ranks == TRUE) {
  print("As specified in advanced settings, printing taxon ranks template for use with advanced plotting.")
  
  # Check that Consensus.Lineage exists
  if (Consensus.Lineage_exists == FALSE) {
    print("Cannot output template file for advanced plotting because no Consensus.Lineage was provided in input OTU table")
  } else {
    plotted_ranks <- unique(OTU_subsetDF[,plotting_rank_col])
    ranks_template_table <- data.frame("Ranks" = plotted_ranks, 
                                       "Functional_group" = character(length = length(plotted_ranks)),
                                       "HTML_colour_code" = character(length = length(plotted_ranks)), stringsAsFactors = FALSE)
    colnames(ranks_template_table)[1] <- plotting_rank
    
    # Knowing that taxonomic rank rows begin with Kingdom and move to Species... find where Kingdom begins as a reference
    kingdom_rank_col <- match("Kingdom", colnames(OTU_subsetDF))
    
    # Find one row in original table matching each unique rank value. Any row matching the rank should have the same preceeding ranks, in theory.
    plotted_ranks_rows <- match(plotted_ranks, OTU_subsetDF[,plotting_rank_col])
    
    # Combine row and column info to extract the ranks preceeding the rank of interest from the original table
    preceeding_ranks <- OTU_subsetDF[plotted_ranks_rows, c(kingdom_rank_col:plotting_rank_col)]
    
    # Join tables
    ranks_template_table <- dplyr::left_join(ranks_template_table, preceeding_ranks, by = plotting_rank)
    
    # Modify OTU table name as the output name for the table
    output_template_table_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_template_plotting_info.tsv", sep = "")
    
    # Save the template table
    write.table(ranks_template_table, output_template_table_name, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}


# Making advanced plot by loading user-modified template file
if (import_taxa_order == TRUE) {
  # Check that Consensus.Lineage exists
  if (Consensus.Lineage_exists == FALSE) {
    print("Cannot make advanced plot because no Consensus.Lineage was provided in input OTU table")
  } else {
    # Import the table
    plotting_info_table <- read.table(plotting_info_file_name, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    
    # Subset table to remove extra ranks info
    plotting_info_table_subset <- plotting_info_table[,1:3]
    
    # Merge with main plotting table
    suppressWarnings(suppressMessages(library(plyr)))
    suppressWarnings(suppressMessages(library(dplyr)))
    OTU_subsetDF_advanced <- dplyr::left_join(OTU_subsetDF, plotting_info_table_subset, by = plotting_rank)
    
    # Apply sort order
    OTU_subsetDF_advanced[,plotting_rank_col] <- factor(OTU_subsetDF_advanced[,plotting_rank_col], levels = plotting_info_table_subset[,1], ordered = TRUE)
    
    # Get colour order
    plotting_colours <- plotting_info_table_subset[,3]
    
    # Plot
    suppressWarnings(suppressMessages(library(ggplot2)))
    library(grid)
    
    taxa_plot_advanced <- ggplot(OTU_subsetDF_advanced, aes(SampleName)) +
      geom_bar(aes_string(weight = "Percentage_abundance", fill = plotting_rank), position = "stack") +
      # geom_bar(aes(weight = Percentage_abundance, fill = Family), position = "stack") +
      scale_fill_manual(values = plotting_colours) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 7), 
            panel.border = element_rect(colour = "transparent"),
            axis.text = element_text(size = 5, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
            legend.text = element_text(size = 5, face = "italic"), legend.title = element_blank(),
            legend.key = element_rect(colour = "transparent"), legend.key.size = unit(3, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", legend.position = "right", plot.margin = unit(c(2,2,2,2), "mm")) +
      guides(fill = guide_legend(ncol = legend_ncol)) +
      xlab("Sample") +
      ylab("Cumulative percentage abundance")
    
    print(taxa_plot_advanced)
    
    # Modify OTU table name as the output name for the plot
    version <- "1.0a"
    output_PDF_name_advanced <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_taxa_plot_advanced_", version, ".pdf", sep = "")
    
    # Save the plot
    ggsave(output_PDF_name_advanced, width = 150, height = 100, units = c("mm"))
  }
}
    

##############################################################################
#### 7. Print designer plot (using advanced functions)                ########
####      **Requires custom modification of the code here by the user ########
##############################################################################

if (custom_plot_with_metadata == TRUE) {
  # Check that Consensus.Lineage exists
  if (Consensus.Lineage_exists == FALSE) {
    print("Cannot generate custom plot because no Consensus.Lineage was provided in input OTU table")
  } else {
  
    # Check if advanced plotting is already on
    if (import_taxa_order == FALSE) {
      print("Making designer plot (custom_plot_with_metadata): import_taxa_order is set to FALSE, but will still import order from plotting_info_file_name as required...")
      
      # Import the table
      plotting_info_table <- read.table(plotting_info_file_name, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
      
      # Subset table to remove extra ranks info
      plotting_info_table_subset <- plotting_info_table[,1:3]
      
      # Merge with main plotting table
      # Here and elsewere: got the suppress warnings idea from http://stackoverflow.com/a/1893172, accessed April 21, 2017
      suppressWarnings(suppressMessages(library(plyr)))
      suppressWarnings(suppressMessages(library(dplyr)))
      OTU_subsetDF_advanced <- dplyr::left_join(OTU_subsetDF, plotting_info_table_subset, by = plotting_rank)
      
      # Apply sort order
      OTU_subsetDF_advanced[,plotting_rank_col] <- factor(OTU_subsetDF_advanced[,plotting_rank_col], levels = plotting_info_table_subset[,1], ordered = TRUE)
      
      # Get colour order
      plotting_colours <- plotting_info_table_subset[,3]
    }
    
    # Load metadata
    metadata <- read.table(metadata_file_name, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    
    # Fix Sample_ID column
    colnames(metadata)[1] <- "SampleName"
    
    # Keep columns of interest in subset
    metadata_subset <- metadata[,c(1,9:12)]
    
    # Merge with main plotting table
    suppressWarnings(suppressMessages(library(plyr)))
    suppressWarnings(suppressMessages(library(dplyr)))
    OTU_subsetDF_custom <- dplyr::left_join(OTU_subsetDF_advanced, metadata_subset, by = "SampleName")
    
    ### Custom changes ####
    # Set depth to factor (already sorted correctly, because was numeric before)
    OTU_subsetDF_custom$Depth_m <- factor(OTU_subsetDF_custom$Depth_m)
    
    # Add order to lake numbers
    OTU_subsetDF_custom$Lake <- factor(OTU_subsetDF_custom$Lake, 
                                       levels = c("L239", "L626", "L373", "L224", "L222", "L221", "L304", "L442", "L227"), ordered = TRUE)
    
    # Plot
    suppressWarnings(suppressMessages(library(ggplot2)))
    library(grid)
    
    taxa_plot_custom <- ggplot(OTU_subsetDF_custom, aes(Depth_m)) +
      geom_bar(aes_string(weight = "Percentage_abundance", fill = plotting_rank), position = "stack", width = 0.8) +
      facet_grid(Sampling_season ~ Lake, scales = "free_x", switch = "y", space = "free_x") +
      scale_fill_manual(values = plotting_colours) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 7), 
            panel.border = element_rect(colour = "black", size = 0.3), panel.spacing.x = unit(0.5, "mm"),
            strip.text = element_text(size = 6), strip.background = element_rect(fill = "#e6e6e6", size = 0.3),
            axis.text = element_text(size = 5, colour = "black"),
            axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.2),
            legend.text = element_text(size = 5, face = "italic"), legend.title = element_blank(),
            legend.key = element_rect(colour = "transparent"), legend.key.size = unit(3, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", legend.position = "right", plot.margin = unit(c(2,2,2,2), "mm")) +
      guides(fill = guide_legend(ncol = legend_ncol)) +
      xlab("Depth (m)") +
      ylab("Cumulative percentage abundance")
    
    print(taxa_plot_custom)
    
    # Modify OTU table name as the output name for the plot
    version <- "1.0a"
    output_PDF_name_custom <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_taxa_plot_custom_", version, ".pdf", sep = "")
    
    # Save the plot
    ggsave(output_PDF_name_custom, width = 150, height = 100, units = c("mm"))
  }
}

#####################################################
#### 8. Print subset OTU table, if desired ##########
#####################################################
if (print_subset_otu_table == TRUE) {
  ## Make subset OTU table
  # Get top OTU IDs
  top_OTUs <- unique(OTU_subsetDF$OTU_ID)
  
  # Make a data frame
  OTU_table_top <- data.frame("OTU_ID" = top_OTUs, stringsAsFactors = FALSE)
  
  # Match OTU IDs to original OTU table
  suppressWarnings(suppressMessages(library(plyr)))
  suppressWarnings(suppressMessages(library(dplyr)))
  
  # Make copy of OTU table to modify
  OTU_table_copy <- OTU_table
  
  # Add OTU_ID column back to OTU table and delete row names
  OTU_table_copy$OTU_ID <- as.character(rownames(OTU_table_copy))
  rownames(OTU_table_copy) <- NULL
  
  OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_copy, by = "OTU_ID")
  
  # Add consensus lineage and representative sequences back on
  if (Consensus.Lineage_exists == TRUE) {
    OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_lineage, by = "OTU_ID")
  }
  if (ReprSequence_exists == TRUE) {
    OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_ReprSequence, by = "OTU_ID")
  }
  
  # Modify OTU table name as the output name for the table
  OTU_table_top_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_top_", subset_value, "percent.tsv", sep = "")
  
  # Write the table
  write.table(OTU_table_top, file = OTU_table_top_name, sep = "\t", row.names = FALSE, col.names = TRUE)
}

#####################################################
#### 9. Print subset log, if desired ################
#####################################################
if (print_subset_log == TRUE) {
  # If print_subset_otu_table was not turned on, need to run some necessary code for log file (duplicated from above section):
  if (print_subset_otu_table == FALSE) {
    ## Make subset OTU table
    # Get top OTU IDs
    top_OTUs <- unique(OTU_subsetDF$OTU_ID)
    
    # Make a data frame
    OTU_table_top <- data.frame("OTU_ID" = top_OTUs, stringsAsFactors = FALSE)
    
    # Match OTU IDs to original OTU table
    suppressWarnings(suppressMessages(library(plyr)))
    suppressWarnings(suppressMessages(library(dplyr)))
    
    # Make copy of OTU table to modify
    OTU_table_copy <- OTU_table
    
    # Add OTU_ID column back to OTU table and delete row names
    OTU_table_copy$OTU_ID <- as.character(rownames(OTU_table_copy))
    rownames(OTU_table_copy) <- NULL
    
    OTU_table_top <- dplyr::left_join(OTU_table_top, OTU_table_copy, by = "OTU_ID")
  }
    
  ## Get some final output stats for the user
  original_OTUs_num <- nrow(OTU_table)
  top_OTUs_num <- length(top_OTUs)
  
  num_samples <- length(read_count)
  hits_kept_overall <- as.numeric(colSums(OTU_table_top[,c(2:(num_samples + 1))]))
  percent_kept_overall <- round(hits_kept_overall / read_count * 100, 1)
  percent_kept_overall_min <- min(percent_kept_overall)
  percent_kept_overall_max <- max(percent_kept_overall)
  percent_kept_overall_mean <- round(mean(percent_kept_overall), 1)
  
  # Make different log text depending on subset method chosen:
  if (subset_type == "rank") {
    log_text <- paste("Output stats:\n\nNumber of samples: ", num_samples, 
                      "\nNumber of OTUs in original OTU table: ", original_OTUs_num,
                      "\nNumber of OTUs in subset OTU table (top ", subset_value, " for each sample): ", top_OTUs_num,
                      "\n\nOTUs kept in subset OTU table represent between ", percent_kept_overall_min,
                      " and ", percent_kept_overall_max, " of total reads for each sample (", 
                      percent_kept_overall_mean, " percent on average)\n",
                      sep = "")
  } else if (subset_type == "proportion") {
    log_text <- paste("Output stats:\n\nNumber of samples: ", num_samples, 
                      "\nNumber of OTUs in original OTU table: ", original_OTUs_num,
                      "\nNumber of OTUs in subset OTU table (>= ", subset_value, " percent): ", top_OTUs_num,
                      "\n\nOTUs kept in subset OTU table represent between ", percent_kept_overall_min,
                      " and ", percent_kept_overall_max, " of total reads for each sample (", 
                      percent_kept_overall_mean, " percent on average)\n",
                      sep = "")
  }
    
  # Modify OTU table name as the output name for the log
  log_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_subset_log.txt", sep = "")
  
  # Write the table
  write(log_text, file = log_name, ncolumns = 1)
}
  

#####################################################
#### 10. Print subset read stats, if desired ########
#####################################################
# Make a data frame of sample read counts and print for the user, if desired
if (print_total_read_counts == TRUE) {
  # Get percent hits kept for each sample with given abundance criterion
  percentage_kept <- unlist(lapply(1:length(OTU_subset), function(x) { sum(OTU_subset[[x]]$Percentage_abundance) } ))
  hits_kept <- percentage_kept / 100 * read_count
  
  # If desired, get percent hits kept for each sample in the actual OTU table (where abundant OTUs from one sample are also included for the others)
  if (print_subset_otu_table == TRUE) {
    num_samples <- length(read_count)
    hits_kept_overall <- as.numeric(colSums(OTU_table_top[,c(2:(num_samples + 1))]))
    
    # Make data frame
    total_read_counts <- data.frame("Sample" = names(OTU_subset), "Read count before subset" = read_count, 
                                    "Hits with abundance over threshold in sample" = hits_kept,
                                    "Hits kept in final OTU table" = hits_kept_overall)
  } else {
    total_read_counts <- data.frame("Sample" = names(OTU_subset), "Read count before subset" = read_count, 
                                    "Hits with abundance over threshold in sample" = hits_kept)
  }
    
  # Modify OTU table name as the output name for the table
  total_read_counts_name <- paste(substr(otutable_name, 1, nchar(otutable_name)-4), "_initial_read_counts.tsv", sep = "")
  
  # Write the table
  write.table(total_read_counts, file = total_read_counts_name, sep = "\t", row.names = FALSE, col.names = TRUE)
}
