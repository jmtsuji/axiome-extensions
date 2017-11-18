#!/usr/bin/env Rscript
# Enhanced mesas-pcoa.R script from "Neufeld MESaS"
# Copyright Jackson M. Tsuji, 2017
# Neufeld lab, University of Waterloo, Canada

########################################################
########################################################
# Variables for local machine use - USER INPUT HERE (only)
# After adding desired input, run the script for output
setwd("/Users/JTsuji/Documents/Research_General/PhD/04e_Biogeography/06_iTag_analysis/01_ELA2016_survey_AXIOME/02_round1_all/convert_otu_table/rarefied_biom_to_tab0/02_noStreams/")
outName <- "pcoa"
otuTable <- "otu_table_auto_noStreams.tab"
mappingFile <- "171116_metadata_ELA2016_noStreams.tsv"

# Method for calulating distance matrix
dmethod <- "bray"

# Confidence threshold (e.g., 0.95) for data ellipses. Set to 0 to not plot ellipses.
ellipsoidConf <- 0
########################################################
########################################################

# Dump a text file in case of error, for debugging
options(error = quote(dump.frames("mesas-pcoa-debug", TRUE)))

# Function to automatically install package if missing
pkgTest <- function(x)
{
        if (!require(x,character.only = TRUE))
        {
                install.packages(x,dep=TRUE)
                if(!require(x,character.only = TRUE)) stop("Package not found")
        }
}

# Load libraries
pkgTest('getopt')
pkgTest("vegan")
pkgTest("ape")
pkgTest("scales")
pkgTest("ggplot2")
pkgTest("scales")
pkgTest("glue")
pkgTest("RColorBrewer")
if (!is.null(ellipsoidConf) ) {
  pkgTest("car")
}

# Make sure we have a valid distance method
dlist = c('manhattan','euclidean','canberra','bray','kulczynski','jaccard','gower','altGower','morisita','horn','mountford','raup','binomial','chao','cao')

if (! dmethod %in% dlist ) {
	stop(print(paste("Invalid distance method:", dmethod)))
} 

# Read in and format otu table
rawtable <- read.table(otuTable, skip = 1, comment.char = "", header = TRUE, row.names = 1, sep = "\t")
otutable.first <- t(rawtable[1:(ncol(rawtable) - 1)])

# Read metadata mapping file
mapping.first <- read.table(mappingFile, header = TRUE, comment.char = "", row.names = 1, sep = "\t") 

# Adjust numeric sample IDs in metadata file to match OTU table ("X" gets added during AXIOME)
numeric_names <- suppressWarnings(sapply(rownames(mapping.first),as.numeric))
numeric <- sum(sapply(numeric_names,is.na)) == 0
if (numeric) {
  rownames(mapping.first) <- paste("X", rownames(mapping.first), sep = "")
}
rownames(mapping.first) <- gsub("-",".",rownames(mapping.first),fixed=TRUE)

# Check if sample IDs are now the same
if (length(base::setdiff(rownames(mapping.first), rownames(otutable.first))) >= 1) {
  stop("Sample IDs do not match between metadata file and OTU table. Exiting...")
}

# Sort samples by sample ID
otutable <- otutable.first[order(rownames(otutable.first), na.last = TRUE, decreasing = FALSE),]
mapping <- mapping.first[order(rownames(mapping.first), na.last = TRUE, decreasing = FALSE),]

# Calculate distance matrix
print("Computing PCoA")
d <- vegan::vegdist(otutable, method = dmethod)

# Perform PCoA calculations
print("Making PCoA Plot")
p <- ape::pcoa(d)

# Print out plotting vectors for user's reference
p_print <- data.frame(Axis.1 = p$vectors[,1], Axis.2 = p$vectors[,2], Axis.3 = p$vectors[,3])
write.table(p_print, file = "plotting_values.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

### Prepare to join by sample_ID
p_vectors <- as.data.frame(p$vectors)[,c(1:2)]
p_vectors <- cbind("sample_ID" = rownames(p_vectors), p_vectors)
mapping_join <- cbind("sample_ID" = rownames(mapping), mapping)

# Check if sample IDs are the same
if (length(base::setdiff(p_vectors$sample_ID, mapping_join$sample_ID)) >= 1) {
  stop("Sample IDs do not match between PCoA output file and metadata file. Exiting...")
}

# Merge by sample_ID
plotting_table <- dplyr::left_join(p_vectors, mapping_join, by = "sample_ID")

# Function to make PCoA plot
make_plot <- function(plotting_data_table, mapping_col) {
  # # Test vars
  # plotting_data_table <- plotting_table
  # mapping_col <- 12
  
  md_name <- colnames(plotting_data_table)[mapping_col]
  metadata <- plotting_data_table[,mapping_col]
  fac.len <- length(unique(metadata))
  
  # Skip plotting if all values are identical in the metadata category
  if (fac.len == 1) {
    
    print(paste("Ignoring", md_name, "because all values are identical."))
    
  } else if (fac.len > 1) {
    
    # Make general plot
    ordination_general <- ggplot(plotting_data_table, aes(x = Axis.1, y = Axis.2)) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 14), 
            strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 8), strip.background = element_rect(fill = "#e6e6e6"),
            panel.border = element_rect(colour = "black", size = 2), 
            axis.text = element_text(size = 9, colour = "black"),
            axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
            legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
            legend.key = element_rect(colour = "transparent", size = 0.5), legend.key.size = unit(5, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
      ggtitle(paste("PCoA ordination: ", md_name, "\nMethod: '", dmethod, "'", sep="")) +
      xlab(paste("Axis 1 (", round(p$values[1,3]*100, digits = 1), "%)", sep = "")) +
      ylab(paste("Axis 2 (", round(p$values[2,3]*100, digits = 1), "%)", sep = ""))
      
    # Finish the plot differently if data is numeric vs. non-numeric
    if (is.numeric(metadata) == TRUE) {
      
      # If numeric, use a gradient scale
      gradient_palette <- rev(RColorBrewer::brewer.pal(name = "YlGnBu", n = 5))
      
      ordination_final <- ordination_general + 
        geom_point(aes_string(fill = md_name), shape = 21, size = 3.5, alpha = 0.8) +
        scale_fill_gradientn(colours = gradient_palette, na.value = "grey50") +
        # Can change legend title from NULL to name if desired
        guides(fill = guide_colourbar(title = NULL))
      
    } else {
      
      # If non-numeric, treat as a factor
      metadata <- factor(metadata)
      
      # Choose simpler colour palette if number of colours is low, and a gradient if higher
      if (fac.len == 2) {
        # Special case because RColorBrewer can give a minimum of 3 colours as output
        metadata_colours <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)[c(1,2)]
      } else if (fac.len <= 8) {
        metadata_colours <- RColorBrewer::brewer.pal(name = "Dark2", n = fac.len)
      } else if (fac.len <= 12) {
        metadata_colours <- RColorBrewer::brewer.pal(name = "Set3", n = fac.len)
      } else {
        # Can adjust start and end below to change spectrum of colours plotted
        metadata_colours <- hue_pal(h = c(20,290))(fac.len)
      }
      
      ordination_final <- ordination_general + 
        geom_point(aes_string(fill = md_name), shape = 21, size = 3.5, alpha = 0.8) +
        scale_fill_manual(values = metadata_colours, na.value="white") +
        scale_colour_manual(values = metadata_colours, na.value="white") +
        # Can change legend title from NULL to name if desired
        guides(fill = guide_legend(title=NULL))
      
    }
    
    if (!is.null(ellipsoidConf) && ellipsoidConf != 0) {
      # ellipse <- car::dataEllipse(plotting_data_table$Axis.1, plotting_data_table$Axis.2, 
      #                             groups = metadata, levels = as.numeric(ellipsoidConf), add=TRUE, 
      #                             plot.points=FALSE, grid=FALSE, center.pch=FALSE, col=rainbow(fac.len))
      
      ordination_final <- ordination_final +
        stat_ellipse(aes_string(colour = md_name), level = ellipsoidConf, na.rm = TRUE, show.legend = FALSE,
                     alpha = 0.6, size = 0.75)
    }
    
    return(ordination_final)
    
  } else {
    
    # If fac.len is less than 1, it means something is wrong with the metadata.
    stop(paste("Something is wrong with metadata category, '", md_name, "'. Appears to contain < 1 unique value. Exiting...", sep = ""))
    
  }

}

# Make plots
plot_list <- lapply(4:ncol(plotting_table), function(x) { make_plot(plotting_table, x) })

# Remove ones where all values are identical
plots_to_keep <- logical(length = 0)
for (i in 1:length(plot_list)) {
  if (is.character(plot_list[[i]])) {
    plots_to_keep <- append(plots_to_keep, FALSE)
  } else {
    plots_to_keep <- append(plots_to_keep, TRUE)
  }
}
plot_list <- subset(plot_list, plots_to_keep)

# Restore mapping file to its original state when the PCoA was computed
# mapping <- mapping.first[order(rownames(mapping.first), na.last = TRUE, decreasing = FALSE),]

print("Making biplot")
# First, we take the mapping file and we coerce the columns to numeric
numeric_mapping <- suppressWarnings(apply(mapping,2,as.numeric))
numeric <- colSums(apply(numeric_mapping,2,is.na)) == 0
numeric_mapping <- numeric_mapping[,numeric, drop = FALSE]
rownames(numeric_mapping) <- rownames(mapping)
bi_plot <- biplot(p, apply(numeric_mapping, 2, scale, center=TRUE, scale=TRUE))

# Prepare Eigenvalues for output
eigenval_df <- data.frame("Eigenvalues" = p$values$Eigenvalues, "Relative_Eigenvalues" = p$values$Relative_eig)

# Write all output
pcoa_filename <- paste(dirname(outName),"/pcoa-", dmethod, ".pdf", sep="")
pdf(pcoa_filename, width = 8, height = 7)
print(plot_list)
dev.off()

biplot_filename <- paste(dirname(outName),"/pcoa-biplot.pdf", sep="")
pdf(biplot_filename)
print(biplot(p, apply(numeric_mapping, 2, scale, center=TRUE, scale=TRUE)))
dev.off()

eigenval_filename <- paste(dirname(outName),"/eigenvalues.tsv",sep="")
write.table(eigenval_df, eigenval_filename, col.names = TRUE, row.names = FALSE, sep = "\t")
