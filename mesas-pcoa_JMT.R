#!/usr/bin/env Rscript
# Based off mesas-pcoa.R script in "Neufeld MESaS"
# Created June 3, 2015 by Jackson Tsuji, to improve plot legend and axis titles, and to swtich to using ggplot
# Last Updated Sept 17, 2015
# Still in progress... **biplot NOT necessarily reliable!
# Version: 3.1.3

########################################################
########################################################
# Variables for local machine use - USER INPUT HERE (only)
# After adding desired input, run the script for output
setwd("/Users/JTsuji/Documents/Research_General/PhD/04e_Biogeography/06_iTag_analysis/01_ELA2016_survey_AXIOME/02_round1_all/convert_otu_table/rarefied_biom_to_tab0/02_noStreams/")
outName <- "pcoa"
otuTable <- "otu_table_auto_noStreams.tab"
mappingFile <- "171116_metadata_ELA2016_noStreams.tsv"

# Next two are currently set to defaults and often do not need to be changed
dmethod <- "bray"
ellipsoidConf <- 0
########################################################
########################################################

options(error = quote(dump.frames("mesas-pcoa-debug", TRUE)))

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

#Make sure we have a valid distance method
dlist = c('manhattan','euclidean','canberra','bray','kulczynski','jaccard','gower','altGower','morisita','horn','mountford','raup','binomial','chao','cao')

if (! dmethod %in% dlist ) {
	print(paste("Invalid distance method:", dmethod))
	q(status=1)
} 

#Read in and format otu table
print("Reading OTU table")

rawtable <- read.table(otuTable, skip = 1,
comment.char = "", header = TRUE, row.names = 1, sep = "\t")
otutable.first <- t(rawtable[1:(ncol(rawtable) - 1)])
#Sort numerically the samples
#otutable <- otutable[order(as.integer(sub("X","", rownames(otutable)))),]

print("Reading mapping")
mapping.first <- read.table(mappingFile, header = TRUE, comment.char = "", row.names = 1, sep = "\t") 
                            #colClasses = "character")
# For non-numeric sample names, comment out the following line
#colnames(mapping) <- paste("X", colnames(mapping), sep = "")
numeric_names <- suppressWarnings(sapply(rownames(mapping.first),as.numeric))
numeric <- sum(sapply(numeric_names,is.na)) == 0
if (numeric) {
  rownames(mapping.first) <- paste("X", rownames(mapping.first), sep = "")
}
rownames(mapping.first) <- gsub("-",".",rownames(mapping.first),fixed=TRUE)

print("Sorting tables")
# Sorting samples
otutable <- otutable.first[order(rownames(otutable.first), na.last = TRUE, decreasing = FALSE),]
mapping <- mapping.first[order(rownames(mapping.first), na.last = TRUE, decreasing = FALSE),]

print("Computing PCoA")

d <- vegdist(otutable, method=dmethod)

print("Making PCoA Plot")
p <- pcoa(d)

# Uncomment to reverse the orientations of Axis 1 or Axis 2 points, respectively (just for visualization purposes; does not affect the data)
# p$vectors[,1]=-p$vectors[,1]
# p$vectors[,2]=-p$vectors[,2]

# Preparing to sort by metadata category
p.sort <- p$vectors
p.sort <- cbind("sort"=c(1:nrow(p.sort)),p.sort)
mapping.sort = cbind("sort"=c(1:nrow(mapping)),mapping)
if (identical(row.names(p.sort), row.names(mapping.sort)) == FALSE) {
  stop(paste("Row names in sorting vectors do not match.\nOTU table: ", glue::collapse(head(row.names(p.sort), n = 3), sep = ", "), "...",
             "\n", "Mapping file: ", glue::collapse(head(row.names(mapping.sort), n = 3), sep = ", "), "...",
             "\n", "Results will be unreliable. Exiting...", sep = ""))
}
# Consider impoving this sanity check (above) to exit out of the script if it fails.

# Print out vectors for user
p.print <- data.frame(Axis.1 = p$vectors[,1], Axis.2 = p$vectors[,2])
write.table(p.print, file = "plotting_values.tsv", sep = "\t", col.names = NA)

make_plot <- function(mapping_col) {
  # Sorting by metadata category
  mapping.sort <- mapping.sort[order(mapping.sort[,mapping_col+1],na.last = TRUE,decreasing = FALSE),]
  mapping.sort[,1] = c(1:nrow(mapping.sort))
  mapping <- mapping.sort[,-1]
  mapping.sort <- mapping.sort[order(rownames(mapping.sort),na.last = TRUE,decreasing = FALSE),]
  
  p.sort <- p.sort[order(rownames(p.sort),na.last = TRUE,decreasing = FALSE),]
  p.sort[,1] <- mapping.sort[,1]
  p.sort <- p.sort[order(p.sort[,1],na.last = TRUE,decreasing = FALSE),]
  p$vectors <- p.sort[,-1]
  
  # Check row order names are identical
  if (identical(rownames(p$vectors), rownames(mapping)) == FALSE) {
    # Find specific non-matching entries
    non_matching_entries_num <- which(!(rownames(p$vectors) == rownames(mapping)))
      
    stop(paste("Metadata column ", mapping_col, ": row names of vegan output and mapping file do not match after sorting. See entries ", 
               head(glue::collapse(non_matching_entries_num, sep = ", ")), "...", sep = ""))
  }
  
  name <- colnames(mapping)[mapping_col]
  metadata <- mapping[,mapping_col]
  fac.len <- length(unique(metadata))
  
  if (fac.len < 2 ) {
    print(paste("Ignoring", name, "because all values are identical."))
  } else if (fac.len >= 2) {
   
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
    
    p.frame = data.frame("Axis.1" = p$vectors[,1], "Axis.2" = p$vectors[,2], "Metadata" = metadata)
    
    # General plot
    z = ggplot(p.frame, aes(x = Axis.1, y = Axis.2)) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 14), 
            strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 8), strip.background = element_rect(fill = "#e6e6e6"),
            panel.border = element_rect(colour = "black", size = 2), 
            axis.text = element_text(size = 9, colour = "black"),
            axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
            legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
            legend.key = element_rect(colour = "transparent", size = 0.5), legend.key.size = unit(5, "mm"),
            legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
      ggtitle(paste("PCoA ordination: ", name, "\nMethod: '", dmethod, "'", sep="")) +
      xlab(paste("Axis 1 (", round(p$values[1,3]*100,digits = 1), "%)", sep = "")) +
      ylab(paste("Axis 2 (", round(p$values[2,3]*100,digits = 1), "%)", sep = ""))
      
    # Change plotting type if numeric vs. non-numeric
    if (is.numeric(metadata) == TRUE) {
      gradient_palette <- rev(RColorBrewer::brewer.pal(name = "YlGnBu", n = 5))
      
      z_final <- z + 
        geom_point(aes(fill = Metadata), shape = 21, size = 3.5, alpha = 0.8) +
        scale_fill_gradientn(colours = gradient_palette, na.value = "grey50") +
        # Can change legend title from NULL to name if desired
        guides(fill = guide_colourbar(title=NULL))
      
    } else {
      metadata <- factor(metadata)
      
      z_final <- z + 
        geom_point(aes(fill = metadata), shape = 21, size = 3.5, alpha = 0.8) +
        scale_fill_manual(values = metadata_colours, na.value="white") +
        # Can change legend title from NULL to name if desired
        guides(fill = guide_legend(title=NULL))
    }
    
    if (!is.null(ellipsoidConf) && ellipsoidConf != 0) {
      ellipse <- dataEllipse(p$vectors[,1],p$vectors[,2],groups=metadata,levels=as.numeric(ellipsoidConf),add=TRUE,plot.points=FALSE,grid=FALSE,center.pch=FALSE,col=rainbow(fac.len))
    }
    
    return(z_final)
    
  } else {
    stop("Something is fishy with one of the metadata categories. Having a hard time determining the number of factors. Exiting...")
  }

}

# Make plots
plot_list <- lapply(1:ncol(mapping), make_plot)

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
