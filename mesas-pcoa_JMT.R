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
setwd("/Users/JTsuji/Documents/Research_General/PhD/04e_Biogeography/06_iTag_analysis/02_Sep2017_ELA_CyaMig_AXIOME/CyaMig_01/convert_otu_table/rarefied_biom_to_tab0/03_PCoA_noCTRL/")
outName <- "pcoa"
otuTable <- "otu_table_auto_noCtrl.tab"
mappingFile <- "171115_metadata_CyaMig_MOD_noCtrl.tsv"

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

pkgTest('getopt')

#Make sure we have a valid distance method
dlist = c('manhattan','euclidean','canberra','bray','kulczynski','jaccard','gower','altGower','morisita','horn','mountford','raup','binomial','chao','cao')

if (! dmethod %in% dlist ) {
	print(paste("Invalid distance method:", dmethod))
	q(status=1)
} 

pkgTest("vegan")
pkgTest("ape")
if (!is.null(ellipsoidConf) ) {
	pkgTest("car")
}
#Read in and format otu table
print("Reading OTU table")

rawtable <- read.table(otuTable, skip = 1,
comment.char = "", header = TRUE, row.names = 1, sep = "\t")
otutable.first <- t(rawtable[1:(ncol(rawtable) - 1)])
#Sort numerically the samples
#otutable <- otutable[order(as.integer(sub("X","", rownames(otutable)))),]

print("Reading mapping");
mapping.first <- read.table(mappingFile, header = TRUE, comment.char = "", row.names = 1, sep = "\t") 
                            #colClasses = "character")
# For non-numeric sample names, comment out the following line
#colnames(mapping) <- paste("X", colnames(mapping), sep = "");
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

print("Computing PCoA");

d <- vegdist(otutable, method=dmethod);

print("Making PCoA Plot");
#pdf(outName,useDingbats=F);
p <- pcoa(d);

# Uncomment to reverse the orientations of Axis 1 or Axis 2 points, respectively (just for visualization purposes; does not affect the data)
# p$vectors[,1]=-p$vectors[,1]
# p$vectors[,2]=-p$vectors[,2]

# Preparing to sort by metadata category
p.sort <- p$vectors
p.sort <- cbind("sort"=c(1:nrow(p.sort)),p.sort)
mapping.sort = cbind("sort"=c(1:nrow(mapping)),mapping)
if (identical(row.names(p.sort), row.names(mapping.sort)) == FALSE) {
  print("ERROR: Row names in sorting vectors do not match. Results will be unreliable.")
}
# Consider impoving this sanity check (above) to exit out of the script if it fails.

# Print out vectors for user
p.print <- data.frame(Axis.1 = p$vectors[,1], Axis.2 = p$vectors[,2])
write.table(p.print, file = "plotting_values.tsv", sep = "\t", col.names = NA)

pkgTest('ggplot2')

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
  
  name <- colnames(mapping)[mapping_col]
  metadata <- factor(mapping[rownames(mapping),mapping_col])
  fac.len <- length(levels(metadata))
  if (fac.len < 2 ) {
    print(paste("Ignoring", name, "because all values are identical."))
  } else if (fac.len >= 2) {
    metadata_colours <- metadata
    # Can adjust start and end below to change spectrum of colours plotted
    levels(metadata_colours) <- rainbow(fac.len,start=1/7,end=3/4);
    #metadata_colours<-unlist(lapply(metadata_colours, function(mapping_col) { substr(mapping_col, 0, 7); }));
    p.frame = data.frame(Axis.1 = p$vectors[,1], Axis.2 = p$vectors[,2])
    z = ggplot(p.frame) +
      geom_point(aes(x = Axis.1, y = Axis.2, colour = factor(metadata))) +
      theme_bw() +
      ggtitle(paste("PCoA ordination: ", name, "\nMethod: '", dmethod, "'", sep="")) +
      xlab(paste("Axis 1 (", round(p$values[1,3]*100,digits = 1), "%)", sep = "")) +
      ylab(paste("Axis 2 (", round(p$values[2,3]*100,digits = 1), "%)", sep = "")) +
      # Can change legend title from NULL to name if desired
      guides(colour = guide_legend(title=NULL)) +
      scale_colour_manual(values=levels(metadata_colours), na.value="white")
    # print(z)
    
    if (!is.null(ellipsoidConf) && ellipsoidConf != 0) {
      ellipse <- dataEllipse(p$vectors[,1],p$vectors[,2],groups=metadata,levels=as.numeric(ellipsoidConf),add=TRUE,plot.points=FALSE,grid=FALSE,center.pch=FALSE,col=rainbow(fac.len))
    }
    
    return(z)
    
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

print("Making biplot");
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
