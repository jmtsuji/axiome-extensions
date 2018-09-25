# Functions for parsing silva 7-rank taxonomy in an OTU table
# Copyright Jackson M. Tsuji, 2018
# Load these into an R script to process taxonomy.

# Load libraries
library(glue)
library(parallel)

# Function: Recursively identify ambiguity hits (i.e., occurances of TRUE) starting from the back of the taxonomy_ranks
# Inputs: 'ambiguous_ranks' - logical (length 7) of whether or not each position in a taxonomy_ranks vector is considered ambiguous
        # 'rank_number' - numeric (length 1). Must start at the end of the ambiguous_ranks, and then the function will propogate back to the start
# Return: 'ambiguous_rank_position' - numeric (length 1) of the 'last' ambiguous position in the ranks when moving back from the end of the ranks
check_ambiguity_pattern <- function(ambiguous_ranks, rank_number) {
  if ( ambiguous_ranks[rank_number] == TRUE ) {
    if ( rank_number == 1 ) {
      # End point 3: this means that you've reached all the way to the front of the taxonomy ranks and all have been ambiguous.
      # Return "1" to the user
      ambiguous_rank_position <- rank_number
      return(ambiguous_rank_position)
      
    } else if ( rank_number > 1 && rank_number <= 7 ) {
      # Recursion: If you're still at the far right or the middle of the taxonomic ranks, then all is well. Keep moving back until you hit the end.
      message("Found ambiguity at position ", rank_number, ". Moving back one.")
      check_ambiguity_pattern(ambiguous_ranks, (rank_number - 1))
      
    } else {
      stop("rank_number must be a number <= 7. You provided '", rank_number, "'. Exiting...")
    }
    
  } else if ( ambiguous_ranks[rank_number] == FALSE ) {
    
    if ( rank_number == length(ambiguous_ranks) ) {
      # End point 1: this means you've just started running from the end of the ranks and found no ambiguities on the far right side
      # Thus, no ambiguous ranks of the proper pattern are here.
      ambiguous_rank_position <- NA
      return(ambiguous_rank_position)
      
    } else if ( rank_number < length(ambiguous_ranks) ) {
      # End point 2: you've propagated up the taxonomy ranks and eventually hit a rank that is non-ambiguous. Flag this.
      ambiguous_rank_position <- rank_number + 1
      return(ambiguous_rank_position)
      
    } else {
      stop("rank_number must be a number <= 7. You provided '", rank_number, "'. Exiting...")
    }
    
  }
}

# Function: given an unparsed input string of silva taxonomy, outputs the ranks in a character vector. Adds meaning to "ambiguous" ones.
# Input: a character vector (length 1) of unparsed silva taxonomy. MUST be 7-rank taxonomy.
# Return: a character vector (length 7) of parsed silva taxonomy. 'D_' header is removed. Ambiguous features have proceeding meaningful taxonomic classifications appended.
parse_silva_taxonomy <- function(taxonomy_string) {
  # Divide into ranks based on the ';' separator
  taxonomy_ranks <- strsplit(x = taxonomy_string, split = "; ")[[1]]
  
  # Confirm output is 7 long
  if ( length(taxonomy_ranks) != 7 ) {
    stop("Entry should have 7 taxonomic ranks but has ", length(taxonomy_ranks), ". Here are the ranks: '",
         glue::glue_collapse(taxonomy_ranks, ", "), "'. Exiting...")
  }
  
  # Remove leading 'D_' (or just first five characters, apparently)
  # TODO - look up pattern
  silva_leading_pattern <- "" # TODO
  taxonomy_ranks <- gsub(pattern = silva_leading_pattern, replacement = "", x = taxonomy_ranks)
  
  # HARD-CODED ambigous taxon names in silva
  ambiguous_terms <- c("uncultured bacterium", "uncultured", "Ambiguous_taxa")
  
  # Find all occurances of the ambiguous_terms among the taxonomy_ranks (output is TRUE/FALSE for each taxonomy rank)
  ambiguous_ranks <- taxonomy_ranks %in% ambiguous_terms
  
  # Get the position in the taxonomic ranks where the chain of ambiguous names ends (starting from the back)
  ambiguous_rank_position <- check_ambiguity_pattern(ambiguous_ranks, rank_number = 7)
  
  # Check for orphaned ambiguous ranks (i.e., ones that don't propagate in a chain from the back)
  # Don't modify orphaned entries, but throw a warning
  first_ambiguous_position <- match(TRUE, ambiguous_ranks)
  if ( is.na(first_ambiguous_position) == TRUE ) {
    # Don't do anything -- there are no ambiguous ranks for this entry.
  } else if ( first_ambiguous_position < ambiguous_rank_position ) {
    
    warning("Ambiguous taxonomy terms don't occur in a natural chain moving up from 'species' for this taxonomy entry: '",
         glue::glue_collapse(taxonomy_ranks, ", "), "'. Will only correct for the entries in a chain at the end.")
    
  }
  
  # Append taxonomic information to ambiguous entries
  if ( is.na(ambiguous_rank_position) == TRUE ) {
    # All done -- no ambiguous entries
  } else if ( ambiguous_rank_position > 1 && ambiguous_rank_position <=7 ) {
    
    # Get the taxon rank immediately above the ambiguous_rank_position
    last_informative_taxon <- taxonomy_ranks[ambiguous_rank_position - 1]
    
    # Append to end of ambiguous ranks
    for ( i in ambiguous_rank_position:7 ) {
      taxonomy_ranks[i] <- paste(taxonomy_ranks[i], "_", last_informative_taxon, sep = "")
    }
    
  }
  
  return(taxonomy_ranks)
  
}

# Function: adds the parsed silva taxonomy onto an OTU table loaded as a data frame
# Inputs: 'OTU_table' - data frame of the OTU table, with 'sites' in the columns and OTUs in the rows
        # 'threads' - numeric (length 1) of the number of processor threads to run this function with
# Output: 'OTU_table' - data frame with added parsed taxonomy (if Consensus.Lineage was present)
add_silva_taxonomy_to_OTU_table <- function(OTU_table, threads) {
  
  # # Uncomment this if you wanted to read in the OTU table from a file.
  # OTU_table <- read.table(OTU_table_filename, header = TRUE, sep="\t", comment.char = "", skip = 1, stringsAsFactors = FALSE)
  
  # Look for consensus lineage column
  if ( "Consensus.Lineage" %in% colnames(OTU_table) ) {
    
    # TODO - add a check that the taxonomy is silva and not greengenes
    
    # Parse the taxonomy (multi-threaded support)
    parsed_taxonomy <- mclapply(OTU_table$Consensus.Lineage, parse_silva_taxonomy, mc.cores = threads)
    
    # Bind into a data frame
    parsed_taxonomy <- as.data.frame(matrix(unlist(parsed_taxonomy), ncol = 7, byrow = TRUE))
    colnames(parsed_taxonomy) <- c("t_Domain", "t_Phylum", "t_Class", "t_Order", 
                                   "t_Family", "t_Genus", "t_Species")
    
    # Attach to the OTU table
    OTU_table <- cbind(OTU_table, parsed_taxonomy)
  } else {
    warning("The 'Consensus.Lineage' column is not in the provided OTU table, so can't parse taxonomy. Returning the original OTU table.")
  }
  
  return(OTU_table)
  
}

