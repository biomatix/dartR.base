#' @name gl.read.csv 
#' 
#' @title Reads SNP data from a csv file into a genlight object 
#' @family io 
#' @description 
#' This script takes SNP genotypes from a csv file, combines them with 
#' individual and locus metrics and creates a genlight object. 
#' The SNP data need to be in one of two forms. SNPs can be coded 0 for 
#' homozygous reference, 2 for homozygous alternate, 1 for heterozygous, and NA 
#' for missing values; or the SNP data can be coded A/A, A/C, C/T, G/A etc, 
#' and -/- for missing data. In this format, the reference allele is the most 
#' frequent allele, as used by DArT. Other formats will throw an error. 
#' The SNP data need to be individuals as rows, labeled, and loci as columns, 
#' also labeled. If the orientation is individuals as columns and loci by rows, 
#' then set transpose=TRUE. 
#' The individual metrics need to be in a csv file, with headings, with a 
#' mandatory id column corresponding exactly to the individual identity labels 
#' provided with the SNP data and in the same order. 
#' The locus metadata needs to be in a csv file with headings, with a mandatory 
#' column headed AlleleID corresponding exactly to the locus identity labels 
#' provided with the SNP data and in the same order. 
#' Note that the locus metadata will be complemented by calculable statistics 
#' corresponding to those that would be provided by Diversity Arrays Technology 
#' (e.g. CallRate). 
#' 
#' @param filename Name of the csv file containing the SNP genotypes [required]. 
#' @param transpose If TRUE, rows are loci and columns are individuals 
#' [default FALSE]. 
#' @param ind.metafile Name of the csv file containing the metrics for 
#' individuals [optional]. 
#' @param loc.metafile Name of the csv file containing the metrics for 
#' loci [optional]. 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log; 3, progress and results summary; 5, full report 
#' [default 2 or as specified using gl.set.verbosity]. 
#' 
#' @author Custodian: Luis Mijangos -- Post to 
#' \url{https://groups.google.com/d/forum/dartr} 
#' # @examples 
#' #Taken out, will not build 
#' # csv_file <- system.file('extdata','platy_test.csv', package='dartR.data') 
#' # ind_metadata <- system.file('extdata','platy_ind.csv', package='dartR.data') 
#' # gl <- gl.read.csv(filename = csv_file, ind.metafile = ind_metadata) 
#' @export 
#' @return A genlight object with the SNP data and associated metadata included.
#' 
gl.read.csv <- function(filename,
                        transpose = FALSE,
                        ind.metafile = NULL,
                        loc.metafile = NULL,
                        verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (is.null(loc.metafile) & verbose > 0) {
    cat(
      warn(
        "Warning: Locus metafile not provided, locus metrics will be calculated where this is possible\n"
      )
    )
  }
  
  if (is.null(ind.metafile) & verbose > 0) {
    cat(
      warn(
        "Warning: Individual metafile not provided, pop set to 'A' for all individuals\n"
      )
    )
  }
  
  # DO THE JOB
  # FIRST THE SNP DATA
  
  # Read genotype CSV – avoid factors so that numeric conversion is correct
  df0 <-
    read.csv(file = filename,
             header = FALSE,
             stringsAsFactors = FALSE)
  
  if (transpose) {
    df0 <- t(df0)
  }
  
  numrows <- dim(df0)[1]
  numcols <- dim(df0)[2]
  
  if (verbose > 0) {
    cat(
      report(
        "Input data should be a csv file with individuals as rows, loci as columns\n"
      )
    )
    
    showL <- min(5, numcols - 1)
    showI <- min(5, numrows - 1)
    
    cat("  ", numcols - 1, "loci, confirming first:",
        as.matrix(df0[1, 2:(1 + showL)]), "\n")
    cat("  ", numrows - 1, "individuals, confirming first:",
        as.matrix(df0[2:(1 + showI), 1]), "\n")
    
    cat(important("    If these are reversed, re-run the script with transpose=TRUE\n"))
  }
  
  # Extract SNP matrix
  data <- as.matrix(df0[2:numrows, 2:numcols])
  
  loci <- df0[1, 2:numcols]
  loci <- as.character(as.matrix(loci))
  
  individuals <- df0[2:numrows, 1]
  individuals <- as.character(individuals)
  
  # Uniqueness checks
  if (length(unique(individuals)) != length(individuals)) {
    cat(error("Fatal Error: Individual labels are not unique, check and edit your input file\n"))
    stop()
  }
  if (length(unique(loci)) != length(loci)) {
    cat(error("Fatal Error: AlleleID not unique, check and edit your input file\n"))
    stop()
  }
  
  # -----------------------------------------------------------
  # VALIDATE AND CONVERT SNP DATA
  # -----------------------------------------------------------
  
  test <- paste0(data, collapse = "")
  test <- gsub("NA", "9", test)
  test <- gsub(" ", "", test)
  
  # -----------------------------------------------------------------
  # CHARACTER-CODED GENOTYPES (A/T, G/C, A/A, -/-)
  # -----------------------------------------------------------------
  if (nchar(test) > nrow(data) * ncol(data)) {
    
    if (verbose >= 2) {
      cat(report("Character data detected, assume genotypes are of the form C/C, A/T, C/G, -/- etc\n"))
    }
    
    # Global check: only A/C/G/T/- are allowed
    s1 <- paste(data, collapse = " ")
    s1 <- gsub("/", " ", s1)
    s1 <- toupper(s1)
    s2 <- unlist(strsplit(s1, " "))
    tmp <- table(s2)
    
    if (!all(names(tmp) %in% c("A", "C", "G", "T", "-"))) {
      cat(error("Fatal Error: Genotypes must be defined by letters A, C, G, T or missing -\n"))
      stop()
    }
    
    # For each locus: confirm biallelic and convert to 0/1/2
    for (i in 1:dim(data)[2]) {
      v1 <- data[, i]
      v1 <- paste(v1, collapse = " ")
      v1 <- gsub("/", " ", v1)
      v1 <- toupper(v1)
      v1 <- unlist(strsplit(v1, " "))
      
      # Remove missing codes and empty strings
      v1 <- v1[v1 != "" & v1 != "-"]
      
      tmp <- table(v1)
      tmp <- tmp[order(as.numeric(tmp), decreasing = TRUE)]
      
      if (length(names(tmp)) > 2) {
        cat(error("Fatal Error: Loci are not bi-allelic\n"))
        stop()
      }
      
      homRef <- paste0(names(tmp)[1], "/", names(tmp)[1])
      homAlt <- paste0(names(tmp)[2], "/", names(tmp)[2])
      het1   <- paste0(names(tmp)[1], "/", names(tmp)[2])
      het2   <- paste0(names(tmp)[2], "/", names(tmp)[1])
      
      data[, i] <- gsub(homRef, "0", data[, i])
      data[, i] <- gsub(homAlt, "2", data[, i])
      data[, i] <- gsub(het1,   "1", data[, i])
      data[, i] <- gsub(het2,   "1", data[, i])
      data[, i] <- gsub("-/-",  NA, data[, i])
    }
    
    if (verbose >= 2) {
      cat(report("  Data confirmed as biallelic\n"))
      cat(report("  SNP coding converted to 0, 1, 2 and NA\n"))
    }
    
    data <- apply(data, 2, as.numeric)
    
  } else {
    # -----------------------------------------------------------------
    # NUMERIC-CODED GENOTYPES (0/1/2/NA)
    # -----------------------------------------------------------------
    
    if (verbose >= 2) {
      cat(report("  Numeric data detected, assume genotypes are 0, 1, 2, NA\n"))
    }
    
    # Convert safely
    data <- apply(data, 2, function(x) as.numeric(x))
    
    # Validate codes
    s1 <- paste(data, collapse = " ")
    s2 <- unlist(strsplit(s1, " "))
    tmp <- table(s2)
    
    valid_codes <- c("0", "1", "2", "NA")
    
    if (!all(names(tmp) %in% valid_codes)) {
      cat(error("Fatal Error: Genotypes must be defined by 0, 1, 2 or NA\n"))
      stop()
    }
  }
  
  # -----------------------------------------------------------------
  # CREATE GENLIGHT OBJECT
  # -----------------------------------------------------------------
  
  gl <-
    new("genlight",
        data,
        ploidy = 2,
        loc.names = loci,
        ind.names = individuals)
  
  pop(gl) <- array("A", nInd(gl))
  gl <- gl.compliance.check(gl, verbose = verbose)
  
  # Correct initialisation of individual metrics
  gl@other$ind.metrics <-
    data.frame(id = indNames(gl),
               pop = array("A", nInd(gl)),
               stringsAsFactors = TRUE)
  
  # -----------------------------------------------------------------
  # LOCUS METADATA
  # -----------------------------------------------------------------
  
  if (!is.null(loc.metafile)) {
    loc.metrics <-
      read.csv(file = loc.metafile,
               header = TRUE,
               stringsAsFactors = TRUE)
    
    if (!("AlleleID" %in% names(loc.metrics))) {
      cat(error("Fatal Error: mandatory AlleleID column absent from locus metrics file\n"))
      stop()
    }
    
    for (i in 1:nLoc(gl)) {
      if (loc.metrics[i, 1] != gl@other$loc.metrics$AlleleID[i]) {
        stop(error(
          "Fatal Error: AlleleID in locus metrics file does not match input data, or is out of order\n"
        ))
      }
    }
    
    gl@other$loc.metrics <- loc.metrics
  }
  
  gl <- gl.recalc.metrics(gl, verbose = 0)
  
  if (verbose >= 2) {
    cat(report("  Added or updated locus metrics in other$loc.metrics\n"))
  }
  
  # -----------------------------------------------------------------
  # INDIVIDUAL METADATA
  # -----------------------------------------------------------------
  
  if (!is.null(ind.metafile)) {
    ind.metrics <-
      read.csv(file = ind.metafile,
               header = TRUE,
               stringsAsFactors = TRUE,
               fileEncoding = "UTF-8-BOM")
    
    if (!("id" %in% names(ind.metrics))) {
      cat(error("Fatal Error: mandatory id column absent from individual metadata file\n"))
      stop()
    }
    
    for (i in 1:nInd(gl)) {
      if (ind.metrics[i, 1] != gl@other$ind.metrics$id[i]) {
        cat(error(
          "Fatal Error: id in individual metrics file does not correspond to SNP input file, or is out of order\n"
        ))
        stop()
      }
    }
    
    if (!("pop" %in% names(ind.metrics))) {
      cat(warn("  Warning: pop column absent in individual metadata; setting pop='A'\n"))
      gl@other$ind.metrics <- ind.metrics
      gl@other$ind.metrics$id <- individuals
      gl@other$ind.metrics$pop <- array("A", nInd(gl))
      pop(gl) <- gl@other$ind.metrics$pop
    } else {
      gl@other$ind.metrics <- ind.metrics
      gl@other$ind.metrics$id <- individuals
      gl@other$ind.metrics$pop <- ind.metrics$pop
      pop(gl) <- gl@other$ind.metrics$pop
    }
    
    if (verbose >= 2) {
      cat(report("  Added individual metadata to other$ind.metrics\n"))
    }
  }
  
  # MAKE COMPLIANT
  gl <- gl.compliance.check(gl, verbose = verbose)
  
  # ADD HISTORY
  gl@other$history <- list()
  gl@other$history[[1]] <- match.call()
  
  # END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(gl)
}
