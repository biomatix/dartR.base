#' @name gl.smearplot
#' @title A plot of individuals against loci, showing the character state [0, 1, 2, NA]
#' @family graphics

#' @description
#' Each locus is color coded for scores of 0, 1, 2 and NA for SNP data and 0, 1
#' and NA for presence/absence (SilicoDArT) data. Individual labels can be added.

#' Plots can become cluttered when there are many individuals. Use
#' \code{ind.labels = FALSE} to fall back to numeric ticks, or
#' \code{ind.labels = "auto"} to let the function decide based on
#' \code{label.limit}.
#'
#' Individuals can be reordered via \code{ind.order}; note that any
#' accompanying dendrogram should not be drawn in that case, as its tip order
#' will no longer correspond to the smear plot y-axis. Loci can be reordered
#' via \code{loc.order}, including by a named locus metric such as
#' \code{"AvgPIC"} or \code{"CallRate"}.
#'
#' Works with both SNP data and P/A data (SilicoDArT).

#' @param x Name of the genlight object [required].
#' @param ind.labels One of TRUE (show individual IDs), FALSE (numeric ticks),
#' or \code{"auto"} (show IDs when \code{nInd(x) <= label.limit}, otherwise
#' numeric ticks) [default FALSE].
#' @param label.size Size of the individual labels [default 10].
#' @param label.limit Maximum \code{nInd(x)} for which \code{ind.labels = "auto"}
#' will display individual IDs [default 50].
#' @param ind.order Controls y-axis ordering. One of: \code{NULL} (default
#' order, matching prior behaviour); \code{"pop"} (order individuals by
#' population without faceting); or a character vector that is a permutation
#' of \code{indNames(x)}. When an explicit order is supplied, any accompanying
#' dendrogram should be suppressed because its tip order will no longer match.
#' [default NULL].
#' @param loc.order Controls x-axis ordering. One of: \code{NULL} (retain the
#' locus order in \code{x}); the name of a column in
#' \code{x@other$loc.metrics} (e.g. \code{"AvgPIC"}, \code{"CallRate"}); a
#' character vector that is a permutation of \code{locNames(x)}; or an integer
#' permutation of \code{1:nLoc(x)}. [default NULL].
#' @param loc.order.decreasing When \code{loc.order} names a locus metric, sort
#' in decreasing order if TRUE [default TRUE].
#' @param group.pop If TRUE, facet the plot by population [default FALSE].
#' @param het.only If TRUE, show only the heterozygous state (SNP data only)
#' [default FALSE].
#' @param plot.display If TRUE, the plot is displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default NULL].
#' @param plot.colors List of four color names for the column fill for homozygous reference,
#' heterozygous, homozygous alternate, and missing value (NA) [default c("#0000FF","#00FFFF","#FF0000","#e0e0e0")].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param legend Position of the legend: "left", "top", "right", "bottom" or
#'  "none" [default = "bottom"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#'
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' gl.smearplot(testset.gl,ind.labels=FALSE)
#' gl.smearplot(testset.gs,ind.labels=FALSE)
#' gl.smearplot(testset.gl[1:10,],ind.labels=TRUE)
#' gl.smearplot(testset.gs[1:10,],ind.labels=TRUE)
#' # Order loci by average PIC (SNP datasets with AvgPIC in loc.metrics):
#' # gl.smearplot(testset.gl, loc.order = "AvgPIC")
#' # Custom individual order (dendrogram cannot accompany):
#' # gl.smearplot(testset.gl[1:10,], ind.labels=TRUE,
#' #              ind.order = rev(indNames(testset.gl[1:10,])))
#' 
# TEST SCRIPTS
# gl.smearplot(testset.gl)                                  # default — should match old output
# gl.smearplot(testset.gl, ind.labels = TRUE)               # alphabetical names — should match old
# gl.smearplot(testset.gl, ind.labels = "auto")             # decides based on nInd
# gl.smearplot(testset.gl, loc.order = "AvgPIC")            # loci ordered by AvgPIC desc
# gl.smearplot(testset.gl, ind.order = rev(indNames(testset.gl)))
# gl.smearplot(testset.gs)                                   # SilicoDArT legend should now say Absence/Presence
# gl.smearplot(testset.gs, het.only = TRUE)                 # should warn and still render SilicoDArT correctly

#' @export
#' @return Returns the ggplot object
#'
# TEST
# ddd <- matrix(data=0,nrow=10,ncol=10)
# ddd[8,10] <- NA
# ddd[9,10] <- 2
# ddd[10,10] <- 2
# ddd
# ddd <- as.genlight(ddd)
# ploidy(ddd) <- 2
# ddd <- gl.compliance.check(ddd)
# gl.smearplot(ddd)

gl.smearplot <- function(x,
                        plot.display = TRUE,
                        ind.labels = FALSE,
                        label.size = 10,
                        label.limit = 50,
                        ind.order = NULL,
                        loc.order = NULL,
                        loc.order.decreasing = TRUE,
                        group.pop = FALSE,
                        plot.theme = NULL,
                        plot.colors = NULL,
                        plot.file = NULL,
                        plot.dir = NULL,
                        het.only = FALSE,
                        legend = "bottom",
                        verbose = NULL) {

    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "reshape2"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if (verbose == 0) { plot.display <- FALSE }

    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir, verbose = 0)

    # SET COLOURS
    if (is.null(plot.colors)) {
      plot.colors <- c("#0000FF","#00FFFF","#FF0000","#e0e0e0")
    } else {
      if (length(plot.colors) > 4) {
        if (verbose >= 2) cat(warn("  Specified plot colours exceed 4, first 4 only are used\n"))
        plot.colors <- plot.colors[1:4]
      }
    }

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # Apply het.only (SNP only). Moved above the colour-count munging so that
    # for SilicoDArT we warn and leave plot.colors untouched (previous code
    # recoloured then tried to recover, corrupting the SilicoDArT palette).
    if (het.only) {
      if (datatype == "SilicoDArT") {
        if (verbose >= 2) cat(warn("  het.only applies to SNP data only; ignoring\n"))
        het.only <- FALSE
      } else {
        plot.colors <- c("#d3d3d3","#00FFFF","#d3d3d3","#e0e0e0")
      }
    }

    # RESOLVE ind.labels -> show_names (accepts TRUE / FALSE / "auto")
    if (is.character(ind.labels) && length(ind.labels) == 1 &&
        ind.labels == "auto") {
      show_names <- nInd(x) <= label.limit
      if (verbose >= 2) {
        cat(report(sprintf(
          "  ind.labels='auto': nInd=%d, label.limit=%d -> %s individual names\n",
          nInd(x), label.limit, if (show_names) "showing" else "hiding")))
      }
    } else if (isTRUE(ind.labels)) {
      show_names <- TRUE
    } else if (isFALSE(ind.labels)) {
      show_names <- FALSE
    } else {
      stop(error("Fatal error: ind.labels must be TRUE, FALSE, or 'auto'\n"))
    }

    # RESOLVE ind.order. Defaults below match the original as.factor() behaviour:
    #   show_names=TRUE  -> alphabetical (as.factor(character) sorts levels)
    #   show_names=FALSE -> positional   (as.factor(1:N) sorts numerically 1..N)
    if (is.null(ind.order)) {
      ind_order <- if (show_names) sort(indNames(x)) else indNames(x)
    } else if (is.character(ind.order) && length(ind.order) == 1 &&
               ind.order == "pop") {
      ind_order <- indNames(x)[order(as.character(pop(x)))]
    } else if (is.character(ind.order)) {
      if (!setequal(ind.order, indNames(x))) {
        miss <- setdiff(indNames(x), ind.order)
        ext  <- setdiff(ind.order, indNames(x))
        stop(error(sprintf(
          "Fatal error: ind.order must be a permutation of indNames(x). Missing: %s; Unknown: %s\n",
          paste(head(miss, 5), collapse = ", "),
          paste(head(ext,  5), collapse = ", "))))
      }
      ind_order <- ind.order
    } else {
      stop(error("Fatal error: ind.order must be NULL, 'pop', or a character vector of individual names\n"))
    }

    # RESOLVE loc.order. Stored as an integer permutation of 1:nLoc(x) that
    # the factor-level reordering below uses verbatim.
    loc_axis_label <- "Loci"
    if (is.null(loc.order)) {
      loc_perm <- seq_len(nLoc(x))
    } else if (is.character(loc.order) && length(loc.order) == 1 &&
               !is.null(x@other$loc.metrics) &&
               loc.order %in% colnames(x@other$loc.metrics)) {
      metric_values <- x@other$loc.metrics[[loc.order]]
      loc_perm <- order(metric_values,
                        decreasing = loc.order.decreasing,
                        na.last = TRUE)
      loc_axis_label <- sprintf("Loci (ordered by %s, %s)",
                                loc.order,
                                if (loc.order.decreasing) "decreasing" else "increasing")
    } else if (is.character(loc.order)) {
      if (!setequal(loc.order, locNames(x))) {
        stop(error("Fatal error: loc.order as character must be a permutation of locNames(x) or a single locus-metric column name\n"))
      }
      loc_perm <- match(loc.order, locNames(x))
    } else if (is.numeric(loc.order)) {
      if (!setequal(loc.order, seq_len(nLoc(x)))) {
        stop(error("Fatal error: loc.order as numeric must be a permutation of 1:nLoc(x)\n"))
      }
      loc_perm <- as.integer(loc.order)
    } else {
      stop(error("Fatal error: loc.order must be NULL, a locus metric name, a locus-name character vector, or an integer permutation\n"))
    }

    # BUILD DATAFRAME FROM GENLIGHT
    df.matrix <- as.data.frame(as.matrix(x))
    colnames(df.matrix) <- as.character(seq_len(nLoc(x)))
    df.matrix$id  <- indNames(x)
    df.matrix$pop <- pop(x)

    df.listing <- reshape2::melt(df.matrix, id.vars = c("pop", "id"))
    df.listing$value <- as.character(df.listing$value)
    df.listing$value <- ifelse(df.listing$value == "NA", NA, df.listing$value)
    colnames(df.listing) <- c("pop", "id", "locus", "genotype")

    # Apply orderings via factor levels. With the defaults above, this
    # reproduces the previous visual output exactly.
    df.listing$id    <- factor(df.listing$id,    levels = ind_order)
    df.listing$locus <- factor(df.listing$locus, levels = as.character(loc_perm))

    # Tick locations for discrete axes
    loc_labels <- pretty(1:nLoc(x), 5)
    loc_labels <- loc_labels[loc_labels >= 1 & loc_labels <= nLoc(x)]
    id_labels  <- pretty(1:nInd(x), 5)
    id_labels  <- id_labels[id_labels  >= 1 & id_labels  <= nInd(x)]

    locus <- id <- genotype <- NA  # silence R CMD check for aes() refs

    # Assign colours and labels for genotypic data
    labels_genotype <- as.character(unique(df.listing$genotype))
    labels_genotype <- labels_genotype[!is.na(labels_genotype)]
    labels_genotype <- labels_genotype[order(labels_genotype)]
    plot.colors.hold <- plot.colors
    tmp <- NULL
    if (length(labels_genotype) < 3) {
      if ("0" %in% labels_genotype) { tmp[1] <- plot.colors[1] }
      if ("1" %in% labels_genotype) {
        tmp <- if (is.null(tmp)) plot.colors[2] else c(tmp, plot.colors[2])
      }
      if ("2" %in% labels_genotype) {
        tmp <- if (is.null(tmp)) plot.colors[3] else c(tmp, plot.colors[3])
      }
      tmp <- c(tmp, plot.colors[4])
      plot.colors <- tmp
    }
    n.colors <- length(plot.colors)

    labels_genotype[which(is.na(labels_genotype))] <- "Missing data"
    labels_genotype[labels_genotype == "0"] <- "Homozygote reference"
    labels_genotype[labels_genotype == "1"] <- "Heterozygote"
    labels_genotype[labels_genotype == "2"] <- "Homozygote alternate"

    # SilicoDArT legend labels. Positional: first slot is "0" -> Absence,
    # second is "1" -> Presence. (Previous code assigned via character
    # indices on an unnamed vector, which appended rather than replaced.)
    labels_silicodart <- c("Absence", "Presence")

    # BUILD PLOT: single code path. Datatype selects the fill scale; show_names
    # selects between named and numeric-tick y-axes.
    p3 <- ggplot(df.listing, aes(x = locus, y = id, fill = genotype)) +
      geom_raster()

    if (datatype == "SNP") {
      p3 <- p3 + scale_fill_discrete(
        type     = plot.colors,
        na.value = plot.colors[n.colors],
        name     = "Genotype",
        labels   = labels_genotype
      )
    } else if (datatype == "SilicoDArT") {
      # The <3-genotype-class munging above is for SNP only; restore the
      # unmunged palette before indexing into it.
      plot.colors <- plot.colors.hold
      p3 <- p3 + scale_fill_discrete(
        type     = plot.colors[c(1, 3)],
        na.value = plot.colors[4],
        name     = "Sequence Tag",
        labels   = labels_silicodart
      )
    }

    # X-axis: discrete with positional ticks. breaks map through loc_perm so
    # tick "200" sits at display position 200 in the (possibly reordered) axis.
    p3 <- p3 + scale_x_discrete(
      breaks   = as.character(loc_perm[loc_labels]),
      labels   = as.character(loc_labels),
      name     = loc_axis_label,
      position = "bottom"
    )

    # Y-axis: either label every individual or fall back to numeric ticks.
    if (show_names) {
      p3 <- p3 + scale_y_discrete(name = "Individuals")
    } else {
      p3 <- p3 + scale_y_discrete(
        breaks   = ind_order[id_labels],
        labels   = as.character(id_labels),
        name     = "Individuals",
        position = "left"
      )
    }

    if (!is.null(plot.theme)) {
      p3 <- p3 + plot.theme
    }

    p3 <- p3 + theme(
      legend.position = legend,
      axis.text.y     = element_text(size = label.size)
    )

    if (group.pop) {
      p3 <- p3 + facet_wrap(~ pop,
                            ncol   = 1,
                            dir    = "v",
                            scales = "free_y")
    }

    # PRINTING OUTPUTS
    print(p3)

    # Optionally save the plot
    if (!is.null(plot.file)) {
      tmp <- utils.plot.save(p3,
                             dir     = plot.dir,
                             file    = plot.file,
                             verbose = verbose)
    }

    # FLAG SCRIPT END
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }

    # RETURN
    invisible(p3)
}
