suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(logging)))
suppressMessages(suppressWarnings(library(scater)))

filterCells <- function(kallisto.dirs, filtered.json) {
    all.dirs <- list.dirs(path = kallisto.dirs, full.names = TRUE, recursive = FALSE)
    loginfo(sprintf("Reading kallisto subdirs from %s", kallisto.dirs))

    sce <- readKallistoResults(directories = all.dirs, samples = all.dirs, read_h5=TRUE) 
    #sce <- readKallistoResultsOneSample(kallisto.dir, read_h5 = FALSE, kallisto_version = "current")
    sce <- calculateQCMetrics(sce)
    print(sce)
    colnames(colData(sce))
    plotQC(sce, type = "highest-expression", exprs_values = "counts")
    pdf("qc.pdf")
        plotQC(sce, type = "highest-expression", exprs_values = "counts")
    dev.off()
    return(sce)
    #return(rownames(sce))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    warning("qc.R: Not exact number of arguments. Expected 2 args: kallisto directory for one cell, json filename for cell separation")
} else {
    kallisto.dir  <- args[1]
    filtered.json <- args[2]
    filterCells(kallisto.dir, filtered.json)
}
