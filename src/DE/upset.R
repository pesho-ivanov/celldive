suppressWarnings(library(logging))
basicConfig()  # for logging

suppressWarnings(library(UpSetR))
suppressWarnings(library(ggplot2))

PlotUpset <- function(gene.sets, name, outfile) {
    if (length(gene.sets) == 0) {
        cat("Warning: ", name, " zero lenght\n")
        return
    } else {
        df <- fromList(gene.sets)  # an UpSet function

        #gene.sets <- gene.sets[1:2]
        #df <- do.call(rbind.data.frame, gene.sets)
        #print(length(gene.sets[["DE_DESeq2_WaG_SE"]]))
        #gene.sets <- gene.sets[ lapply(gene.sets, length) > 0 ]
        #print(length(gene.sets))
        #print(df)
        cat(name, dim(df), "\n")
        p <- upset(df,
              order.by = c("degree", "freq"),
              nsets=nrow(df),
              #nintersects = 60,
              #scale.intersections = "log2",
              keep.order = TRUE,
              sets = names(df),
              mainbar.y.label = "Common genes",
              sets.x.label = "#DE genes Clonal vs Nonclonal ",
    #          empty.intersections = "on",
              #mb.ratio = c(0.3, 0.7),
        )
        #ggsave(filename=paste0(dir.out, "DE_", name, ".png"), plot=p)
        #ggsave(filename=outfile, plot=p) #, device="png")
        loginfo(sprintf("UpSet diagram saving in %s...", outfile))
        plot(p)
        #png(outfile)
        #  plot(p)
        #ggplot(plot=p)#, device=png)
        #dev.print(file=outfile, device=png, width=800)
        dev.off()
    }
}

GetSubset <- function(set, grep.rule) {
    subset <- list()
    for (name in names(set)) {
        if (grepl(grep.rule, name)) {
            subset[[name]] <- set[[name]]
        }
    }
    return(subset)
}

SubsetAndUpset <- function(gene.sets, regexps, outfile) {
    for (regexp in regexps) { 
        subset <- GetSubset(gene.sets, regexp)
        PlotUpset(subset, name=regexp, outfile=outfile)
    }
}

# TODO make ReadGeneset
# TODO limit p_val and p_val_adj
ReadGenesets <- function(dir, recursive=FALSE) {
    files <- list.files(path=dir, pattern="*.csv", full.names=T, recursive=recursive)
    loginfo(sprintf("Found %d DE tables in %s", length(files), dir))
    gene.sets <- list() #data.frame()
    for (f in files) {
        t <- read.csv(f, header=T, stringsAsFactors=FALSE) # load file
        #print(t)
        if (nrow(t) > 0) {
            genes <- t$symbol
            name <- tools::file_path_sans_ext(basename(f))
            gene.sets[[name]] <- genes
        }
    }
    loginfo(sprintf("Sourced %d DE tables:", length(gene.sets)))
    print(names(gene.sets))
    return(gene.sets)
}

Dir2UpSet <- function(dir, outfile) {
    gene.sets <- ReadGenesets(dir)
    SubsetAndUpset(gene.sets, c(""), outfile)
}
        
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    warning("de.R: Not exact number of arguments. Expected 2 args: directory with gene tables, outfile(for a DE table)")
} else {
    dir_with_de_tables <- args[1]
    outfile <- args[2]
    Dir2UpSet(dir_with_de_tables, outfile)
}
