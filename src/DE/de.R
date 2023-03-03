# 1) https://hemberg-lab.github.io/scRNA.seq.course/
# *) Why Kallisto:
#    https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/
# 2) SCDE: http://hms-dbmi.github.io/scde/diffexp.html
# 3) tximport (transcripts->genes):
#      http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 4) Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences (2016):
#      https://f1000research.com/articles/4-1521/v2
# 5) DE methods comparison: SCDE and ROTS are good but work on counts data; MAST, DESeq, Limma are shit
#      https://academic.oup.com/bib/article/18/5/735/2562772

options(warn=-1)

suppressPackageStartupMessages({
  library(logging)
  library(SingleCellExperiment)
  #library(Seurat)
  #library(mclust)
  library(scater)
  library(rjson)
  library(scde)
  library(ROTS)
  library(edgeR)
  library(dplyr)
  library(tibble)
})
basicConfig()  # for logging

#suppressMessages(source("seurat.R"))
source("DE/kallisto2genetable.R")
source("DE/vis.R")

#filterUnsupervised_scater <- function(sce, min.lib.size, min.reads, min.detected) {
#    # min.lib.size -- cells with less genes
#    # min.reads -- genes with less reads are removed
#    # min.detected -- genes seen in less cells are removed
#    old_ma3x <- counts(sce)
#    loginfo(sprintf("Unfiltered matrix with %d features in %d cells.", 
#                    nrow(old_ma3x), ncol(old_ma3x)))
#    loginfo(sprintf("Filtering with min.lib.size=%d, min.reads=%d, min.detected=%d.",
#                    min.lib.size, min.reads, min.detected))
#    ma3x <- clean.counts(old_ma3x, min.lib.size=min.lib.size, min.reads=min.reads, min.detected=min.detected)  # SCDE, filter genes
#    sce_new <- sce
#    counts(sce_new) <- ma3x
#    loginfo(sprintf("Filtered to %d features in %d cells.", nrow(ma3x), ncol(ma3x)))
#    return(sce_new)
#}

filterUnsupervised <- function(sce) {
    #keep.total <- sce$total_counts > 1e5
    #keep.n <- sce$total_features_by_counts > 50
    #sce <- sce[, keep.total]
    is.sce <- tpm(sce) >= 10
    rows <- rowSums(is.sce) >= 0.25 * ncol(is.sce) #5   # each gene should be seen in at least N cells
    #rows <- rowMeans(counts(sce)) >= 10
    #filter <- means >= 10
    #cols <- colSums(is.sce) > 5   # each cell should include at least N expressed genes
    #sce <- sce[rows, cols]
    sce <- sce[rows, ]

    loginfo(sprintf("Filtered to %d features in %d cells.", nrow(sce), ncol(sce)))
    return(sce)
}

sceSubset <- function(sce, gene.list.fn) {
      loginfo(sprintf("Reading the predefined gene list for DE tests from %s", gene.list.fn))
      genes <- read.csv(gene.list.fn)
      m <- match(genes$symbol, rownames(sce))
      idxs <- m[!is.na(m)]
      sce <- sce[ idxs, ] ## left_join(results, genes$symbol, by="symbol")
}

#TMM_normalize <- function(sce) {
#    loginfo(sprintf("Normalizing by TMM."))
#    counts <- counts(sce)
#    cds <- calcNormFactors(counts, method="TMM")
#    scale <- cds["samples"]["lib.size"] * cds["samples"]["norm.factors"]
#    normCounts <- round(t(t(counts)/scale) * mean(scale))
#
#    #sce_new <- sce
#    normcounts(sce) <- normCounts
#    return(sce)
#}

substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
}

checkDuplicates <- function(bla) {
    n_occur <- data.frame(table(bla))
    #print(length(n_occur[n_occur$Freq > 1,]))
    #print(head(n_occur[n_occur$Freq > 1,]))
    if (sum(n_occur$Freq > 1) > 0) {
        logerror("Duplicating genes!")
        stop()
    }
    #print(head(bla[bla %in% n_occur$Var1[n_occur$Freq > 1],]))
    # /tmp
}
                                           
DE_run_and_output <- function(sce, method, sample_name, outdir, de_table_file, gene.list.fn) {
    loginfo(sprintf("Differential expression on %d cells with %d features using method=%s",
            ncol(sce), nrow(sce), method))

    dir.create(outdir, recursive=TRUE)
    outfile <- file.path(outdir, paste0("/DE_", sample_name, "_", method, ".pdf"))
    groups <- colData(sce)$groups
    print(groups)
    #groups <- factor(substrRight(colnames(sce), 1), levels = c("N", "C"))
    #loginfo(sprintf("DE between %d clonal and %d non-clonal cells.",
    #                sum(groups == "C"), sum(groups == "N")))
    #loginfo(table(groups))

    if (length(groups) != ncol(sce)) {
        logerror("Wrong sizes.")
        stop()
    }
    if (gene.list.fn != 'all') {
        sce <- sceSubset(sce, gene.list.fn)
    }
    loginfo(sprintf("Running %s DE method...", method))
    if (method == "ROTS") {
        # runtime ~1-2min on all de samples
        # works on integer counts
        # produces "logFC" which is log2

        # Conquer pipeline
        # grp <- L$condt
        # dge <- DGEList(counts = L$count)
        # dge <- edgeR::calcNormFactors(dge)
        # cpms <- cpm(dge)
        # rots <- ROTS(data = cpms, groups = grp, B = 1000, K = 1000, log = FALSE, seed = 123)

        fdr <- 0.05
        B <- 300
        K <- 300

        #sce <- TMM_normalize(sce)  ## TODO: move before subsetting
        dge <- DGEList(counts = counts(sce))
        dge <- edgeR::calcNormFactors(dge)
        ma3x <- cpm(dge)
        # ma3x <- counts(sce)
        rots <- ROTS(data=ma3x, groups=(levels(groups)[[2]]), B=B, K=K, seed=1234, log=FALSE, progress=FALSE)

        #summary(rots, fdr=fdr)
        results <- data.frame(summary(rots, fdr=fdr))
        results <- rownames_to_column(results, "symbol") #print(head(rots$logfc))
        logfc <- as.data.frame(rots$logfc)
        colnames(logfc) <- "logFC"
        logfc <- rownames_to_column(logfc, "symbol")
        results <- left_join(results, logfc, by="symbol")
        results["pvalue_adj"] <- results["FDR"]
        
        # visualize
        loginfo(sprintf("Writing ROTS DE graphs to %s", outfile))
        #outfile <- paste0(de_table_file, ".pdf")
        pdf(outfile)
            plot(rots, fdr=fdr, type="volcano")
            plot(rots, fdr=fdr, type="heatmap")
            plot(rots, fdr=fdr, type="ma")
            plot(rots, fdr=fdr, type="reproducibility")
            plot(rots, fdr=fdr, type="pvalue")
            plot(rots, fdr=fdr, type="pca")
        dev.off()

    } else if (method == "SCDE") {
        # 22min on all de samples
        # TODO: fdr control
        # TODO: cZ -> pvalue
        # TODO: speed-up?
        # requires counts; developed for 8+8 cells

        # TODO: move filtering earlier
        #cd <- clean.counts(ma3x, min.lib.size=1000, min.reads = 1, min.detected = 1)

        cd <- counts(sce)
        storage.mode(cd) <- "integer" #mode(cd) <- "integer" 
        randomizations <- 10 # 100
        cores <- 25

        # fit a model independently for both groups
        # tutorial: http://hms-dbmi.github.io/scde/diffexp.html
        o.ifm <- scde.error.models(
            counts = cd,
            groups = groups,
            n.cores = cores,
            threshold.segmentation = TRUE, save.crossfit.plots = FALSE, 
            save.model.plots = FALSE, verbose = 0)
        o.prior <- scde.expression.prior(
            models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
        ediff <- scde.expression.difference(
            o.ifm, cd, o.prior, groups = groups, n.randomizations = randomizations,
            n.cores = cores, verbose = 0)
        results <- ediff[order(ediff$Z, decreasing = TRUE), ] 
        results <- rownames_to_column(results, "symbol")
        results <- results[results$ce != 0, ]
        results$logFC <- results$ce
    } else if (method == "edgeR") {
        # 1-2min on all de samples
        # Normalization is only necessary for sample-specific effects"
        # produces "logFC" which is log2

        # Conquer pipeline:
        # dge <- DGEList(L$count, group = L$condt)
        # dge <- calcNormFactors(dge)
        # design <- model.matrix(~L$condt)
        # dge <- estimateDisp(dge, design = design)
        # fit <- glmQLFit(dge, design = design)
        # qlf <- glmQLFTest(fit)
        # tt <- topTags(qlf, n = Inf)

        cd <- counts(sce)
        storage.mode(cd) <- "integer" #mode(cd) <- "integer" 
        y <- DGEList(counts=cd, group=groups)
        y <- calcNormFactors(y)
#        print("groups:")
#        print(groups)
        design <- model.matrix(~groups)
#        print("design:")
#        print(design)
#        print(y)
        #print(nrow(design))
        #print(ncol(y))
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design)

        if (length(levels(groups)) == 2) {
            qlf <- glmQLFTest(fit) #, coef=2)
        } else if (length(levels(groups)) == 4) {
# earlyclonal is a reference
            contrasts <- makeContrasts((groupslateclonal)-(groupslatesingle-groupsearlysingle), levels=design)
            qlf <- glmQLFTest(fit, contrast=contrasts) #=c(1,-1,-1,1)) #, coef=2)
        } else {
            logerror("Wrong number of groups: %d", length(levels(groups)))
            stopifnot(FALSE)
        }
        #quasi-likelihood F-tests
        tt <- topTags(qlf, n=Inf, p=0.05)
        results <- as.data.frame(tt)
        if (ncol(results)) {
          results <- rownames_to_column(results, "symbol")
          #colnames(results)[colnames(results) == 'PValue'] <- 'pvalue'
          results["pvalue"] <- results["PValue"]
          #results$logFC <- -results$logFC  # so that more expression in clonal is positive
        }

        edger.fn <- file.path(outdir, paste0("/DE_", sample_name, "_", method, "_DGEList.rds"))
        #saveRDS(y, file=edger.fn)
        loginfo(sprintf("Writing EdgeR DE graphs to %s", outfile))
        outfile <- paste0(de_table_file, ".pdf")

        pdf(outfile)
            #pch <- c(0,1,2,15,16,17)
            plotMDS(y, col=c("red", "blue")) #, pch=pch[group])

            plotMD(y, column=1)
            abline(h=0, col="red", lty=2, lwd=2)
            plotMD(y, column=2)
            abline(h=0, col="red", lty=2, lwd=2)

            plotBCV(y)
            plotQLDisp(fit)
            #plotMD(y, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
            #legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)
        dev.off()
        
        #plotVidger(y, sample_name, outdir)
    } else {
        logerror("Wrong DE method name: %s", method)
        stopifnot(FALSE)
    }

    loginfo(sprintf("Found %d differentially expressed genes.", ncol(results)))

    #print(head(results))
    cnames <- colnames(results)
    #if (!"pvalue_adj" %in% cnames) {
    #    logerror(sprintf("pvalue_adj column not among the produced columns for %s", method))
    #    print(cnames)
    #}

    #results <- results[results$FDR <= 0.05, ]
    checkDuplicates(results$ensg)

    grch38_unique <- grch38[ !duplicated(grch38$symbol), ]

    if (ncol(results)) {
      results <- left_join(results, grch38_unique, by="symbol") %>%
            select(symbol, everything())
            #select(ensg, symbol, description)
    }
    loginfo(sprintf("DE table written to %s", de_table_file))
    #results["tested genes"] <- rownames(sce)
    write.csv(results, de_table_file, row.names = FALSE)

    checkDuplicates(results$ensg)

    loginfo("Differential expression done.")
    return(results)
}

printGenes <- function(sce, fn, sample_name) {
    fn <- file.path(out_dir, sprintf("genes_after_filtering_%s.csv", sample_name))
    df <- data.frame(rownames(counts(sce)))
    colnames(df) <- sample_name
    loginfo(sprintf("Writing %d gene names to %s", nrow(df), fn))
    write.csv(df, fn, row.names = FALSE)
    return(sce)
}
                                           
DE_e2e <- function(sample_name, clonality_file, method, out_dir, de_table_file, clear.groups, gene.list.fn, config.fn) {
    loginfo(sprintf("Differential expression for %s...", sample_name))
    if (!file.exists(config.fn)) {
        warning(paste0("File config.fn=", config.fn, "doesn't exist."))
    }

    sce <- clonality2sce(clonality_file, sample_name, clear.groups, config.fn) %>%
                        #sce2file(paste0(de_table_file, ".matrix")) %>%
    # TODO: filter out SpikeIns
    #sce <- readSamples_old(kallisto_dir, clonality) %>%
    #                   readSamples(clonality) %>%
    #sce_tx    <- pair[["sce"]]
    #groups    <- pair[["groups"]]
    #                   tx2geneid %>%
                       calculateQCMetrics %>%
    #                   geneid2genesymbol %>%

    #sce_gene        <- geneid2genesymbol(sce_gene)   #sce_tx    <- addAnnotsBM(sce_tx)
    #sce_filtered    <- filterUnsupervised_scater(sce_gene, 100, 5, 3)
                       filterUnsupervised()
    # printGenes(sce, out_dir, sample_name) %>%

    results <- DE_run_and_output(sce, method, sample_name, out_dir, de_table_file, gene.list.fn)
    return(results)
}
                                           
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    warning("de.R: Not exact number of arguments. Expected 7 args: sample_name(01_SE), clonality_file(json), method(ROTS, SCDE), out_dir(all all kinds of pictures), outfile(for a DE table), clear.groups(1 or 0), file with gene lists as a first column or 'all'")
} else {
    sample_name    <- args[1]
    clonality_file <- args[2]
    method         <- args[3]
    out_dir        <- args[4]
    de_table_file  <- args[5]
    clear.groups   <- args[6]
    gene.list.fn   <- args[7]
    config.fn      <- args[8]
    print(config.fn)
    result         <- DE_e2e(sample_name, clonality_file, method, out_dir, de_table_file, clear.groups, gene.list.fn, config.fn)
}
