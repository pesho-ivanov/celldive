library(Seurat)
#library(annotables)

#library(cowplot)  # cowplot enables side-by-side ggplots
library(dplyr)

source("annots.R")

options(stringsAsFactors = FALSE)

###########################################
############# batches #####################
###########################################
SchH.good <- c('SchH_SE', 'SchH_SL')
WaG.good <- c('WaG_BL-C', 'WaG_SE')
FrK.good <- c('FrK_BL', 'FrK_SE', 'FrK_SL')
#batches.bad <- c('SchH_BE', 'WaG_BE', 'WaG_BL' 'WaG_SL')

nGene.limits        = list(SchH_SE = c(4000, 7000), SchH_SL = c(), "WaG_BL-C" = c(0, Inf), WaG_SE = c(4000, 7000), FrK_BL = c(2500, 7000), FrK_SE = c(5000, 9000), FrK_SL = c(4000, 7000))
nUMI.limits         = list(SchH_SE = c(986000, 1000000), SchH_SL = c(), "WaG_BL-C" = c(0, Inf), WaG_SE = c(),      FrK_BL = c(990000, Inf), FrK_SE = c(), FrK_SL = c(985000, Inf))
percent.mito.limits = list(SchH_SE = c(0,1), SchH_SL = c(0, 0.15), "WaG_BL-C" = c(0, 0.3), WaG_SE = c(),           FrK_BL = c(0, 0.12), FrK_SE = c(0, 0.17), FrK_SL = c(0, 0.11))

batches.good <- c(SchH.good, WaG.good, FrK.good)  

#########################################
##### genes  ############################
#########################################
Park.genes = c("TP53", "POT1", "ATM", "ZNG365", "PLCG1", "CARD11", "CD28", "RLTPR", "PRKG1", "PTPRN2", "IRF4", "PDCD1", "PRKCQ", "ARID1A", "TRRAP",
              "DNMT3A", "TET2", "CREBBP", "KMT2D", "NCOR1", "SMARCB1", "BCOR", "KMT2C", "CTCF", "KDM6A", "SETDB2", "RHOA", "ARHGEF3", "ZEB1",
              "RARA", "STAT5B", "JAK3", "JAK1", "STAT3", "JAK2", "CCR4", "FAS", "RFXAP", "CD58", "TNFRSF1B", "NFKB2", "CSNK1A1",
              "TNFAIP3", "PRKCB", "VAV1", "PIK3R1", "U2AF1", "NF1", "MAPK1", "BRAF", "KRAS", "MAPK2K1", "NRAS", "CDKN2A", "LATS1")
Park.genes.df <- data.frame(Park.genes, stringsAsFactors=FALSE) %>%
    inner_join(grch38, by=c("Park.genes"="symbol")) %>%
    select(Park.genes, ensgene, description) 

# by hand
clonal.vs.nonclonal.DESeq2.findings <- c(
    'HPGDS',   # FrK_SE, FrK_SL, SchH_SE, SchH_SL
    'PLS3',    # FrK_BL, FrK_SE, FrK_SL
    'CCL5',    # FrK_SL, SchH_SL
    'GZMK',    # FrK_SL, SchH_SL
    'HLA-DRB1',# FrK_SL, SchH_SL
    'MTG1',    # FrK_BL, FrK_SL
    'GPR15',   # FrK_BL, FrK_SE
    'CCL3'     # FrK_SL, SchH_SL
    )




###########################################
#######  I/O  #############################
###########################################
ReadTable <- function(fn) {
    print(paste0("Reading ", fn, "..."))
    table <- read.table(file = fn, sep = "\t", header=TRUE, row.names="target_id")
    return(table)
}

isSingleString <- function(input) {
    is.character(input) & length(input) == 1
}

ReadBatchPart <- function(fn) {
    #if (!b.data) {
    #    stop(paste0(fn, "does not exist"))
    #}
    if (isSingleString(fn)) {
        table <- ReadTable(fn)
    } else {
        table <- fn
    }
    geneids <- tids2geneids(table)
    symbols <- geneids2symbols(geneids)
    cat(nrow(table), "transcripts", "->", nrow(symbols), "genes", "\n")
    batch <- CreateSeuratObject(raw.data = symbols)
    return(batch)
}

#plot_violin <- function(dir.out, batch, name) {
#    p <- VlnPlot(object = batch, features.plot = c("nGene"), nCol = 1)
#    dir.create(dir.out, showWarnings = FALSE)
#    ggsave(filename=paste0(dir.out, "violin_", name, "clonal.png"), plot=p)
#}

AddMeta <- function(batch) {
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = batch@data), value = TRUE)
    percent.mito <- Matrix::colSums(batch@raw.data[mito.genes, ])/Matrix::colSums(batch@raw.data)

    # AddMetaData adds columns to object@meta.data, and is a great place to
    # stash QC stats
    batch <- AddMetaData(object = batch, metadata = percent.mito, col.name = "percent.mito")
    batch
}

GetLimits <- function(limits, name) {
    if (name %in% names(limits))
        if (length(limits[[name]]) > 0)
            return(limits[[name]])
    return(c(-Inf,+Inf))
}

FilterCellsInBatches <- function(batches, dir.out, nGene.lim, percent.mito.lim, nUMI.lim) {
    new.batches <- list()
    for (batch.name in names(batches)) {
        batch <- batches[[batch.name]] %>% AddMeta()
        nGene <- GetLimits(nGene.lim, batch.name)
        percent.mito <- GetLimits(percent.mito.lim, batch.name)
        nUMI <- GetLimits(nUMI.lim, batch.name)

        #print(percent.mito)
        #print(batches[[batch.name]]@meta.data)
        #head(x = FetchData(object = nGene, vars.all = 'percent.mito' ))

        print(batch)
        new.batch <- FilterCells(object = batch,
                             subset.names = c("nGene"), 
                             low.thresholds = c(nGene[[1]]),
                             high.thresholds = c(nGene[[2]]))
#        new.batch <- FilterCells(object = batch,
#                             subset.names = c("nGene", "percent.mito", "nUMI"), 
#                             low.thresholds = c(nGene[[1]], percent.mito[[1]], nUMI[[1]]),
#                             high.thresholds = c(nGene[[2]], percent.mito[[2]], nUMI[[2]]))
#
        PlotQc(dir.out, new.batch, batch.name)
        new.batches[[batch.name]] <- new.batch
    }
    return(new.batches)
}

ReadBatch <- function(dir.in, dir.out, batch_name, unknown.as.clonal=FALSE, unknown.as.nonclonal=FALSE, nonreconstructed.as.nonclonal=FALSE) {
    cl <- ReadBatchPart(paste0(dir.in, batch_name, "_clonal.csv"))
    noncl <- ReadBatchPart(paste0(dir.in, batch_name, "_nonclonal.csv"))
    
    if (unknown.as.clonal) {
        unk <- ReadBatchPart(paste0(dir.in, batch_name, "_unknown.csv"))
        cl <- MergeSeurat(object1 = cl, object2 = unk, project = batch_name)
    }
    
    if (unknown.as.nonclonal) {
        unk <- ReadBatchPart(paste0(dir.in, batch_name, "_unknown.csv"))
        noncl <- MergeSeurat(object1 = noncl, object2 = unk, project = batch_name)
    }
    
    if (nonreconstructed.as.nonclonal) {
        nonrec <- ReadBatchPart(paste0(dir.in, batch_name, "_nonreconstructed.csv"))
        noncl <- MergeSeurat(object1 = noncl, object2 = nonrec, project = batch_name)
    }
    
    batch <- MergeSeurat(object1 = cl, object2 = noncl, project = batch_name)
    batch <- SetIdent(batch, cells.use = WhichCells(object = cl), id = 'clonal')
    batch <- SetIdent(batch, cells.use = WhichCells(object = noncl), id = 'nonclonal')
    
    batch <- NormalizeData(object = batch)
    batch <- ScaleData(object = batch, display.progress = FALSE)
    batch <- FindVariableGenes(object = batch, do.plot = FALSE)
    
    PlotQc(dir.out, batch, batch_name)
    return(batch)
}

ReadBatches <- function(batch.names, dir.in, dir.out, unknown.as.clonal=FALSE, unknown.as.nonclonal=FALSE, nonreconstructed.as.nonclonal=FALSE) {
    batches <- list()
    for (batch.name in batch.names) {
        batches[[batch.name]] <- ReadBatch(dir.in, dir.out, batch.name,
                                           unknown.as.clonal=unknown.as.clonal,
                                           unknown.as.nonclonal=unknown.as.nonclonal,
                                           nonreconstructed.as.nonclonal=nonreconstructed.as.nonclonal)
        print(batches[[batch.name]])
    }
    return(batches)
}

read_early_late_clonal <- function(dir.in, dir.out, early.batch, late.batch) {
    early.cl <- ReadBatchPart(paste0(dir.in, early.batch, "_clonal.csv"))
    late.cl <- ReadBatchPart(paste0(dir.in, late.batch, "_clonal.csv"))
    
    batch <- MergeSeurat(object1 = early.cl, object2 = late.cl, add.cell.id1 = "early", add.cell.id2 = "late", project = paste0(early.batch, " ", late.batch))
    #print(levels(batch@ident))
    return#(batch)
    batch <- SetIdent(batch, cells.use = WhichCells(object = early.cl), id = 'early')
    batch <- SetIdent(batch, cells.use = WhichCells(object = late.cl), id = 'late')
    
    plot_qc(dir.out, batch, paste0(early.batch, "_", late.batch))
    return(batch)
}


        
##################################################################
###################  Plotting  ###################################
##################################################################
PlotQc <- function(dir.out, batch, name) {
    batch <- AddMeta(batch)
    p <- VlnPlot(object = batch, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

    #par(mfrow = c(1, 2))
    #reads.mito.p <- GenePlot(object = batch, gene1 = "nUMI", gene2 = "percent.mito")
    #reads.genes.p <- GenePlot(object = batch, gene1 = "nUMI", gene2 = "nGene")
    
    dir.create(dir.out, showWarnings = FALSE)
    ggsave(filename=paste0(dir.out, "qc_", name, ".png"), plot=p)
    #ggsave(filename=paste0(dir.out, "qc_reads_mito_", name, ".png"), plot=reads.mito.p)
    #ggsave(filename=paste0(dir.out, "qc_reads_genes_", name, ".png"), plot=reads.genes.p)
}

PlotGenes <- function(batch, genes, fn) {
    cat(length(genes)," genes to be plot.\n")
    if (length(genes) > 0) {
        p <- VlnPlot(object = batch, nCol = 4, features.plot = genes, do.return=TRUE) 
        cat("Writing violins: ", fn, "\n")
        ggsave(filename=fn, plot=p)
    }
}




#################################################################
######################### DE ####################################
#################################################################
# •  "wilcox" : Wilcoxon rank sum test (default)
# •  "bimod" :  Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)
# •  "roc" : Standard AUC classifier
# •  "t" : Student’s t-test
# •  "tobit" :  Tobit-test for differential gene expression (Trapnell et al., Nature Biotech, 2014)
# •  "poisson" :  Likelihood ratio test assuming an underlying poisson distribu- tion. Use only for UMI-based datasets
# •  "negbinom" :  Likelihood ratio test assuming an underlying negative bino- mial distribution. Use only for UMI-based datasets
# •  "MAST : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015)
# •  "DESeq2 :  DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014)
de.tests <- c("wilcox", "DESeq2", "bimod", "t", "tobit", "poisson")

DE <- function(batch, test, id1="clonal", id2="nonclonal", genes.use=NULL, max.p_val = 0.05, max.p_val_adj = 0.9999) {
    options(warn=-1)
    print(paste0(test, " on ", batch@project.name, "..."))
    markers <- FindMarkers(object = batch,
                           test.use = test,
                           ident.1 = id1,
                           ident.2 = id2,
                           min.pct = 0.25,
                           genes.use = genes.use)
    #markers <- add_rownames(markers, "enst")
    #markers <- add_symbol_name(markers)
    #markers <- markers %>%
    #    tibble::rownames_to_column("symbol") %>%
    #    filter(p_val <= max.p_val) %>%
    #    filter(p_val_adj <= max.p_val_adj)
    #print(x = head(x = markers, n = 10))
    file_path <- paste0(dir.out, "DE_", test, "_", batch@project.name, ".csv")
    write.csv(markers, file = file_path, row.names = TRUE)
    options(warn=0)
    return(markers)
}

RunOnAllBatches <- function(batches, de.tests, name, genes.use = NULL, dir.out) {
    for (test in de.tests) {
        for (batch.name in names(batches)) {
            batch <- batches[[batch.name]]
            markers <- DE(batch, test)
            filtered.markers <- FilterGenes(markers, genes.use, max.p_val=0.3)
            fn <- paste0(dir.out, "violins_", name, "_", batch.name, ".png")
            PlotGenes(batch, genes=filtered.markers$symbol, fn=fn)
        }
    }
}



#add_symbol_name <- function(markers) {
#    markers %>%
#        add_rownames("enst") %>%
#        arrange(p_val_adj) %>% 
#        head(20) %>% 
#        inner_join(grch38_tx2gene, by=c("enst"="enstxp")) %>%
#        inner_join(grch38, by=c("ensgene"="ensgene")) %>% 
#        select(enst, p_val, avg_logFC, pct.1, pct.2, p_val_adj, symbol) %>% 
#        as.data.frame()
#}

FilterOutTCRGenes <- function(markers) {
    df <- data.frame(symbol=markers) %>%
        #mutate_if(is.factor, as.character) %>%
        #tibble::rownames_to_column("symbol") %>%
        filter(!grepl("^(TRA|TRB)", symbol))
    cat(length(markers), "->", nrow(df), ", ")
    return(df$symbol)
}
        
FilterGenes <- function(markers, genes, max.p_val=0.05, max.p_val_adj=1.0) {
    if (! "symbol" %in% names(markers)) {
        markers <- markers %>% tibble::rownames_to_column("symbol")
    }
    filtered.markers <- markers %>%
        filter(symbol %in% genes) %>%
        filter(p_val <= max.p_val) %>%
        filter(p_val_adj <= max.p_val_adj)
    #print(filtered.markers)
    return(filtered.markers)

#    df <- table %>%
#        #mutate_if(is.factor, as.character) %>%
#        #tibble::rownames_to_column("symbol") %>%
#        filter(symbol %in% genes)
#    cat(length(markers), "->", nrow(df), ", ")
#    return(df)
}
        
#
#dir.in <- "/home/pesho/repos/bioviz/in/clonality_tables/"
#dir.out <- "/home/pesho/repos/bioviz/out/seurat/"
#
#batch <- ReadBatch(dir.in, dir.out, 'FrK_BL')
##batch <- FilterCellsInBatches(batches, dir.out, nGene.limits, nUMI.limits, percent.mito.limits)
#markers <- DE_all(batch, tests=c("DESeq2"))
#
#source("tooling_scripts/seurat.R")
#
#tests <- de.tests
## "roc",  # no "p_val" column
## "negbiom", # Error in rownames(x = to.return): object 'to.return' not found
##tests <- c("DESeq2")
#DE_all(batches, tests=tests)  #'DESeq2'
#
