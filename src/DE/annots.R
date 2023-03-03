suppressMessages(library(annotables))
suppressMessages(library(scater))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

library(logging)


# DON'T USE
# using Scater, BioMart -- ONLINE
#addAnnotsBioMart <- function(sce) {
#    loginfo("Loading BioMart annotations...")
#    #flush.console()
#    sce <- getBMFeatureAnnos(
#        sce,
#        filters = "ensembl_transcript_id", 
#        attributes = c(
#            "ensembl_transcript_id",
#            "ensembl_gene_id",
#            "hgnc_symbol",
#            "chromosome_name",
#            "start_position",
#            "end_position"
#        ), 
#        feature_symbol = "hgnc_symbol",
#        feature_id = "hgnc_symbol",  # was ensembl_gene_id
#        biomart = "ENSEMBL_MART_ENSEMBL", 
#        dataset = "hsapiens_gene_ensembl",
#        host = "www.ensembl.org")
#    sce
#}


# changes the rownames from ENSTs to ENSGs
#tids2geneids <- function(table) {
#    cells <- names(table)
#    table %>%
#        tibble::rownames_to_column("enstxp") %>%
#        inner_join(grch38_tx2gene, by=c("enstxp")) %>%
#        #arrange(ensgene) %>%
#        #head(20) %>%
#        group_by(ensgene) %>%
#        summarise_if(is.numeric, sum) %>% #summarise_each(funs(sum), cells)
#        tibble::column_to_rownames('ensgene') %>%
#        as.data.frame()
#}
#
## changes the row names from ENSG to gene symbols
## ex: ENSG00000000003 -> TSPAN6
#geneids2symbols <- function(table) {
#    cells <- names(table)
#    table %>%
#        tibble::rownames_to_column("ensgene") %>%
#        inner_join(grch38, by="ensgene") %>%
#        select(symbol, cells) %>%
#        group_by(symbol) %>%
#        summarise_if(is.numeric, sum) %>%
#        tibble::column_to_rownames('symbol') %>%
#        as.data.frame()
#}

# using Scater
aggregate_tx2genes <- function(sce) {
    #isSpike(sce, "Spike-") <- grepl("Spike-", rownames(sce)) # which(isSpike(sce, "Spike-"))
    
    loginfo("Aggregate transcripts into genes...")

    # TODO: not aggregate over counts?!
    # https://www.rdocumentation.org/packages/scater/versions/1.0.4/topics/summariseExprsAcrossFeatures
    sce_gene <- summariseExprsAcrossFeatures(sce, exprs_values="counts", summarise_by="feature_id")  # scater
    loginfo(sprintf("%d transcripts aggregated into %d genes.", nrow(sce), nrow(sce_gene)))
    if (nrow(sce_gene) == nrow(sce))
        logwarn("The aggregation did not lower the number of features!")
    return(sce_gene)
}

tx2geneid <- function(sce) {
    loginfo("tx2geneid")
    
    # tid -> gid
    rn <- rownames(sce)
    df <- data.frame(rn) %>%
          left_join(grch38_tx2gene, by=c("rn"="enstxp")) %>%
          select(ensgene)
    #df <- replace_na(c(df), rn)
    nas <- which( is.na(df$ensgene) )
    loginfo(sprintf("%d transcripts were not converted to genes but left as they are.", length(nas)))
    df$ensgene[nas] <- rn[nas]
    #rownames(sce) <- df$ensgene #[["symbol"]]
    rowData(sce)[["feature_id"]] <- df[["ensgene"]]
    
    # aggregate same genes
    sce_gids_unique <- aggregate_tx2genes(sce)
    print("after tx2geneid")
    print(sce_gids_unique)
    return(sce_gids_unique)
}

geneid2genesymbol <- function(sce) {
    loginfo("geneid2genesymbol")
    if (sum(is.na(rownames(sce))) > 0) {
        logerror(sprintf("%d NA geneid", sum(is.na(rownames(sce)))))
        stop()
    }
    
    # unique gid -> symbol
    rn <- rownames(sce)
    df <- data.frame(rn) %>%
          left_join(grch38, by=c("rn"="ensgene")) #%>%
          #select(symbol) #%>%
          #distinct(symbol) #, .keep_all = TRUE)
    if (sum(is.na(df$symbol)) > 0) {
        logerror(sprintf("tmp: %d NA genesymbols", sum(is.na(df$symbol))))
        stop()
    }
    nas <- is.na(df$symbol)
    df$symbol[nas] <- df$rn[nas]
    loginfo(sprintf("%d geneids were not converted to gene symbols but left as they are.", length(nas)))
    #df$symbol[nas] <- rn[nas]
    #rownames(sce) <- df$symbol #[["symbol"]]
    if (sum(is.na(rownames(sce))) > 0) {
        logerror(sprintf("%d NA genesymbols", sum(is.na(rownames(sce)))))
        stop()
    }

    #print("after geneid2genesymbol")
    #print(sce)
    return(sce)
}

