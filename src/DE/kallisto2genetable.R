suppressMessages(suppressWarnings(library(tximport)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(library(annotables))
suppressMessages(library(dplyr))

clonality2triple <- function(clonality.file, sample.name, clear.groups, config.fn) {
    #getSuffix <- function(x, n) {
    #    return(substr(x, nchar(x)-n+1, nchar(x)))
    #}
    
    removeSamplePrefix <- function(names) {
        prefix <- paste0("^", sample.name, "_")
        return(sub(prefix, "", names))
    }
    
    loginfo(sprintf("Reading the %s's clonality file %s...", sample.name, clonality.file))
    json_data      <- fromJSON(file=clonality.file)
    if (length(json_data) != 1) {
      logerror("The upper level of the clonality.file JSON should have exactly one element with a cell information.")
      stop()
    }

    cell_partition <- json_data[[1]]

    loginfo(sprintf("Reading info about the compared classes from config %s...", config.fn))
    json.config    <- fromJSON(file=config.fn)
    #loginfo(json.config)
    #if (length(json.config) != 2) {
    #  logerror("The upper level of the config.fn JSON should have an array with exactly two elements: two arrays with classes for each of the two DE groups.")
    #  stop()
    #}

    first_classes <- json.config[[1]]  # list of classes
    second_classes <- json.config[[2]]  # list of classes

    #loginfo(sprintf("Reading expressions for first group [%s] and from the second group [%s]...", first_classes, second_classes))
    loginfo("Reading expressions for groups: ", json.config)
    #print(json.config)

    getCells <- function(group) {
      group <- cell_partition[[group]]
      cells <- removeSamplePrefix(names(group))
      return(cells)
    }

    getKallistoDirs <- function(group) {
      group <- cell_partition[[group]]
      dirs <- lapply(group, function(cell) file.path(cell[["kallisto_subdir"]], "abundance.h5"))
      return(dirs)
    }

    #cell_partition <- json_data[["cell_partition"]]
    #X_cells <- getCells(first_classes[[1]])
    #X_dirs <- getKallistoDirs(first_classes[[1]])
    #Y_cells <- getCells(second_classes[[1]])
    #Y_dirs <- getKallistoDirs(second_classes[[1]])

    cells <- sapply(json.config, getCells)
    dirs <- sapply(json.config, getKallistoDirs)
    #print(cells)
    #print(dirs)

    # TODO
    #if (clear.groups == 0) {
    #    X_plus_cells <- getCells("main+")
    #    X_plus_dirs <- getKallistoDirs("main+")
    #    X_cells <- c(X_cells, X_plus_cells)
    #    X_dirs <- c(X_dirs, X_plus_dirs)

    #    Y_plus_cells <- getCells("small")
    #    Y_plus_dirs <- getKallistoDirs("small")
    #    Y_cells <- c(Y_cells, Y_plus_cells)
    #    Y_dirs <- c(Y_dirs, Y_plus_dirs)

    #    loginfo(sprintf("Considering extended groups: %d additional cells to the main and %d additional cells to the bystanders.", length(X_plus_cells), length(Y_plus_cells)))
    #}

    #chosen_cells <- c(X_cells, Y_cells)
    #chosen_dirs <- c(X_dirs, Y_dirs)

    #chosen_groups <- c( rep("C", length(X_cells)), rep("N", length(Y_cells)) )
    #chosen_groups <- c( rep(first_classes[[1]], length(X_cells)), rep(second_classes[[1]], length(Y_cells)) )
    groups <- unlist(sapply(1:length(cells), function(i) rep(json.config[[i]], length(cells[[i]]))))
    #print(dirs)
    #print(chosen_dirs)
    #print(groups)
    #print(chosen_groups)
    cells <- unname(unlist(cells))
    dirs <- unname(unlist(dirs))
    #print("Cells")
    #print(cells)
    #print("ChosenCells")
    #print(chosen_cells)
    #print("dirs")
    #print(dirs)
    #print("Chosendirs")
    #print(chosen_dirs)
    #print("groups")
    #print(groups)
    #print("Chosen_groups")
    #print(chosen_groups)
    if (length(cells) != length(dirs)) {
    #if (length(X_cells) != length(X_dirs) || length(Y_cells) != length(Y_dirs)) {
        logerror("Different number of cells and dirs.")
        stop()
    }

    if (length(cells) == 0 || length(dirs) == 0) {
    #if (length(X_cells) == 0 || length(Y_cells) == 0) {
        logerror("Empty group(s). Stopping the DE.")
        logerror(json.config)
        #logerror(first_classes[[1]])
        #logerror(second_classes[[1]])
        stop()
    }
    #loginfo(sprintf("Looking at two cell groups: %d clonal and %d non-clonal", length(X_cells), length(Y_cells)))

    return(list("cells" = cells, "dirs" = dirs, "groups" = groups, "group.levels" = make.names(json.config)))
    #return(list("cells" = chosen_cells, "dirs" = chosen_dirs, "groups" = chosen_groups))
}

kallisto2genetable <- function(cells, dirs, countsFromAbundance="no", gene.table.fn="no") {
    loginfo(sprintf("Loading Kallisto quantifications for %d cells...", length(cells)))

    # scater
    #sce <- readKallistoResults(directories = unlist(chosen_dirs), samples = chosen_cells, read_h5=TRUE) 

    # tximport
    #files <- unlist(dirs)
    files <- dirs
    #loginfo("files: ", files)
    names(files) <- cells
    loginfo(sprintf("Reading expression data from %d files", length(names(files))))
    loginfo(sprintf("First filename: %s", files[[1]]))
    #b <- file.exists(files)
    { 
        # Silently read Kallisto.
        sink("/dev/null");  
        tx2symbol <- left_join(grch38_tx2gene, grch38, by="ensgene") %>% select(enstxp, symbol)
        txi <- tximport(files, type = "kallisto", tx2gene = tx2symbol, countsFromAbundance=countsFromAbundance) #tx2gene = grch38_tx2gene)
        sink();
    }

    checkDuplicates(rownames(txi$counts))
    loginfo("\nKallisto quantifications loaded.")

    if (gene.table.fn != "no") {
        loginfo(sprintf("Writing gene table to %s", out.fn))
        write.csv(txi$counts, out.fn)
    }
    return (txi$counts)
}

clonality2sce <- function(clonality.file, sample.name, clear.groups, config.fn) {
    clonality <- clonality2triple(clonality.file, sample.name, clear.groups, config.fn)
    counts <- kallisto2genetable(cells=clonality$cells, dirs=clonality$dirs, countsFromAbundance="no")
    tpm <- kallisto2genetable(cells=clonality$cells, dirs=clonality$dirs, countsFromAbundance="scaledTPM")

    #loginfo(clonality)
    #loginfo(rownames(counts))
    loginfo("kallisto read.")
    #loginfo(counts)

    sce <- SingleCellExperiment(assays = list(counts = counts, tpm = tpm))
    #traceback()
    #loginfo("sce created")
    #print(dim(counts))
    #print(dim(tpm))
    #print(length(colnames(sce)))
    #print(length(clonality$groups))
    colnames(sce) <- mapply(function(cell, group) paste0(cell, "_", group),
                            colnames(sce), clonality$groups)
    colData(sce)$groups <- factor(clonality$groups, levels=clonality$group.levels)
    #print(colData(sce))
    #sce$group.levels <- clonality$group.levels
    return(sce)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    warning("de.R: Not exact number of arguments. Expected 6 args: sample.name(01_SE), clonality.file(json), method(ROTS, SCDE), out_dir(all all kinds of pictures), outfile(for a DE table), clear.groups(1 or 0), config.fn(json)")
} else {
    sample.name         <- args[1]
    clonality.file      <- args[2]
    countsFromAbundance <- args[3]
    gene.table.fn       <- args[4]
    clear.groups        <- args[5]
    config.fn           <- args[6]

    clonality <- clonality2triple(clonality.file, sample.name, clear.groups, config.fn)
    ma3x <- kallisto2genetable(cells=clonality$cells, dirs=clonality$dirs, countsFromAbundance=countsFromAbundance, gene.table.fn=gene.tablefn)
}

