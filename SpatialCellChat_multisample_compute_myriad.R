#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialCellChat)
  library(Matrix)
  library(Seurat)
})

options(
  stringsAsFactors = FALSE,
  future.globals.maxSize = 16 * 1024^3
)
future::plan("sequential")

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
input_rds <- Sys.getenv(
  "INPUT_RDS",
  unset = file.path(project_dir, "Data", "Spatial_Xenium_six_donors_Myometrium.rds")
)
output_dir <- Sys.getenv(
  "OUTPUT_DATA_DIR",
  unset = file.path(project_dir, "Output", "Data", "multi_sample_cellchat_scale5")
)
condition_id <- tolower(Sys.getenv("CONDITION_ID", unset = ""))
signal_type <- tolower(Sys.getenv("SIGNALING_TYPE", unset = ""))
pathway_thresh <- as.numeric(Sys.getenv("PATHWAY_THRESH", unset = "0.05"))

condition_donors <- list(
  control = c("FVB", "FVQ"),
  disease = c("FCM", "FEP", "FVS", "FAM")
)
stopifnot(file.exists(input_rds))
stopifnot(condition_id %in% names(condition_donors))
stopifnot(signal_type %in% c("contact", "secreted"))
stopifnot(is.finite(pathway_thresh))

donors <- condition_donors[[condition_id]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_db <- function(signal_type) {
  custom_lr <- data.frame(
    ligand = "CLEC11A",
    receptor = "KIT",
    pathway_name = "CLEC11A",
    interaction_name = "CLEC11A_KIT",
    interaction_name_2 = "CLEC11A - KIT",
    annotation = "Secreted Signaling",
    evidence = paste(
      "NicheNet-predicted EVT-to-artery interaction;",
      "Greenbaum et al., Nature 2023, Fig. 6h;",
      "doi:10.1038/s41586-023-06298-9"
    ),
    stringsAsFactors = FALSE
  )
  custom_db <- SpatialCellChat::updateCellChatDB(
    db = custom_lr,
    merged = TRUE,
    species_target = "human"
  )

  switch(
    signal_type,
    contact = SpatialCellChat::subsetDB(
      CellChatDB.human,
      search = "Cell-Cell Contact",
      key = "annotation"
    ),
    secreted = SpatialCellChat::subsetDB(
      custom_db,
      search = "Secreted Signaling",
      key = "annotation"
    )
  )
}

output_path <- function() {
  suffix <- switch(
    signal_type,
    contact = "contactRange20",
    secreted = "secreted_range250"
  )
  file.path(
    output_dir,
    paste0("cellchat_", condition_id, "_multi_", suffix, "_final.rds")
  )
}

save_rds <- function(object, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, tmp)
  stopifnot(file.exists(tmp), file.info(tmp)$size > 0)
  file.rename(tmp, path)
}

add_group_centrality <- function(chat) {
  slot_data <- methods::slot(chat, "netP")
  prob <- slot_data$prob
  pval <- slot_data$pval

  if (is.null(prob) || length(dim(prob)) != 3L || dim(prob)[3] == 0L) {
    return(chat)
  }
  if (!is.null(pval) && identical(dim(pval), dim(prob))) {
    prob[pval >= pathway_thresh] <- 0
  }

  centrality_names <- c(
    "outdeg_unweighted",
    "indeg_unweighted",
    "outdeg",
    "indeg",
    "page_rank"
  )
  centrality <- lapply(
    seq_len(dim(prob)[3]),
    function(i) {
      SpatialCellChat:::computeCentralityLocal(
        prob[, , i, drop = TRUE],
        degree.only = TRUE
      )
    }
  )

  slot_data$centr <- array(
    unlist(centrality, use.names = FALSE),
    dim = c(length(centrality_names), dim(prob)[1], dim(prob)[3]),
    dimnames = list(centrality_names, dimnames(prob)[[1]], dimnames(prob)[[3]])
  )
  methods::slot(chat, "netP") <- slot_data
  chat
}

obj <- readRDS(input_rds)
required_cols <- c("donor", "CellTypeManual.l2", "centroid_X_1", "centroid_X_2")
stopifnot(all(required_cols %in% colnames(obj@meta.data)))
stopifnot("Xenium" %in% Seurat::Assays(obj))
stopifnot(all(donors %in% unique(as.character(obj@meta.data$donor))))

cells <- rownames(obj@meta.data)[obj@meta.data$donor %in% donors]
obj <- subset(obj, cells = cells)

DefaultAssay(obj) <- "Xenium"
obj <- Seurat::NormalizeData(
  obj,
  assay = "Xenium",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
obj <- Seurat::JoinLayers(obj, assay = "Xenium")
obj$samples <- factor(as.character(obj$donor), levels = donors)
obj$labels <- factor(obj$CellTypeManual.l2)

expr <- SeuratObject::LayerData(obj, assay = "Xenium", layer = "data")
meta <- obj@meta.data[colnames(expr), , drop = FALSE]
coords <- meta[colnames(expr), c("centroid_X_1", "centroid_X_2")]
colnames(coords) <- c("x", "y")

chat <- SpatialCellChat::createSpatialCellChat(
  object = expr,
  meta = meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = coords,
  spatial.factors = data.frame(
    ratio = rep(1, length(donors)),
    tol = rep(5, length(donors)),
    row.names = donors
  )
)

chat@DB <- make_db(signal_type)
chat <- SpatialCellChat::subsetData(chat)

expressed_cells <- Matrix::rowSums(chat@data.signaling > 0)
keep_genes <- names(expressed_cells[expressed_cells >= 10])
stopifnot(length(keep_genes) > 0)

chat <- SpatialCellChat::preProcessing(chat)
chat@var.features$features <- keep_genes
chat@var.features$features.info <- data.frame(
  features = names(expressed_cells),
  n_cells_expressed = as.numeric(expressed_cells),
  selected = expressed_cells >= 10,
  row.names = names(expressed_cells),
  stringsAsFactors = FALSE
)
chat <- SpatialCellChat::identifyOverExpressedInteractions(
  chat,
  variable.both = FALSE
)

chat <- switch(
  signal_type,
  contact = SpatialCellChat::computeCommunProb(
    chat,
    distance.use = TRUE,
    scale.distance = 5.0,
    contact.dependent = TRUE,
    interaction.range = 250,
    contact.range = 20
  ),
  secreted = SpatialCellChat::computeCommunProb(
    chat,
    distance.use = TRUE,
    scale.distance = 5.0,
    contact.dependent = FALSE,
    interaction.range = 250
  )
)

chat <- SpatialCellChat::filterCommunication(chat, min.cells = 10)
chat <- SpatialCellChat::computeCommunProbPathway(
  chat,
  thresh = pathway_thresh,
  do.group = TRUE,
  do.cell = TRUE
)
chat <- add_group_centrality(chat)

stopifnot("pathways.cell" %in% names(chat@netP))
stopifnot("prob.cell" %in% names(chat@netP))
stopifnot("tmp" %in% names(chat@netP))
stopifnot("prob.cell" %in% names(chat@netP$tmp))
stopifnot("centr" %in% names(chat@netP))

save_rds(chat, output_path())
