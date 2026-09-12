#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialCellChat)
  library(Matrix)
  library(Seurat)
})

options(stringsAsFactors = FALSE)
future::plan("sequential")

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
input_rds <- Sys.getenv(
  "INPUT_RDS",
  unset = file.path(project_dir, "Data", "Spatial_Xenium_six_donors_Myometrium.rds")
)
output_dir <- Sys.getenv(
  "OUTPUT_DATA_DIR",
  unset = file.path(project_dir, "Output", "Data", "single_sample_cellchat")
)
donor_id <- Sys.getenv("DONOR_ID", unset = "")
pathway_thresh <- as.numeric(Sys.getenv("PATHWAY_THRESH", unset = "0.05"))

valid_donors <- c("FCM", "FEP", "FVQ", "FVS", "FVB", "FAM")
stopifnot(file.exists(input_rds))
stopifnot(donor_id %in% valid_donors)
stopifnot(is.finite(pathway_thresh))

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
      non_protein = FALSE
    ),
    secreted = SpatialCellChat::subsetDB(
      custom_db,
      search = "Secreted Signaling",
      non_protein = FALSE
    )
  )
}

output_path <- function(signal_type) {
  suffix <- switch(
    signal_type,
    contact = "contactRange20",
    secreted = "secreted_range250"
  )
  file.path(output_dir, paste0("chat_", donor_id, "_", suffix, "_final.rds"))
}

save_rds <- function(object, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, tmp)
  stopifnot(file.exists(tmp), file.info(tmp)$size > 0)
  file.rename(tmp, path)
}

run_one_signal <- function(seurat_obj, signal_type) {
  DefaultAssay(seurat_obj) <- "Xenium"
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    assay = "Xenium",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  expr <- SeuratObject::LayerData(seurat_obj, assay = "Xenium", layer = "data")
  meta <- seurat_obj@meta.data[colnames(expr), , drop = FALSE]
  coords <- meta[colnames(expr), c("centroid_X_1", "centroid_X_2")]
  colnames(coords) <- c("x", "y")

  chat <- SpatialCellChat::createSpatialCellChat(
    object = expr,
    meta = meta,
    group.by = "CellTypeManual.l2",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = list(ratio = 1, tol = 5)
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
      scale.distance = 0.8,
      contact.dependent = TRUE,
      interaction.range = 250,
      contact.range = 20
    ),
    secreted = SpatialCellChat::computeCommunProb(
      chat,
      distance.use = TRUE,
      scale.distance = 0.8,
      contact.dependent = FALSE,
      interaction.range = 250
    )
  )

  chat <- SpatialCellChat::filterProbability(chat)
  chat <- SpatialCellChat::filterCommunication(
    chat,
    min.cells = NULL,
    min.links = 10,
    min.cells.sr = 10
  )
  chat <- SpatialCellChat::computeAvgCommunProb(
    chat,
    avg.type = "avg",
    min.percent = 0.1,
    min.cells.sr = 10,
    do.permutation = TRUE,
    nboot = 100,
    seed.use = 1L
  )
  chat <- SpatialCellChat::filterCommunication(
    chat,
    min.cells = 10,
    min.links = NULL,
    min.cells.sr = NULL
  )
  chat <- SpatialCellChat::computeCommunProbPathway(
    chat,
    thresh = pathway_thresh,
    do.group = TRUE,
    do.cell = TRUE
  )
  chat <- SpatialCellChat::aggregateNet(chat)

  stopifnot("pathways.cell" %in% names(chat@netP))
  stopifnot("prob.cell" %in% names(chat@netP))
  stopifnot("tmp" %in% names(chat@netP))
  stopifnot("prob.cell" %in% names(chat@netP$tmp))

  save_rds(chat, output_path(signal_type))
}

obj <- readRDS(input_rds)
required_cols <- c("donor", "CellTypeManual.l2", "centroid_X_1", "centroid_X_2")
stopifnot(all(required_cols %in% colnames(obj@meta.data)))
stopifnot("Xenium" %in% Seurat::Assays(obj))

cells <- rownames(obj@meta.data)[obj@meta.data$donor == donor_id]
stopifnot(length(cells) > 0)
obj <- subset(obj, cells = cells)

for (signal_type in c("contact", "secreted")) {
  run_one_signal(obj, signal_type)
  invisible(gc())
}
