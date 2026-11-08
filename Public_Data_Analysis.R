# Original R scripts: verbatim compilation
# Only section boundary comments were added; source-script bytes are unchanged.
# Read README.md before use. This file is not a one-command pipeline.
# Original paths, source() calls, historical alternatives and environment resets remain.


# ==============================================================================
# BEGIN ORIGINAL SCRIPT P01: 00_pub_theme_colors.R
# ==============================================================================
############################################################
## 00_pub_theme_colors.R
## Unified publication-style colors, paths, and helper functions
## For single-cell re-plotting: annotation / proportion / CellChat / CytoTRACE2 / hdWGCNA
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
})

############################################################
## 1. Project paths
## 建议：只需要修改 PROJECT_DIR 即可。
## 如果你在服务器上运行，通常可设为：/home/data/t210636/New_Analysis
############################################################

PROJECT_DIR <- Sys.getenv("SC_PROJECT_DIR", unset = "")

if (!nzchar(PROJECT_DIR)) {
  candidate_project_dirs <- c(
    "/home/data/t210636/New_Analysis",
    "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL",
    getwd()
  )
  PROJECT_DIR <- candidate_project_dirs[dir.exists(candidate_project_dirs)][1]
}

if (is.na(PROJECT_DIR) || !dir.exists(PROJECT_DIR)) {
  PROJECT_DIR <- getwd()
  warning("未自动识别项目总目录，暂时使用当前工作目录：", PROJECT_DIR,
          "\n如路径不对，请手动修改 00_pub_theme_colors.R 中的 PROJECT_DIR。")
}

PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = FALSE)

first_existing_dir <- function(paths, default = paths[1]) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[dir.exists(paths)][1]
  if (is.na(hit)) default else normalizePath(hit, winslash = "/", mustWork = FALSE)
}

first_existing_file <- function(paths, default = paths[1]) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) default else normalizePath(hit, winslash = "/", mustWork = FALSE)
}

ANNO_DIR <- first_existing_dir(c(
  file.path(PROJECT_DIR, "2.Sc_annotation", "new"),
  file.path(PROJECT_DIR, "2.Sc_annotation"),
  file.path(PROJECT_DIR, "2.Sc_annotation_new")
))

RATIO_DIR <- first_existing_dir(c(
  file.path(PROJECT_DIR, "3.Sc_cell_ratio", "new_pub"),
  file.path(PROJECT_DIR, "3.Sc_cell_ratio"),
  file.path(PROJECT_DIR, "3.SC_cell_ratio")
))

COMMUNI_DIR <- first_existing_dir(c(
  file.path(PROJECT_DIR, "4.Sc_cell_communication", "final"),
  file.path(PROJECT_DIR, "4.Sc_cell_communication")
))

TRACE_DIR <- first_existing_dir(c(
  file.path(PROJECT_DIR, "12.CytoTRACE2", "0.new"),
  file.path(PROJECT_DIR, "12.CytoTRACE2"),
  file.path(PROJECT_DIR, "16.CytoTRACE2")
))

# HDWGCNA_DIR <- first_existing_dir(c(
#   file.path(PROJECT_DIR, "6.hdWGCNA", "0.HSC-like_30hub_new"),
#   file.path(PROJECT_DIR, "6.hdWGCNA"),
#   file.path(PROJECT_DIR, "6.SChdWGCNA")
# ))
# 
# ENRICH_DIR <- first_existing_dir(c(
#   file.path(PROJECT_DIR, "7.enrich_analysis", "new"),
#   file.path(PROJECT_DIR, "7.enrich_analysis")
# ))

SEURAT_ANNO_FILE <- first_existing_file(c(
  file.path(ANNO_DIR, "2.resolution_0.5_annotation.qs"),
  file.path(ANNO_DIR, "2.resolution_0.5_annotated.qs"),
  file.path(PROJECT_DIR, "2.Sc_annotation", "new", "2.resolution_0.5_annotation.qs"),
  file.path(PROJECT_DIR, "2.Sc_annotation", "2.resolution_0.5_annotation.qs"),
  file.path(PROJECT_DIR, "2.Sc_annotation", "2.resolution_0.5_annotated.qs"),
  file.path(PROJECT_DIR, "2.Sc_annotation_new", "2.resolution_0.5_annotation.qs")
))

############################################################
## 2. Cell-type order and name mapping
## 你的最终 21 种细胞类型，以短名称为统一展示名称。
############################################################

celltype_levels_short <- c(
  # Stem / Progenitor
  "HSC-like",
  "GranuleProg",
  "Promyelo",
  "MyeloidProg",
  "Progenitor",
  "CyclingProg",
  "GMP-like",
  "LymphoidProg",

  # Monocyte / Macrophage
  "cMono",
  "intMono",
  "ncMono",
  "Macrophage",

  # DC
  "cDC",
  "moDC",
  "pDC",

  # Lymphoid
  "NaiveT",
  "NK_T",
  "Bcell",
  "Plasma",

  # Other
  "Mast",
  "Erythrocyte"
)

## 原始长名称 -> 统一短名称。
## 如果你的 meta.data$celltype 已经是短名称，也会自动保留。
celltype_map <- c(
  "HSC-like cell"          = "HSC-like",
  "HSC-like"               = "HSC-like",
  "MBP-like cell"          = "GranuleProg",
  "GranuleProg"            = "GranuleProg",
  "ProMyelocyte"           = "Promyelo",
  "Promyelo"               = "Promyelo",
  "MDP-like cell"          = "MyeloidProg",
  "MyeloidProg"            = "MyeloidProg",
  "Progenitor-like cell"   = "Progenitor",
  "Progenitor"             = "Progenitor",
  "Cycling Progenitor"     = "CyclingProg",
  "CyclingProg"            = "CyclingProg",
  "GMP-like cell"          = "GMP-like",
  "GMP-like"               = "GMP-like",
  "LMPP-like cell"         = "LymphoidProg",
  "LymphoidProg"           = "LymphoidProg",
  "cMono"                  = "cMono",
  "intMono"                = "intMono",
  "ncMono"                 = "ncMono",
  "Macrophage"             = "Macrophage",
  "cDC"                    = "cDC",
  "moDC"                   = "moDC",
  "pDC"                    = "pDC",
  "Naive T cell"           = "NaiveT",
  "NaiveT"                 = "NaiveT",
  "NK cell"                = "NK_T",
  "NK_T"                   = "NK_T",
  "B cell"                 = "Bcell",
  "Bcell"                  = "Bcell",
  "Plasma cell"            = "Plasma",
  "Plasma"                 = "Plasma",
  "Mast cell"              = "Mast",
  "Mast"                   = "Mast",
  "Erythrocyte"            = "Erythrocyte"
)

celltype_levels_full <- names(celltype_map)

############################################################
## 3. Cell-type color palette
## 原则：同一大类细胞使用相近色系。
## Stem/Progenitor: blue-cyan family
## Monocyte/Macrophage: green-brown family
## DC: cyan-turquoise family
## Lymphoid: purple-magenta family
## Other: yellow-red family
############################################################

celltype_colors_short <- c(
  ##########################################################
  ## Stem / Progenitor: cool blue-cyan family
  ##########################################################
  "HSC-like"      = "#1F4E79",
  "GranuleProg"  = "#2F6DAE",
  "Promyelo"     = "#4C8FD1",
  "MyeloidProg"  = "#6FAEDB",
  "Progenitor"   = "#3E7CB1",
  "CyclingProg"  = "#7B68A6",
  "GMP-like"     = "#00A6A6",
  "LymphoidProg" = "#66C2A5",

  ##########################################################
  ## Monocyte / Macrophage: green-olive-brown family
  ##########################################################
  "cMono"         = "#1B9E77",
  "intMono"       = "#4DAF4A",
  "ncMono"        = "#A6D854",
  "Macrophage"    = "#8C6D31",

  ##########################################################
  ## Dendritic cells: cyan-turquoise family
  ##########################################################
  "cDC"           = "#17A2B8",
  "moDC"          = "#4DBBD5",
  "pDC"           = "#00A087",

  ##########################################################
  ## Lymphoid cells: purple-magenta family
  ##########################################################
  "NaiveT"        = "#6A4C93",
  "NK_T"          = "#8E63A9",
  "Bcell"         = "#AA4499",
  "Plasma"        = "#E7298A",

  ##########################################################
  ## Other: yellow-red family
  ##########################################################
  "Mast"          = "#E6AB02",
  "Erythrocyte"   = "#C44E52"
)

## 长名称颜色：自动继承对应短名称颜色。
celltype_colors_full <- setNames(
  celltype_colors_short[unname(celltype_map)],
  names(celltype_map)
)

fallback_30_colors <- c(
  "#3B6FB6", "#56B4E9", "#009E73", "#5DBB63", "#A6D854",
  "#E69F00", "#F4A261", "#D55E00", "#CC6677", "#C44E52",
  "#AA4499", "#8E63A9", "#6A4C93", "#7F7F7F", "#4D4D4D",
  "#8C6D31", "#B07AA1", "#17A2B8", "#1B9E77", "#66A61E",
  "#E7298A", "#7570B3", "#A6761D", "#E6AB02", "#B3B3B3",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD"
)

get_named_colors <- function(levels_use,
                             palette_named = celltype_colors_short,
                             fallback = fallback_30_colors) {
  levels_use <- as.character(levels_use)
  out <- palette_named[levels_use]
  miss <- is.na(out)
  if (any(miss)) {
    out[miss] <- rep(fallback, length.out = sum(miss))
  }
  names(out) <- levels_use
  out
}

get_celltype_colors <- function(celltypes, use_short = TRUE) {
  if (use_short) {
    get_named_colors(celltypes, palette_named = celltype_colors_short)
  } else {
    get_named_colors(celltypes, palette_named = celltype_colors_full)
  }
}

add_celltype_short <- function(seurat_obj,
                               celltype_col = "celltype",
                               short_col = "celltype_short") {
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("meta.data lacks column: ", celltype_col)
  }
  ct <- as.character(seurat_obj@meta.data[[celltype_col]])
  short <- unname(celltype_map[ct])
  short[is.na(short)] <- ct[is.na(short)]
  seurat_obj@meta.data[[short_col]] <- factor(short, levels = celltype_levels_short)
  seurat_obj
}

############################################################
## 4. Other palettes
############################################################

group_colors <- c(
  "Healthy" = "#3B6FB6",
  "Normal"  = "#3B6FB6",
  "NORMAL"  = "#3B6FB6",
  "Control" = "#3B6FB6",
  "AML"     = "#C44E52",
  "Tumor"   = "#C44E52"

)

cytotrace_score_cols <- c(
  "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
  "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"
)

potency_category_colors <- c(
  "Multipotent"    = "#E6AB02",
  "Oligopotent"    = "#A6D854",
  "Unipotent"      = "#1B9E77",
  "Differentiated" = "#6A4C93",
  "Unknown"        = "#B3B3B3"
)

roe_colors <- colorRampPalette(c(
  "#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF",
  "#FDDBC7", "#EF8A62", "#B2182B"
))(100)

module_colors_pub <- c(
  "turquoise"     = "#1B9E77",
  "blue"          = "#3B6FB6",
  "brown"         = "#8C6D31",
  "yellow"        = "#E6AB02",
  "green"         = "#66A61E",
  "red"           = "#C44E52",
  "black"         = "#4D4D4D",
  "pink"          = "#E7298A",
  "magenta"       = "#AA4499",
  "purple"        = "#7570B3",
  "greenyellow"   = "#A6D854",
  "tan"           = "#D8B365",
  "salmon"        = "#F4A261",
  "cyan"          = "#17A2B8",
  "midnightblue"  = "#313695",
  "lightcyan"     = "#80B1D3",
  "grey"          = "#B3B3B3"
)

############################################################
## 5. Plot themes and save helpers
############################################################

theme_sc_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

## UMAP/tSNE 主图专用主题：无四方框、无常规坐标轴。
theme_dim_shortaxis <- function(base_size = 12) {
  theme_void(base_size = base_size) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = base_size + 2),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      legend.key = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8)
    )
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_plot_pub <- function(plot,
                          filename_base,
                          width = 6,
                          height = 5,
                          dpi = 600,
                          bg = "white") {
  ggsave(paste0(filename_base, ".pdf"), plot = plot, width = width, height = height, bg = bg, useDingbats = FALSE)
  ggsave(paste0(filename_base, ".png"), plot = plot, width = width, height = height, dpi = dpi, bg = bg)
  invisible(filename_base)
}

save_base_pdf <- function(file, width = 7, height = 6, expr) {
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

safe_run <- function(label, expr) {
  message("[RUN] ", label)
  tryCatch(
    force(expr),
    error = function(e) {
      message("[SKIP] ", label, " | ", conditionMessage(e))
      invisible(NULL)
    }
  )
}

message("Loaded 00_pub_theme_colors.R")
message("PROJECT_DIR: ", PROJECT_DIR)
message("ANNO_DIR: ", ANNO_DIR)
message("SEURAT_ANNO_FILE: ", SEURAT_ANNO_FILE)

# END ORIGINAL SCRIPT P01

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P02: 02_BULK_DIFF.R
# ==============================================================================
# ============================================================
# GSE30029 | Bulk DEG analysis only
# Input:
#   1) GSE30029_series_matrix.xlsx
#   2) GPL6947-13512.xlsx
# Output:
#   1) GSE30029_limma_all.xlsx
#   2) GSE30029_DEG_sig.xlsx
#   3) GSE30029_DEG_up.xlsx
#   4) GSE30029_DEG_down.xlsx
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(limma)
  library(impute)
})

options(stringsAsFactors = FALSE)

# -----------------------------
# 1. User parameters
# -----------------------------
base_dir   <- "D:/0.DataAnalysis/AML_BULK_RNA/小论文/GSE30029_new"   # <<< 改成你的目录
gpl_file   <- file.path(base_dir, "1.DEG-seq", "GPL6947-13512.xlsx")
expr_file  <- file.path(base_dir, "1.DEG-seq", "GSE30029_series_matrix.xlsx")
out_dir    <- file.path(base_dir, "DEG-seq_result")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

logFC_cutoff <- 1
p_cutoff     <- 0.05
use_adjP     <- FALSE   # 按你现在要求，默认用 P.Value < 0.05

# -----------------------------
# 2. Sample grouping
#    前31个 Healthy，后46个 AML
# -----------------------------
healthy_samples <- c(
  "GSM743074","GSM743080","GSM743091","GSM743093","GSM743095","GSM743096",
  "GSM743100","GSM743110","GSM743113","GSM743119","GSM743121","GSM743123",
  "GSM743124","GSM743129","GSM743144","GSM743145","GSM743147","GSM743148",
  "GSM743149","GSM743150","GSM743151","GSM743153","GSM743154","GSM743155",
  "GSM743156","GSM743157","GSM743167","GSM743177","GSM743180","GSM743181",
  "GSM743190"
)

aml_samples <- c(
  "GSM743073","GSM743075","GSM743081","GSM743083","GSM743085","GSM743087",
  "GSM743089","GSM743090","GSM743097","GSM743101","GSM743102","GSM743105",
  "GSM743106","GSM743107","GSM743108","GSM743112","GSM743117","GSM743120",
  "GSM743122","GSM743125","GSM743126","GSM743128","GSM743130","GSM743131",
  "GSM743135","GSM743138","GSM743139","GSM743141","GSM743146","GSM743158",
  "GSM743159","GSM743161","GSM743164","GSM743169","GSM743170","GSM743171",
  "GSM743173","GSM743175","GSM743178","GSM743179","GSM743182","GSM743183",
  "GSM743184","GSM743185","GSM743186","GSM743189"
)

sample_order_expected <- c(healthy_samples, aml_samples)

# -----------------------------
# 3. Read GPL annotation
# -----------------------------
ann <- read.xlsx(gpl_file)
colnames(ann) <- trimws(colnames(ann))

id_col <- grep("^ID$|probe|id$", colnames(ann), ignore.case = TRUE, value = TRUE)[1]
symbol_col <- grep("symbol", colnames(ann), ignore.case = TRUE, value = TRUE)[1]

if (is.na(id_col) || is.na(symbol_col)) {
  stop(
    paste0(
      "GPL注释文件列名无法识别。\n当前列名：",
      paste(colnames(ann), collapse = " | ")
    )
  )
}

ann <- ann[, c(id_col, symbol_col)]
colnames(ann) <- c("ID", "gene_symbol")

# -----------------------------
# 4. Read expression matrix
# -----------------------------
expr_raw <- read.xlsx(expr_file)
colnames(expr_raw) <- trimws(colnames(expr_raw))

idref_col <- grep("^ID_REF$|^ID$|probe|id_ref", colnames(expr_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(idref_col)) {
  stop(
    paste0(
      "表达矩阵中未识别到探针ID列。\n当前列名：",
      paste(colnames(expr_raw), collapse = " | ")
    )
  )
}

colnames(expr_raw)[colnames(expr_raw) == idref_col] <- "ID_REF"

# 检查样本列
expr_sample_cols <- setdiff(colnames(expr_raw), "ID_REF")

# 确认前31后46顺序与你给的一致
if (!identical(expr_sample_cols, sample_order_expected)) {
  cat("表达矩阵中的样本顺序与预期不完全一致，开始按你提供的顺序重排...\n")
  
  missing_samples <- setdiff(sample_order_expected, expr_sample_cols)
  extra_samples   <- setdiff(expr_sample_cols, sample_order_expected)
  
  if (length(missing_samples) > 0) {
    stop(paste0("表达矩阵缺少这些样本：", paste(missing_samples, collapse = ", ")))
  }
  if (length(extra_samples) > 0) {
    cat("表达矩阵存在额外样本，将忽略：", paste(extra_samples, collapse = ", "), "\n")
  }
  
  expr_raw <- expr_raw[, c("ID_REF", sample_order_expected), drop = FALSE]
}

# -----------------------------
# 5. Merge annotation + expression
# -----------------------------
geneMatrix <- merge(
  ann, expr_raw,
  by.x = "ID",
  by.y = "ID_REF"
)

geneMatrix <- geneMatrix[, c("gene_symbol", sample_order_expected), drop = FALSE]
geneMatrix$gene_symbol <- as.character(geneMatrix$gene_symbol)
geneMatrix$gene_symbol <- trimws(geneMatrix$gene_symbol)
geneMatrix$gene_symbol <- sub("///.*", "", geneMatrix$gene_symbol)

geneMatrix <- geneMatrix[!is.na(geneMatrix$gene_symbol) & geneMatrix$gene_symbol != "" & geneMatrix$gene_symbol != "---", ]

# 所有样本列转数值
sample_cols <- setdiff(colnames(geneMatrix), "gene_symbol")
geneMatrix[sample_cols] <- lapply(geneMatrix[sample_cols], function(x) suppressWarnings(as.numeric(as.character(x))))

# 按 gene_symbol 合并重复基因，取均值
geneMatrix <- geneMatrix %>%
  group_by(gene_symbol) %>%
  summarise(
    across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# -----------------------------
# 6. Build expression matrix
# -----------------------------
expr_mat <- as.data.frame(geneMatrix, stringsAsFactors = FALSE)
rownames(expr_mat) <- expr_mat$gene_symbol
expr_mat$gene_symbol <- NULL

expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"

# 去掉全NA或均值<=0的基因
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) < ncol(expr_mat), , drop = FALSE]
expr_mat <- expr_mat[rowMeans(expr_mat, na.rm = TRUE) > 0, , drop = FALSE]

# 缺失值填补
expr_mat <- impute.knn(expr_mat)$data

# normalize
expr_mat <- normalizeBetweenArrays(expr_mat)

# 自动判断是否需要 log2
qx <- as.numeric(quantile(expr_mat, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
need_log2 <- (qx[5] > 100) || ((qx[6] - qx[1] > 50) && qx[2] > 0)

if (need_log2) {
  expr_mat[expr_mat <= 0] <- NA
  expr_mat <- log2(expr_mat)
}

# 去掉缺失
expr_mat <- expr_mat[complete.cases(expr_mat), , drop = FALSE]

# -----------------------------
# 7. Differential analysis
# -----------------------------
group <- factor(
  c(rep("Healthy", length(healthy_samples)), rep("AML", length(aml_samples))),
  levels = c("Healthy", "AML")
)

stopifnot(ncol(expr_mat) == length(group))

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Healthy", "AML")

fit <- lmFit(expr_mat, design)
cont.matrix <- makeContrasts(AML - Healthy, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P")
allDiff$gene_symbol <- rownames(allDiff)
allDiff <- allDiff[, c("gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# 防止0值导致后续问题
allDiff$P.Value[allDiff$P.Value == 0] <- 1e-300
allDiff$adj.P.Val[allDiff$adj.P.Val == 0] <- 1e-300

# 阈值列
p_col <- if (use_adjP) "adj.P.Val" else "P.Value"

allDiff$threshold <- "non"
allDiff$threshold[allDiff$logFC >  logFC_cutoff & allDiff[[p_col]] < p_cutoff] <- "up"
allDiff$threshold[allDiff$logFC < -logFC_cutoff & allDiff[[p_col]] < p_cutoff] <- "down"

# 导出差异表
diffSig  <- allDiff[allDiff$threshold != "non", ]
diffUp   <- allDiff[allDiff$threshold == "up", ]
diffDown <- allDiff[allDiff$threshold == "down", ]

# -----------------------------
# 8. Save to xlsx
# -----------------------------
write.xlsx(allDiff,  file = file.path(out_dir, "GSE30029_limma_all.xlsx"), overwrite = TRUE)
write.xlsx(diffSig,  file = file.path(out_dir, "GSE30029_DEG_sig.xlsx"), overwrite = TRUE)
write.xlsx(diffUp,   file = file.path(out_dir, "GSE30029_DEG_up.xlsx"), overwrite = TRUE)
write.xlsx(diffDown, file = file.path(out_dir, "GSE30029_DEG_down.xlsx"), overwrite = TRUE)

cat("Done!\n")
cat("All genes   :", nrow(allDiff), "\n")
cat("Significant :", nrow(diffSig), "\n")
cat("Up          :", nrow(diffUp), "\n")
cat("Down        :", nrow(diffDown), "\n")
cat("Results saved in:\n", out_dir, "\n")



# ============================================================
# Volcano plot
# Directly run after differential analysis (object: allDiff)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------
# 1. Parameters
# -----------------------------
volcano_out_pdf <- file.path(out_dir, "GSE30029_volcano_template_style.pdf")
volcano_out_png <- file.path(out_dir, "GSE30029_volcano_template_style.png")

fc_cut <- logFC_cutoff   # 承接前面差异分析参数
p_cut  <- p_cutoff       # 承接前面差异分析参数

# 这里默认与你前面的差异分析保持一致
# 如果你前面 use_adjP <- FALSE，则这里就用 P.Value
# 如果你前面 use_adjP <- TRUE，则这里就用 adj.P.Val
volcano_p_col <- if (use_adjP) "adj.P.Val" else "P.Value"

main_title  <- "Differentially Expressed Genes in GSE30029"
left_label  <- "Higher in Healthy"
right_label <- "Higher in AML"

# 颜色（按你的要求）
col_down_fill <- "#6A9BCB"   # 蓝
col_up_fill   <- "#CB5148"   # 红
col_ns_fill   <- "#B3B3B3"   # 灰

# 边框色略深一点，避免点发虚
col_down_edge <- "#4F7EAA"
col_up_edge   <- "#A13E37"
col_ns_edge   <- "#8E8E8E"

bg_col <- "white"

# -----------------------------
# 2. Build plotting data
# -----------------------------
vol_df <- allDiff
vol_df$Log2FC <- as.numeric(vol_df$logFC)
vol_df$pvalue <- as.numeric(vol_df[[volcano_p_col]])

vol_df <- vol_df[is.finite(vol_df$Log2FC) & is.finite(vol_df$pvalue), , drop = FALSE]
vol_df$pvalue[vol_df$pvalue == 0] <- 1e-300
vol_df$neglogp <- -log10(vol_df$pvalue)

# 分类
vol_df$regulated <- "Normal"
vol_df$regulated[vol_df$Log2FC >  fc_cut & vol_df$pvalue < p_cut] <- "Up"
vol_df$regulated[vol_df$Log2FC < -fc_cut & vol_df$pvalue < p_cut] <- "Down"

# 百分比
total_n  <- nrow(vol_df)
up_n     <- sum(vol_df$regulated == "Up")
down_n   <- sum(vol_df$regulated == "Down")
up_pct   <- round(up_n / total_n * 100, 1)
down_pct <- round(down_n / total_n * 100, 1)

# -----------------------------
# 3. Axis range
# -----------------------------
x_lim_auto <- ceiling(max(abs(vol_df$Log2FC), na.rm = TRUE))
x_lim <- max(6, x_lim_auto)
x_min <- -x_lim
x_max <-  x_lim

y_lim_auto <- ceiling(max(vol_df$neglogp, na.rm = TRUE))
y_max <- max(20, y_lim_auto)

# 为了更接近模板，y轴按5递增
y_break <- 5

# -----------------------------
# 4. Draw
# -----------------------------
pdf(volcano_out_pdf, width = 9.8, height = 6.7, useDingbats = FALSE)

par(
  bg  = bg_col,
  bty = "n",
  las = 1,
  tcl = -0.35,
  mar = c(5.2, 7.2, 3.8, 1.5),   # 👈 再稍微加一点更稳
  mgp = c(2.6, 0.9, 0)
)

# 空白坐标系
plot(
  vol_df$Log2FC, vol_df$neglogp,
  type = "n",
  xlim = c(x_min, x_max),
  ylim = c(0, y_max),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)

# 坐标轴
axis(1, at = seq(x_min, x_max, by = 1), lwd = 2.8, cex.axis = 1.35)
axis(2, at = seq(0, y_max, by = y_break), lwd = 2.8, cex.axis = 1.35)

# 自定义轴标题
mtext(expression(Log[2]~foldchange), side = 1, line = 3.1, cex = 2.2, font = 1)
mtext(expression(-Log[10]~italic(P)~"-value"),
      side = 2,
      line = 3,
      cex = 2.2,
      las = 0)   # 👈 关键：强制竖直方向

# 阈值线
abline(h = -log10(p_cut), col = "grey65", lwd = 1.6, lty = 3)
abline(v = c(-fc_cut, fc_cut), col = "grey65", lwd = 1.6, lty = 3)

# 中央粗黑线
abline(v = 0, col = "black", lwd = 3.0)

# 中央刻度短横线（模仿模板）
tick_y <- seq(0, y_max, by = y_break)
segments(x0 = -0.14, y0 = tick_y, x1 = 0.14, y1 = tick_y, col = "black", lwd = 2.6)

# 先画灰点
idx_ns <- which(vol_df$regulated == "Normal")
if (length(idx_ns) > 0) {
  points(
    vol_df$Log2FC[idx_ns], vol_df$neglogp[idx_ns],
    pch = 21, cex = 0.58, lwd = 0.35,
    col = col_ns_edge, bg = col_ns_fill
  )
}

# 再画蓝点（Healthy高）
idx_down <- which(vol_df$regulated == "Down")
if (length(idx_down) > 0) {
  points(
    vol_df$Log2FC[idx_down], vol_df$neglogp[idx_down],
    pch = 21, cex = 0.82, lwd = 0.45,
    col = col_down_edge, bg = col_down_fill
  )
}

# 再画红点（AML高）
idx_up <- which(vol_df$regulated == "Up")
if (length(idx_up) > 0) {
  points(
    vol_df$Log2FC[idx_up], vol_df$neglogp[idx_up],
    pch = 21, cex = 0.82, lwd = 0.45,
    col = col_up_edge, bg = col_up_fill
  )
}

# 标题（保留）
title(main = main_title, cex.main = 1.0, font.main = 5)

# -----------------------------
# 5. Top annotation (template style)
# -----------------------------
y_anno <- y_max * 0.92

# 左侧箭头（蓝）
arrows(
  x0 = -0.35, y0 = y_anno,
  x1 = x_min, y1 = y_anno,
  length = 0.09, angle = 20, code = 2,
  lwd = 2.8, col = col_down_fill
)
text(
  x = x_min + 0.6, y = y_anno + y_max * 0.01,
  labels = paste0(down_pct, " %"),
  col = col_down_fill, cex = 1.25, font = 2, pos = 3
)
text(
  x = (x_min + (-0.35)) / 2, y = y_anno + y_max * 0.01,
  labels = left_label,
  col = col_down_fill, cex = 1.15, font = 2, pos = 3
)

# 右侧箭头（红）
arrows(
  x0 = 0.35, y0 = y_anno,
  x1 = x_max, y1 = y_anno,
  length = 0.09, angle = 20, code = 2,
  lwd = 2.8, col = col_up_fill
)
text(
  x = x_max - 0.6, y = y_anno + y_max * 0.01,
  labels = paste0(up_pct, " %"),
  col = col_up_fill, cex = 1.25, font = 2, pos = 3
)
text(
  x = (0.35 + x_max) / 2, y = y_anno + y_max * 0.01,
  labels = right_label,
  col = col_up_fill, cex = 1.15, font = 2, pos = 3
)

dev.off()

# 再输出一份高清 PNG
png(volcano_out_png, width = 2400, height = 1950, res = 300, bg = bg_col)
par(
  bg  = bg_col,
  bty = "n",
  las = 1,
  tcl = -0.35,
  mar = c(5.2, 5.2, 3.8, 1.5),
  mgp = c(2.6, 0.9, 0)
)

plot(
  vol_df$Log2FC, vol_df$neglogp,
  type = "n",
  xlim = c(x_min, x_max),
  ylim = c(0, y_max),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)

axis(1, at = seq(x_min, x_max, by = 1), lwd = 2.8, cex.axis = 1.35)
axis(2, at = seq(0, y_max, by = y_break), lwd = 2.8, cex.axis = 1.35)

mtext(expression(Log[2]~foldchange), side = 1, line = 3.1, cex = 2.2, font = 1)
mtext(expression(-Log[10]~italic(P)~"-value"),
      side = 2,
      line = 3,
      cex = 2.2,
      las = 0)   # 👈 关键：强制竖直方向

abline(h = -log10(p_cut), col = "grey65", lwd = 1.6, lty = 3)
abline(v = c(-fc_cut, fc_cut), col = "grey65", lwd = 1.6, lty = 3)
abline(v = 0, col = "black", lwd = 3.0)

tick_y <- seq(0, y_max, by = y_break)
segments(x0 = -0.14, y0 = tick_y, x1 = 0.14, y1 = tick_y, col = "black", lwd = 2.6)

if (length(idx_ns) > 0) {
  points(vol_df$Log2FC[idx_ns], vol_df$neglogp[idx_ns],
         pch = 21, cex = 0.58, lwd = 0.35, col = col_ns_edge, bg = col_ns_fill)
}
if (length(idx_down) > 0) {
  points(vol_df$Log2FC[idx_down], vol_df$neglogp[idx_down],
         pch = 21, cex = 0.82, lwd = 0.45, col = col_down_edge, bg = col_down_fill)
}
if (length(idx_up) > 0) {
  points(vol_df$Log2FC[idx_up], vol_df$neglogp[idx_up],
         pch = 21, cex = 0.82, lwd = 0.45, col = col_up_edge, bg = col_up_fill)
}

title(main = main_title, cex.main = 1.0, font.main = 2)

arrows(x0 = -0.35, y0 = y_anno, x1 = x_min, y1 = y_anno,
       length = 0.09, angle = 20, code = 2, lwd = 2.8, col = col_down_fill)
text(x = x_min + 0.6, y = y_anno + y_max * 0.01,
     labels = paste0(down_pct, " %"), col = col_down_fill, cex = 1.25, font = 2, pos = 3)
text(x = (x_min + (-0.35)) / 2, y = y_anno + y_max * 0.01,
     labels = left_label, col = col_down_fill, cex = 1.15, font = 2, pos = 3)

arrows(x0 = 0.35, y0 = y_anno, x1 = x_max, y1 = y_anno,
       length = 0.09, angle = 20, code = 2, lwd = 2.8, col = col_up_fill)
text(x = x_max - 0.6, y = y_anno + y_max * 0.01,
     labels = paste0(up_pct, " %"), col = col_up_fill, cex = 1.25, font = 2, pos = 3)
text(x = (0.35 + x_max) / 2, y = y_anno + y_max * 0.01,
     labels = right_label, col = col_up_fill, cex = 1.15, font = 2, pos = 3)

dev.off()

cat("Volcano plot saved to:\n", volcano_out_pdf, "\n", volcano_out_png, "\n")









# END ORIGINAL SCRIPT P02

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P03: 03_BULK_GO_KEGG.R
# ==============================================================================

rm(list = ls())

setwd("D:\\0.DataAnalysis\\AML_BULK_RNA\\小论文\\GSE30029_new\\ENRICH")
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
options(stringsAsFactors = F)
library(ggplot2)
library(openxlsx)

#BiocManager::install("ComplexHeatmap")


kegg.buble = 1
kegg.chord = 1
# setwd("F:\\postgraduate\\PNH_plan\\3.自生成\\10.生成数据的富集分析")
enrich_ID <- read.xlsx("D:\\0.DataAnalysis\\AML_BULK_RNA\\小论文\\GSE30029_new\\ENRICH\\GSE30029_DEG_sig.xlsx")
# enrich_ID <- read.table("3.diff.xls",sep = "\t",header = T)
if(str_detect(colnames(enrich_ID)[3],'logFC')){
  colnames(enrich_ID)[3] <- 'log2FoldChange'
}
colnames(enrich_ID)[1] <- "gene_name"



#GO----
gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

enrich_go <- enrichGO(gene = gene_symbol[,2],
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE)

enrich_go <- mutate(enrich_go, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
GO_ALL <- enrich_go@result
if(median(GO_ALL$richFactor) < 0.1){GO_ALL$richFactor<- GO_ALL$richFactor*10}

GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]
paste0("得到",dim(GO_ALL)[[1]],"个结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# [1]"得到297个结果，其中211个生物学过程，35个细胞组分，51个分子功能"
write.table(GO_ALL,file = "9.GO_ALL.csv",sep = ",",quote = F,row.names = F)
write.table(as.data.frame(BP), '10.GO_BP.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), '11.GO_CC.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), '12.GO_MF.csv', sep = ',', row.names = FALSE, quote = FALSE)

display_number = c(10, 10, 10)
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]
paste0(ego_result_BP$Description,collapse = "（）；")
paste0(ego_result_CC$Description,collapse = "（）；")
paste0(ego_result_MF$Description,collapse = "（）；")
col1 <- c("#E64B35B2", "#DC0000B2",'#E63863','#B2182B', '#D6604D','#EB4B17')
col2 <- c("#4DBBD5B2", "#0073C2B2",'#008280B2','#4393C3', '#23452F','#68A180')
col3 <- c("#E4C755", "#9FA3A8",'#D6E7A3','#F4A371FF','#F9C187FF','#FACF94FF')
Col1 <- sample(col1,1)
Col2 <- sample(col2,1)
Col3 <- sample(col3,1)
#GO柱状图-合并
if(T){
  go_enrich_df <- data.frame(
    ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
    Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
    type=factor(c(rep("BP:biological process", display_number[1]), 
                  rep("CC:cellular component", display_number[2]),
                  rep("MF:molecular function", display_number[3])), 
                levels=c("BP:biological process", "CC:cellular component","MF:molecular function" )))
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
  COLS <- c(Col1, Col2,Col3)
  
  p <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
    geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = COLS) + ###颜色
    coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "GO Terms")+
    theme_bw()+theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=12),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "italic", colour = "black")
    )+scale_x_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "13.GO_barplot.pdf", height = 14, width = 12, p)
} 



# KEGG----
#Sys.setenv(all_proxy="socks5://127.0.0.1:7890")
gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$gene_name,  
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

# #若一直显示连接不上，则跑一下以下代码尝试
# # 1) 清环境变量里的代理（当前 R 会话）
# for (v in c("http_proxy","https_proxy","all_proxy","HTTP_PROXY","HTTPS_PROXY","ALL_PROXY")) Sys.unsetenv(v)
# 
# # 2) 确保 httr/libcurl 也不走代理（当前会话）
# httr::reset_config()                       # 清掉之前 set_config 的设置
# httr::set_config(httr::config(proxy=""))   # 显式禁用代理
# 
# # 3) 重试连通性
# httr::status_code(httr::GET("https://rest.kegg.jp/list/pathway/hsa"))
# 
# options(timeout = 600)

KEGG <- enrichKEGG( gene = gene_symbol$ENTREZID,#基因列表 
                    organism = "hsa",  #物种
                    keyType = "kegg", #
                    minGSSize = 1, 
                    maxGSSize = 500,
                    pvalueCutoff = 1,  
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.9
)
KEGG <- mutate(KEGG, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
kk <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result <- kk@result
# if(median(kegg_result$richFactor) < 0.1){kegg_result$richFactor<- kegg_result$richFactor*10}
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05 & !grepl("Human Diseases", hh$category), ]
dim(hh)#15
hh <- hh[order(hh$Count,decreasing = T),]
write.xlsx(hh, file = "14.KEGG.xlsx")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh[c(1:25),] 
paste0(kk$Description,collapse = "（）；")


if(kegg.buble == 1){
  kk <- hh[c(1:25),] 
  p2 <-  ggplot(kk,aes(y=order,x=Count))+
    geom_point(aes(size=Count,color= pvalue))+# 修改点的大小
    scale_color_gradient(high="#EE0000B2",low = "#008B45B2")+
    labs(color=expression(pvalue,size="Count"), 
         x="Gene Number",y="Pathways",title="KEGG Enrichment")+
    theme_bw()+
    theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=14),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "bold", colour = "black")
    )+scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "15.KEGG_dot.pdf", height = 8, width = 9, p2)
}


# END ORIGINAL SCRIPT P03

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P04: 04_BULK_GSEA.R
# ==============================================================================
rm(list = ls())

############################################################
## GSEA for GSE30029
## Analysis logic unchanged + plotting style unified
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
  library(DOSE)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(gridExtra)
})

options(stringsAsFactors = FALSE)

############################################################
## 0) Global style settings
############################################################
COL_AML     <- "#CB5148"   # positive NES / AML enriched
COL_HEALTHY <- "#6A9BCB"   # negative NES / Healthy enriched
COL_BG      <- "white"
COL_GRID    <- "#D9D9D9"
COL_BORDER  <- "black"
COL_GREY    <- "#B3B3B3"

theme_gsea <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = COL_BG, color = NA),
      plot.background  = element_rect(fill = COL_BG, color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.6),
      panel.grid.minor = element_line(color = "white", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      legend.background = element_rect(fill = COL_BG, color = NA),
      legend.key = element_rect(fill = COL_BG, color = NA)
    )
}

save_pdf <- function(file, plot_obj, width = 8, height = 6) {
  ggsave(file, plot = plot_obj, width = width, height = height,
         units = "in", device = cairo_pdf)
}

############################################################
## 1) Workdir & input
############################################################
setwd("D:/0.DataAnalysis/AML_BULK_RNA/小论文/GSE30029_new/3.GSEA")

data <- read.xlsx("GSE30029_limma_all.xlsx")
colnames(data)[1] <- "gene_symbol"

stopifnot(all(c("gene_symbol", "logFC") %in% colnames(data)))

############################################################
## 2) SYMBOL -> ENTREZID
############################################################
ids_symbol_and_entrezid <- bitr(
  geneID = data$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

data$ENTREZID <- ids_symbol_and_entrezid[
  match(data$gene_symbol, ids_symbol_and_entrezid$SYMBOL),
  "ENTREZID"
]

data <- na.omit(data)

############################################################
## 3) Ranked gene list
## Keep original logic unchanged
############################################################
set.seed(123)
data$logFC <- jitter(data$logFC)
data <- data[order(data$logFC, decreasing = TRUE), ]

gene_list <- data$logFC
names(gene_list) <- data$ENTREZID

############################################################
## 4) KEGG GSEA
## Keep original logic unchanged
############################################################
gsea_result <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 1,
  maxGSSize = 10000,
  pvalueCutoff = 1,
  verbose = FALSE,
  eps = 0
)

gsea_result_df <- as.data.frame(gsea_result)

filtered_results <- gsea_result_df[gsea_result_df$pvalue < 0.05, ]

write.table(gsea_result_df, "16.GSEA.xls", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_results, "17.GSEA_filtered.xls", sep = "\t", quote = FALSE, row.names = FALSE)

gsea_up <- filtered_results[filtered_results$NES > 0, ]
write.xlsx(gsea_up, file = "18.GSEA_UP.xlsx", overwrite = TRUE)

gsea_down <- filtered_results[filtered_results$NES < 0, ]
write.xlsx(gsea_down, file = "19.GSEA_DOWN.xlsx", overwrite = TRUE)

############################################################
## 5) Add helper columns
############################################################
if (nrow(filtered_results) == 0) {
  stop("没有找到 pvalue < 0.05 的 KEGG GSEA 通路。")
}

filtered_results <- filtered_results %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      "positive enrichment (AML-upregulated)",
      "negative enrichment (AML-downregulated)"
    ),
    absNES = abs(NES),
    Description_wrap = str_wrap(Description, width = 35)
  )

############################################################
## 6) Ridgeplot (PDF + TIFF)
############################################################
p_ridge <- ridgeplot(
  gsea_result,
  showCategory = min(20, nrow(gsea_result_df)),
  core_enrichment = TRUE,
  fill = "p.adjust"
) +
  labs(
    title = "KEGG GSEA Ridgeplot",
    x = "Running enrichment distribution"
  ) +
  theme_gsea(13)

save_pdf("15.GSEA_KEGG_ridgeplot.pdf", p_ridge, width = 9, height = 12)

ggsave(
  "15.GSEA_KEGG_ridgeplot.tiff",
  plot = p_ridge,
  width = 9, height = 12, units = "in", dpi = 300, compression = "lzw"
)

############################################################
## 7) Dotplot for significant pathways (fixed)
############################################################
show_n <- min(12, nrow(filtered_results))

# 先做标签换行，避免过长
gsea_result_df$Description <- str_wrap(gsea_result_df$Description, width = 28)
gsea_result@result$Description <- str_wrap(gsea_result@result$Description, width = 28)

p_dot <- dotplot(
  gsea_result,
  showCategory = show_n,
  split = ".sign",
  font.size = 11
) +
  facet_grid(. ~ .sign) +
  labs(title = "KEGG GSEA Dotplot") +
  theme_gsea(12) +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

save_pdf("16.GSEA_KEGG_dotplot.pdf", p_dot, width = 12, height = 8)

ggsave(
  "16.GSEA_KEGG_dotplot.tiff",
  plot = p_dot,
  width = 12, height = 8, units = "in", dpi = 300, compression = "lzw"
)

############################################################
## 8) Top pathways enrichment curves
## Keep original logic similar, but limit to top 5 positive/negative
############################################################
positive_nes <- filtered_results %>%
  filter(NES > 0) %>%
  arrange(desc(NES))

negative_nes <- filtered_results %>%
  filter(NES < 0) %>%
  arrange(NES)

top_positive_nes <- head(positive_nes, 5)
top_negative_nes <- head(negative_nes, 5)

# Positive pathways
if (nrow(top_positive_nes) > 0) {
  for (i in seq_len(nrow(top_positive_nes))) {
    p_pos <- gseaplot2(
      gsea_result,
      geneSetID = top_positive_nes$ID[i],
      color = COL_AML,
      base_size = 12,
      rel_heights = c(1.5, 0.4, 1),
      subplots = 1:3,
      pvalue_table = TRUE,
      ES_geom = "line"
    )
    
    pdf(
      file = paste0("20.GSEA_positive_NES_pathway_", i, ".pdf"),
      width = 10, height = 8, useDingbats = FALSE
    )
    print(p_pos)
    dev.off()
    
    tiff(
      filename = paste0("20.GSEA_positive_NES_pathway_", i, ".tiff"),
      width = 10, height = 8, units = "in", res = 300, compression = "lzw"
    )
    print(p_pos)
    dev.off()
  }
}

# Negative pathways
if (nrow(top_negative_nes) > 0) {
  for (i in seq_len(nrow(top_negative_nes))) {
    p_neg <- gseaplot2(
      gsea_result,
      geneSetID = top_negative_nes$ID[i],
      color = COL_HEALTHY,
      base_size = 12,
      rel_heights = c(1.5, 0.4, 1),
      subplots = 1:3,
      pvalue_table = TRUE,
      ES_geom = "line"
    )
    
    pdf(
      file = paste0("21.GSEA_negative_NES_pathway_", i, ".pdf"),
      width = 10, height = 8, useDingbats = FALSE
    )
    print(p_neg)
    dev.off()
    
    tiff(
      filename = paste0("21.GSEA_negative_NES_pathway_", i, ".tiff"),
      width = 10, height = 8, units = "in", res = 300, compression = "lzw"
    )
    print(p_neg)
    dev.off()
  }
}


############################################################
## 8_1) Multi-pathway GSEA plot in one figure
############################################################

multi_pathway_keywords <- c(
  "DNA replication",
  "Cell cycle",
  "Oxidative phosphorylation",
  "Mismatch repair",
  "Glycolysis/Gluconeogenesis",
  "TNF signaling pathway"
)

multi_paths <- gsea_result_df %>%
  filter(Description %in% multi_pathway_keywords)

if (nrow(multi_paths) == 0) {
  multi_paths <- gsea_result_df %>%
    filter(str_detect(
      Description,
      paste(multi_pathway_keywords, collapse = "|")
    ))
}

multi_paths <- multi_paths %>%
  distinct(ID, .keep_all = TRUE)

print(multi_paths[, c("ID", "Description", "NES", "pvalue", "p.adjust")])

if (nrow(multi_paths) < 1) {
  stop("No matched pathways were found. Please check pathway names in gsea_result_df$Description.")
}

multi_ids <- multi_paths$ID

# 自定义颜色，顺序要和 multi_ids 对应
multi_colors <- c(
  "#F8766D",  # Cytokine
  "#A3A500",  # Neuroactive
  "#00BF7D",  # Olfactory
  "#00B0F6",  # Ribosome
  "#E76BF3",  # Spliceosome
  "#58A79C"
)

multi_colors <- multi_colors[seq_along(multi_ids)]

p_multi_gsea <- gseaplot2(
  gsea_result,
  geneSetID = multi_ids,
  color = multi_colors,
  base_size = 12,
  rel_heights = c(1.5, 0.45, 1),
  subplots = 1:3,
  pvalue_table = FALSE,
  ES_geom = "line"
)

pdf(
  file = "25.GSEA_multi_pathway_enrichment_plot.pdf",
  width = 9,
  height = 8,
  useDingbats = FALSE
)
print(p_multi_gsea)
dev.off()

tiff(
  filename = "25.GSEA_multi_pathway_enrichment_plot.tiff",
  width = 9,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(p_multi_gsea)
dev.off()

write.xlsx(
  multi_paths,
  file = "26.GSEA_multi_pathway_selected_pathways.xlsx",
  overwrite = TRUE
)



############################################################
## 9) Final summary plot
## Requirement: top 25 pathways ranked by |NES|
############################################################
top25_paths <- filtered_results %>%
  arrange(desc(absNES)) %>%
  slice_head(n = min(25, nrow(filtered_results))) %>%
  mutate(
    Direction = factor(Direction, levels = c("negative enrichment (AML-downregulated)", "positive enrichment (AML-upregulated)")),
    Description_wrap = factor(
      Description_wrap,
      levels = rev(unique(Description_wrap))
    )
  )

p_summary <- ggplot(
  top25_paths,
  aes(x = NES, y = Description_wrap)
) +
  geom_segment(
    aes(x = 0, xend = NES, y = Description_wrap, yend = Description_wrap, color = Direction),
    linewidth = 1.3
  ) +
  geom_point(
    aes(size = setSize, fill = Direction),
    shape = 21, color = "black", stroke = 0.4
  ) +
  scale_color_manual(values = c("positive enrichment (AML-upregulated)" = COL_AML, "negative enrichment (AML-downregulated)" = COL_HEALTHY)) +
  scale_fill_manual(values = c("positive enrichment (AML-upregulated)" = COL_AML, "negative enrichment (AML-downregulated)" = COL_HEALTHY)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, color = "black") +
  labs(
    title = "Top25 KEGG pathways ranked by |NES|",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "setSize",
    color = "Direction",
    fill = "Direction"
  ) +
  theme_gsea(13) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right"
  )

save_pdf("22.GSEA_summary_top25_absNES.pdf", p_summary, width = 10, height = 9)

ggsave(
  "22.GSEA_summary_top25_absNES.tiff",
  plot = p_summary,
  width = 10, height = 9, units = "in", dpi = 300, compression = "lzw"
)

############################################################
## 10) Alternative bar-style summary (new, PDF)
############################################################
p_bar <- ggplot(
  top25_paths,
  aes(x = reorder(Description_wrap, NES), y = NES, fill = Direction)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c("positive enrichment (AML-upregulated)" = COL_AML, "negative enrichment (AML-downregulated)" = COL_HEALTHY)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.7, color = "black") +
  labs(
    title = "Top25 KEGG pathways ranked by |NES|",
    x = NULL,
    y = "NES",
    fill = "Direction"
  ) +
  theme_gsea(13) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )

save_pdf("23.GSEA_barplot_top25_absNES.pdf", p_bar, width = 10, height = 9)

############################################################
## 11) Save ranked table used in summary plot
############################################################
write.xlsx(top25_paths, "24.GSEA_top25_absNES.xlsx", overwrite = TRUE)

cat("GSEA analysis finished.\n")
cat("Main outputs:\n")
cat("15.GSEA_KEGG_ridgeplot.pdf / .tiff\n")
cat("16.GSEA_KEGG_dotplot.pdf\n")
cat("20.GSEA_positive_NES_pathway_*.pdf / .tiff\n")
cat("21.GSEA_negative_NES_pathway_*.pdf / .tiff\n")
cat("22.GSEA_summary_top25_absNES.pdf / .tiff\n")
cat("23.GSEA_barplot_top25_absNES.pdf\n")
cat("24.GSEA_top25_absNES.xlsx\n")




















rm(list = ls())

############################################################
## GSEA for GSE30029
## Only 4 single-pathway GSEA plots
## Template style: ES curve + hit index + ranked list metric
## Pathways:
##   DNA replication
##   Cell cycle
##   Glycolysis / Gluconeogenesis
##   TNF signaling pathway
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
  library(DOSE)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(cowplot)
  library(grid)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)

############################################################
## 0) Global parameters
############################################################

## 修改为你的工作目录
WORK_DIR <- "D:/0.DataAnalysis/AML_BULK_RNA/小论文/GSE30029_new/3.GSEA"
INPUT_FILE <- "GSE30029_limma_all.xlsx"

## 输出目录
OUT_DIR <- "GSEA_4_single_pathway_template"

## 右上角标注；如果想完全仿照模板中的 celltype，可改成："celltype: AML CD34+"
ANNOTATION_LABEL <- "comparison: AML vs Healthy"

## 只画这4条通路
SINGLE_PATHWAY_TERMS <- c(
  "DNA replication",
  "Cell cycle",
  "Glycolysis / Gluconeogenesis",
  "TNF signaling pathway"
)

## 是否标注曲线上的top genes。模板图中有基因标注，因此默认 TRUE。
SHOW_TOP_GENES <- TRUE
TOP_GENE_N <- 10

## 颜色设置：尽量贴近你给的模板图
COL_ES_LINE <- "#58B957"    # 绿色ES曲线
COL_UP <- "#C0392B"         # 红色：positively correlated
COL_DOWN <- "#3B4CC0"       # 蓝色：negatively correlated
COL_BG_POS <- "#FADBD8"     # 淡粉背景
COL_BG_NEG <- "#E8EEF9"     # 淡蓝背景；负NES时自动使用
COL_GREY_GRID <- "grey75"

############################################################
## 0.1) Helper functions
############################################################

format_pvalue_template <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("<1e-300")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  return(signif(x, 3))
}

save_plot_highres <- function(plot_obj,
                              filename_prefix,
                              width = 6.5,
                              height = 8,
                              dpi = 600) {
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  ggsave(
    filename = paste0(filename_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = pdf_device
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "png",
    bg = "white"
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
}

clean_file_name <- function(x) {
  x %>%
    str_replace_all("[/\\\\:*?\"<>|]", "_") %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_sub(1, 90)
}

normalize_term <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "")
}

resolve_pathway_ids <- function(terms, result_df) {
  resolved_ids <- character(0)
  
  for (one_term in terms) {
    one_term_trim <- str_trim(one_term)
    
    ## 1. KEGG ID精准匹配，例如 hsa04110
    if (one_term_trim %in% result_df$ID) {
      resolved_ids <- c(resolved_ids, one_term_trim)
      next
    }
    
    ## 2. Description精准匹配
    exact_desc <- result_df %>%
      filter(tolower(Description) == tolower(one_term_trim)) %>%
      arrange(pvalue) %>%
      pull(ID)
    
    if (length(exact_desc) > 0) {
      resolved_ids <- c(resolved_ids, exact_desc[1])
      next
    }
    
    ## 3. 去掉空格、斜杠等符号后匹配，适合 Glycolysis / Gluconeogenesis
    norm_input <- normalize_term(one_term_trim)
    norm_desc_match <- result_df %>%
      mutate(.norm_desc = normalize_term(Description)) %>%
      filter(.norm_desc == norm_input) %>%
      arrange(pvalue) %>%
      pull(ID)
    
    if (length(norm_desc_match) > 0) {
      resolved_ids <- c(resolved_ids, norm_desc_match[1])
      next
    }
    
    ## 4. 模糊匹配
    fuzzy_desc <- result_df %>%
      filter(str_detect(tolower(Description), fixed(tolower(one_term_trim)))) %>%
      arrange(pvalue) %>%
      pull(ID)
    
    if (length(fuzzy_desc) > 0) {
      resolved_ids <- c(resolved_ids, fuzzy_desc[1])
      next
    }
    
    stop("No matched KEGG pathway found for: ", one_term_trim)
  }
  
  unique(resolved_ids)
}

############################################################
## 1) Workdir & input
############################################################

setwd(WORK_DIR)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

expr_df <- read.xlsx(INPUT_FILE)
colnames(expr_df)[1] <- "gene_symbol"

stopifnot(all(c("gene_symbol", "logFC") %in% colnames(expr_df)))

############################################################
## 2) SYMBOL -> ENTREZID
############################################################

ids_symbol_and_entrezid <- bitr(
  geneID = expr_df$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

expr_df$ENTREZID <- ids_symbol_and_entrezid[
  match(expr_df$gene_symbol, ids_symbol_and_entrezid$SYMBOL),
  "ENTREZID"
]

expr_df <- expr_df %>%
  filter(!is.na(ENTREZID), !is.na(logFC), !is.na(gene_symbol))

## 去重：一个ENTREZID只保留 |logFC| 最大的一条，避免GSEA排序向量重名
expr_df <- expr_df %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

## 用于把ENTREZ ID转回 gene symbol 进行曲线标注
entrez_to_symbol <- setNames(expr_df$gene_symbol, expr_df$ENTREZID)

############################################################
## 3) Ranked gene list
############################################################

set.seed(123)
expr_df$logFC <- jitter(expr_df$logFC)

expr_df <- expr_df %>%
  arrange(desc(logFC))

gene_list <- expr_df$logFC
names(gene_list) <- expr_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

############################################################
## 4) KEGG GSEA
############################################################

gsea_result <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 1,
  verbose = FALSE,
  eps = 0
)

gsea_result_df <- as.data.frame(gsea_result)

write.table(
  gsea_result_df,
  file.path(OUT_DIR, "01.GSEA_KEGG_all_results.xls"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.xlsx(
  gsea_result_df,
  file = file.path(OUT_DIR, "01.GSEA_KEGG_all_results.xlsx"),
  overwrite = TRUE
)

############################################################
## 5) Single-pathway GSEA plot function
## Template-like style from your uploaded figure
############################################################

plot_single_gsea_like_template <- function(gsea_obj,
                                           geneSetID,
                                           entrez_to_symbol = NULL,
                                           annotation_label = NULL,
                                           show_top_genes = TRUE,
                                           top_gene_n = 10,
                                           title_width = 60,
                                           directional_bg = TRUE) {
  
  result_df <- as.data.frame(gsea_obj)
  one_result <- result_df %>% filter(ID == geneSetID)
  
  if (nrow(one_result) == 0) {
    stop("geneSetID not found in GSEA result: ", geneSetID)
  }
  
  ## enrichplot内部函数，提取running enrichment score和hit position
  gp <- enrichplot:::gsInfo(gsea_obj, geneSetID = geneSetID)
  geneList <- gsea_obj@geneList
  pathway_genes <- gsea_obj@geneSets[[geneSetID]]
  
  gene_rank_df <- data.frame(
    x = seq_along(geneList),
    metric = as.numeric(geneList),
    gene = names(geneList),
    stringsAsFactors = FALSE
  )
  
  hit_data <- gp[gp$position == 1, ]
  
  x_max <- max(gene_rank_df$x)
  n_total <- x_max
  zero_cross <- which.min(abs(gene_rank_df$metric))
  x_pad <- max(120, round(x_max * 0.03))
  x_left_lim <- 1 - x_pad
  x_right_lim <- x_max + x_pad
  
  nes_val <- one_result$NES[1]
  p_val <- one_result$pvalue[1]
  
  pathway_title <- str_wrap(one_result$Description[1], width = title_width)
  
  y_min <- min(gp$runningScore, na.rm = TRUE)
  y_max <- max(gp$runningScore, na.rm = TRUE)
  y_span <- y_max - y_min
  if (y_span == 0) y_span <- 1
  
  ## NES和pvalue注释：正向富集放右上，负向富集放左下，且不再显示comparison
  if (nes_val < 0) {
    anno_df <- data.frame(
      x = c(x_max * 0.06, x_max * 0.06),
      y = c(y_min + y_span * 0.22, y_min + y_span * 0.10),
      label = c(
        paste0("NES: ", sprintf("%.2f", nes_val)),
        paste0("pvalue: ", format_pvalue_template(p_val))
      ),
      fontface = c("bold", "bold.italic"),
      stringsAsFactors = FALSE
    )
    anno_hjust <- 0
  } else {
    anno_df <- data.frame(
      x = c(x_max * 0.94, x_max * 0.94),
      y = c(y_min + y_span * 0.88, y_min + y_span * 0.76),
      label = c(
        paste0("NES: ", sprintf("%.2f", nes_val)),
        paste0("pvalue: ", format_pvalue_template(p_val))
      ),
      fontface = c("bold", "bold.italic"),
      stringsAsFactors = FALSE
    )
    anno_hjust <- 1
  }
  
  ############################################################
  ## 5.1) Prepare gene labels on ES curve
  ############################################################
  
  label_genes_df <- gene_rank_df %>%
    filter(gene %in% pathway_genes)
  
  if (nes_val >= 0) {
    label_genes_df <- label_genes_df %>%
      filter(metric > 0) %>%
      arrange(x) %>%
      slice_head(n = top_gene_n)
  } else {
    label_genes_df <- label_genes_df %>%
      filter(metric < 0) %>%
      arrange(desc(x)) %>%
      slice_head(n = top_gene_n) %>%
      arrange(x)
  }
  
  ## 如果按NES方向找不到，退而求其次用rank metric绝对值最大的通路基因
  if (nrow(label_genes_df) == 0) {
    label_genes_df <- gene_rank_df %>%
      filter(gene %in% pathway_genes) %>%
      arrange(desc(abs(metric))) %>%
      slice_head(n = top_gene_n) %>%
      arrange(x)
  }
  
  if (nrow(label_genes_df) > 0) {
    label_genes_df$es_y <- gp$runningScore[match(label_genes_df$x, gp$x)]
    label_genes_df$label <- label_genes_df$gene
    
    if (!is.null(entrez_to_symbol)) {
      mapped_symbol <- entrez_to_symbol[label_genes_df$gene]
      label_genes_df$label <- ifelse(
        !is.na(mapped_symbol) & mapped_symbol != "",
        mapped_symbol,
        label_genes_df$gene
      )
    }
    
    label_genes_df$nudge_y <- ifelse(
      seq_len(nrow(label_genes_df)) %% 2 == 1,
      ifelse(nes_val >= 0, 0.06 * y_span, -0.06 * y_span),
      ifelse(nes_val >= 0, -0.025 * y_span, 0.025 * y_span)
    )
    label_genes_df$nudge_x <- ifelse(
      seq_len(nrow(label_genes_df)) %% 2 == 1,
      x_max * 0.010,
      -x_max * 0.010
    )
  }
  
  ############################################################
  ## 5.2) Enrichment score curve
  ############################################################
  
  if (directional_bg && nes_val < 0) {
    bg_layer <- annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = 0,
      fill = COL_BG_NEG,
      alpha = 0.65
    )
    es_upper_lim <- max(0.06, y_max * 1.04)
    es_lower_lim <- y_min - 0.05 * y_span
  } else {
    bg_layer <- annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = 0,
      ymax = Inf,
      fill = COL_BG_POS,
      alpha = 0.65
    )
    es_upper_lim <- y_max + 0.05 * y_span
    es_lower_lim <- min(-0.06, y_min * 1.04)
  }
  
  p_es <- ggplot(gp, aes(x = x, y = runningScore)) +
    bg_layer +
    geom_hline(
      yintercept = pretty(c(y_min, y_max), n = 6),
      linetype = "dotted",
      linewidth = 0.35,
      color = COL_GREY_GRID
    ) +
    geom_vline(
      xintercept = pretty(c(0, x_max), n = 8),
      linetype = "dotted",
      linewidth = 0.35,
      color = COL_GREY_GRID
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      linewidth = 0.55,
      color = "grey40"
    ) +
    geom_line(
      color = COL_ES_LINE,
      linewidth = 1.65,
      alpha = 0.98
    ) +
    geom_point(
      color = COL_ES_LINE,
      size = 0.75,
      alpha = 0.65
    )
  
  if (show_top_genes && nrow(label_genes_df) > 0) {
    p_es <- p_es +
      geom_point(
        data = label_genes_df,
        aes(x = x, y = es_y),
        inherit.aes = FALSE,
        size = 1.6,
        shape = 16,
        color = "black"
      ) +
      geom_text_repel(
        data = label_genes_df,
        aes(x = x, y = es_y, label = label),
        inherit.aes = FALSE,
        size = 3.6,
        color = "black",
        fontface = "italic",
        nudge_y = label_genes_df$nudge_y,
        nudge_x = label_genes_df$nudge_x,
        segment.size = 0.35,
        segment.color = "grey45",
        min.segment.length = 0,
        box.padding = 0.18,
        point.padding = 0.08,
        max.overlaps = 100,
        force = 0.8,
        force_pull = 0.8
      )
  }
  
  p_es <- p_es +
    geom_text(
      data = anno_df,
      aes(x = x, y = y, label = label, fontface = fontface),
      inherit.aes = FALSE,
      size = 5.0,
      color = "black",
      hjust = anno_hjust,
      vjust = 1
    ) +
    labs(
      title = pathway_title,
      x = NULL,
      y = "Enrichment Score"
    ) +
    scale_x_continuous(
      limits = c(x_left_lim, x_right_lim),
      breaks = pretty(c(0, x_max), n = 4),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(es_lower_lim, es_upper_lim),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        color = "black",
        size = 17,
        lineheight = 0.95
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 15,
        color = "black",
        margin = margin(r = 10)
      ),
      axis.text.y = element_text(
        size = 10.5,
        color = "black"
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(
        linewidth = 0.90,
        color = "black"
      ),
      plot.margin = margin(8, 2, 2, 2),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.90
      )
    )
  
  ############################################################
  ## 5.3) Hit index and gradient strip
  ############################################################
  
  ## 让红-白-蓝渐变条铺满完整横坐标范围，包括左右留白区
  gradient_n <- 1500
  gradient_x <- seq(x_left_lim, x_right_lim, length.out = gradient_n)
  n_left_grad <- sum(gradient_x <= zero_cross)
  n_right_grad <- gradient_n - n_left_grad
  
  gradient_colors <- c(
    colorRampPalette(c(COL_UP, "white"))(max(n_left_grad, 2)),
    colorRampPalette(c("white", COL_DOWN))(max(n_right_grad, 2))
  )
  gradient_colors <- gradient_colors[seq_len(gradient_n)]
  
  gradient_df <- data.frame(
    x = gradient_x,
    color = gradient_colors,
    stringsAsFactors = FALSE
  )
  
  gradient_tile_width <- (x_right_lim - x_left_lim) / gradient_n * 1.05
  
  p_hit <- ggplot() +
    geom_tile(
      data = gradient_df,
      aes(x = x, y = -0.15, fill = color),
      width = gradient_tile_width,
      height = 0.30
    ) +
    scale_fill_identity() +
    geom_segment(
      data = hit_data,
      aes(x = x, xend = x, y = 1, yend = 0),
      linewidth = 0.55,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(x_left_lim, x_right_lim),
      breaks = pretty(c(0, x_max), n = 4),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(-0.30, 1.10),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.margin = margin(0, 2, 0, 2),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.90
      )
    )
  
  ############################################################
  ## 5.4) Ranked list metric
  ############################################################
  
  max_metric <- 6
  min_metric <- -6
  
  p_rank <- ggplot(gene_rank_df, aes(x = x, y = metric)) +
    geom_hline(
      yintercept = c(-6, -3, 0, 3, 6),
      linetype = "dotted",
      linewidth = 0.35,
      color = COL_GREY_GRID
    ) +
    geom_vline(
      xintercept = pretty(c(0, x_max), n = 8),
      linetype = "dotted",
      linewidth = 0.35,
      color = COL_GREY_GRID
    ) +
    geom_area(
      fill = "grey60",
      color = "black",
      linewidth = 0.45
    ) +
    ## 仅保留：正相关对应左侧上半区、负相关对应右侧下半区
    annotate(
      "rect",
      xmin = x_left_lim,
      xmax = zero_cross,
      ymin = 0,
      ymax = 5,
      fill = COL_UP,
      alpha = 0.16
    ) +
    annotate(
      "rect",
      xmin = zero_cross,
      xmax = x_right_lim,
      ymin = -5,
      ymax = 0,
      fill = COL_DOWN,
      alpha = 0.16
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.55,
      color = "grey30"
    ) +
    geom_vline(
      xintercept = zero_cross,
      linetype = "dashed",
      linewidth = 0.60,
      color = "black"
    ) +
    annotate(
      "text",
      x = zero_cross + x_max * 0.08,
      y = 1.2,
      label = paste0("Zero cross at ", zero_cross),
      size = 3.8,
      color = "black"
    ) +
    annotate(
      "text",
      x = x_max * 0.16,
      y = 3.6,
      label = "Positively correlated",
      size = 4.0,
      color = COL_UP,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = x_max * 0.97,
      y = -3.6,
      label = "Negatively correlated",
      size = 4.0,
      color = COL_DOWN,
      fontface = "bold",
      hjust = 1
    ) +
    scale_x_continuous(
      limits = c(x_left_lim, x_right_lim),
      breaks = pretty(c(0, x_max), n = 4),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(-6, 6),
      breaks = c(-6, -3, 0, 3, 6),
      expand = c(0, 0)
    ) +
    labs(
      x = "Rank in Ordered Dataset",
      y = "Ranked List Metric"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.title.x = element_text(
        face = "bold",
        size = 15,
        color = "black"
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 15,
        color = "black",
        margin = margin(r = 10)
      ),
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.90,
        color = "black"
      ),
      plot.margin = margin(2, 2, 8, 2),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.90
      )
    )
  
  ############################################################
  ## 5.5) Combine
  ############################################################
  
  aligned_plots <- cowplot::align_plots(
    p_es,
    p_hit,
    p_rank,
    align = "v",
    axis = "lr"
  )
  
  p_final <- cowplot::plot_grid(
    plotlist = aligned_plots,
    ncol = 1,
    rel_heights = c(3.30, 0.52, 1.12)
  )
  
  return(p_final)
}

############################################################
## 6) Resolve target pathways and draw plots
############################################################

single_pathway_ids <- resolve_pathway_ids(
  terms = SINGLE_PATHWAY_TERMS,
  result_df = gsea_result_df
)

single_pathway_info <- gsea_result_df %>%
  filter(ID %in% single_pathway_ids) %>%
  arrange(match(ID, single_pathway_ids))

print(single_pathway_info[, c("ID", "Description", "NES", "pvalue", "p.adjust")])

write.xlsx(
  single_pathway_info,
  file = file.path(OUT_DIR, "02.Selected_4_single_pathways.xlsx"),
  overwrite = TRUE
)

## 逐条画图并保存
single_plot_list <- list()

for (i in seq_along(single_pathway_ids)) {
  one_id <- single_pathway_ids[i]
  one_desc <- gsea_result_df$Description[match(one_id, gsea_result_df$ID)]
  
  message("Drawing single GSEA plot: ", one_id, " | ", one_desc)
  
  p_one <- plot_single_gsea_like_template(
    gsea_obj = gsea_result,
    geneSetID = one_id,
    entrez_to_symbol = entrez_to_symbol,
    show_top_genes = SHOW_TOP_GENES,
    top_gene_n = TOP_GENE_N,
    title_width = 60,
    directional_bg = TRUE
  )
  
  single_plot_list[[one_desc]] <- p_one
  
  out_prefix <- file.path(
    OUT_DIR,
    paste0(
      sprintf("03.Single_GSEA_%02d_", i),
      one_id,
      "_",
      clean_file_name(one_desc)
    )
  )
  
  save_plot_highres(
    p_one,
    out_prefix,
    width = 6.8,
    height = 5.6,
    dpi = 600
  )
}

## 可选：把4张拼成2×2总览图，方便检查；如果不需要可以删除这一段
p_combined_4 <- cowplot::plot_grid(
  plotlist = single_plot_list,
  ncol = 2,
  align = "hv"
)

save_plot_highres(
  p_combined_4,
  file.path(OUT_DIR, "04.Combined_4_single_GSEA_template_style"),
  width = 13,
  height = 16,
  dpi = 600
)

############################################################
## 7) Output information
############################################################

cat("GSEA single-pathway plotting finished.\n")
cat("Output directory: ", OUT_DIR, "\n", sep = "")
cat("Main outputs:\n")
cat("01.GSEA_KEGG_all_results.xlsx / .xls\n")
cat("02.Selected_4_single_pathways.xlsx\n")
cat("03.Single_GSEA_01-04_*.pdf / .png / .tiff\n")
cat("04.Combined_4_single_GSEA_template_style.pdf / .png / .tiff\n")

# END ORIGINAL SCRIPT P04

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P05: 05_BULK_WGCNA.R
# ==============================================================================
rm(list = ls())

############################################################
## WGCNA for GSE30029 (AML CD34+ vs HEALTHY CD34+)
## Final integrated version:
## analysis logic unchanged + figure style refined
############################################################

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(WGCNA)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

############################################################
## A. Global style settings
############################################################

COL_HEALTHY <- "#6A9BCB"
COL_AML     <- "#CB5148"
COL_BG      <- "white"
COL_BORDER  <- "black"
COL_GRID    <- "#D9D9D9"
COL_GREY    <- "#A8A8A8"

module_trait_palette <- colorRampPalette(
  c("#2C7FB8", "#F7F7F7", "#D73027")
)(100)

save_pdf_base <- function(file, width = 8, height = 6) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(file, width = width, height = height)
  } else {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  }
}

save_ggplot_all <- function(plot_obj,
                            filename_prefix,
                            width = 8,
                            height = 6,
                            dpi = 600) {
  
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  ggsave(
    filename = paste0(filename_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = pdf_device
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "tiff",
    compression = "lzw"
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "png",
    bg = "white"
  )
}

set_base_par <- function(mar = c(5, 5, 3, 2),
                         bg = COL_BG,
                         cex = 1.05,
                         cex.lab = 1.15,
                         cex.axis = 1.05,
                         cex.main = 1.25) {
  par(
    bg = bg,
    bty = "n",
    las = 1,
    tcl = -0.3,
    mar = mar,
    mgp = c(2.5, 0.8, 0),
    font = 2,
    font.lab = 2,
    font.axis = 2,
    font.main = 2,
    cex = cex,
    cex.lab = cex.lab,
    cex.axis = cex.axis,
    cex.main = cex.main
  )
}

############################################################
## 1) Set workdir
############################################################

setwd("D:/0.DataAnalysis/AML_BULK_RNA/小论文/GSE30029_new/4.WGCNA")

############################################################
## 2) Read expression matrix
############################################################

expr_file <- "GSE30029_logFPKM.xlsx"
expr <- read.xlsx(expr_file)

stopifnot(ncol(expr) >= 3)

############################################################
## 3) Ensure gene column name
############################################################

gene_col <- 1
colnames(expr)[gene_col] <- "gene_symbol"

############################################################
## 4) Clean gene symbols & collapse duplicates
############################################################

sample_cols <- setdiff(colnames(expr), "gene_symbol")

expr[sample_cols] <- lapply(
  expr[sample_cols],
  function(x) as.numeric(as.character(x))
)

expr <- expr %>%
  mutate(gene_symbol = as.character(gene_symbol)) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(
    across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

############################################################
## 5) Transpose -> datExpr0 (samples x genes)
############################################################

expr_mat <- as.data.frame(expr, stringsAsFactors = FALSE)
rownames(expr_mat) <- expr_mat$gene_symbol
expr_mat$gene_symbol <- NULL

datExpr0 <- t(as.matrix(expr_mat))
storage.mode(datExpr0) <- "numeric"

cat("datExpr0 dim (samples x genes):", dim(datExpr0), "\n")

############################################################
## 6) Variance filtering: remove bottom 25%
############################################################

geneVar <- apply(datExpr0, 2, var)
summary(geneVar)

var_cut <- quantile(geneVar, 0.25, na.rm = TRUE)
keep_genes <- geneVar > var_cut
datExpr <- datExpr0[, keep_genes]

cat("Genes retained after variance filtering:", ncol(datExpr), "\n")

############################################################
## 7) Remove bad samples/genes
############################################################

gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

cat("datExpr dim after gsg:", dim(datExpr), "\n")

############################################################
## 8) Build trait data: AML vs Healthy
############################################################

aml_samples <- c(
  "GSM743073","GSM743075","GSM743081","GSM743083","GSM743085","GSM743087",
  "GSM743089","GSM743090","GSM743097","GSM743101","GSM743102","GSM743105",
  "GSM743106","GSM743107","GSM743108","GSM743112","GSM743117","GSM743120",
  "GSM743122","GSM743125","GSM743126","GSM743128","GSM743130","GSM743131",
  "GSM743135","GSM743138","GSM743139","GSM743141","GSM743146","GSM743158",
  "GSM743159","GSM743161","GSM743164","GSM743169","GSM743170","GSM743171",
  "GSM743173","GSM743175","GSM743178","GSM743179","GSM743182","GSM743183",
  "GSM743184","GSM743185","GSM743186","GSM743189"
)

healthy_samples <- c(
  "GSM743074","GSM743080","GSM743091","GSM743093","GSM743095","GSM743096",
  "GSM743100","GSM743110","GSM743113","GSM743119","GSM743121","GSM743123",
  "GSM743124","GSM743129","GSM743144","GSM743145","GSM743147","GSM743148",
  "GSM743149","GSM743150","GSM743151","GSM743153","GSM743154","GSM743155",
  "GSM743156","GSM743157","GSM743167","GSM743177","GSM743180","GSM743181",
  "GSM743190"
)

trait <- data.frame(
  sampleID = rownames(datExpr),
  group = NA_character_,
  stringsAsFactors = FALSE
)

trait$group[trait$sampleID %in% aml_samples] <- "AML"
trait$group[trait$sampleID %in% healthy_samples] <- "Healthy"

if (any(is.na(trait$group))) {
  stop("Some samples are not assigned to AML/Healthy")
}

trait$AML <- ifelse(trait$group == "AML", 1, 0)
trait$Healthy <- ifelse(trait$group == "Healthy", 1, 0)

rownames(trait) <- trait$sampleID
datTraits <- trait[, c("AML", "Healthy")]

############################################################
## 9) Sample clustering
############################################################

sampleTree <- hclust(dist(datExpr), method = "average")

sample_group_colors <- ifelse(
  trait$group[match(rownames(datExpr), trait$sampleID)] == "AML",
  COL_AML,
  COL_HEALTHY
)

sample_group_colors <- matrix(sample_group_colors, ncol = 1)
colnames(sample_group_colors) <- "Sample group"

save_pdf_base("WGCNA_01_sample_clustering_style.pdf", width = 16, height = 9)

set_base_par(
  mar = c(5, 4.5, 3.5, 2),
  cex = 1.05,
  cex.lab = 1.15,
  cex.axis = 1.05,
  cex.main = 1.40
)

plotDendroAndColors(
  dendro = sampleTree,
  colors = sample_group_colors,
  groupLabels = colnames(sample_group_colors),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Sample clustering"
)

dev.off()

############################################################
## 10) Choose soft-threshold power
############################################################

powers <- c(1:10, seq(12, 24, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]

softPower <- ifelse(
  any(fit >= 0.8),
  sft$fitIndices[which(fit >= 0.8)[1], 1],
  sft$fitIndices[which.max(fit), 1]
)

cat("Selected softPower =", softPower, "\n")

scale_df <- data.frame(
  power = sft$fitIndices[, 1],
  value = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
)

conn_df <- data.frame(
  power = sft$fitIndices[, 1],
  value = sft$fitIndices[, 5]
)

save_pdf_base("WGCNA_02_soft_threshold_2panel_style.pdf", width = 8.4, height = 5.0)

par(
  mfrow = c(1, 2),
  bg = COL_BG,
  mar = c(4.8, 4.8, 2.6, 1.4),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  font = 2,
  font.lab = 2,
  font.axis = 1,
  font.main = 2
)

plot(
  scale_df$power,
  scale_df$value,
  type = "n",
  xlab = "Soft Power Threshold",
  ylab = "Scale-free Topology Model Fit",
  main = "Scale independence",
  cex.main = 1.22,
  cex.lab = 1.10,
  cex.axis = 1.00
)

text(
  scale_df$power,
  scale_df$value,
  labels = scale_df$power,
  cex = 0.98,
  col = "black",
  font = 1
)

abline(v = softPower, lty = 2, lwd = 1.4, col = "black")
abline(h = 0.8, lty = 2, lwd = 1.4, col = "black")

sel_y1 <- scale_df$value[scale_df$power == softPower]

points(
  softPower,
  sel_y1,
  pch = 21,
  bg = "black",
  col = "black",
  cex = 1.6
)

text(
  softPower,
  sel_y1,
  labels = softPower,
  col = "white",
  cex = 0.76,
  font = 1
)

plot(
  conn_df$power,
  conn_df$value,
  type = "n",
  xlab = "Soft Power Threshold",
  ylab = "Mean Connectivity",
  main = "Mean connectivity",
  cex.main = 1.22,
  cex.lab = 1.10,
  cex.axis = 1.00
)

text(
  conn_df$power,
  conn_df$value,
  labels = conn_df$power,
  cex = 0.98,
  col = "black",
  font = 1
)

abline(v = softPower, lty = 2, lwd = 1.4, col = "black")

sel_y2 <- conn_df$value[conn_df$power == softPower]

points(
  softPower,
  sel_y2,
  pch = 21,
  bg = "black",
  col = "black",
  cex = 1.6
)

text(
  softPower,
  sel_y2,
  labels = softPower,
  col = "white",
  cex = 0.76,
  font = 1
)

dev.off()

############################################################
## 11) Network construction & module detection
############################################################

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)

moduleColors <- net$colors
MEs <- orderMEs(net$MEs)

############################################################
## 12) Dendrogram + module colors
############################################################

save_pdf_base("WGCNA_03_gene_dendrogram_modules_style.pdf", width = 14, height = 9.2)

par(
  bg = COL_BG,
  mar = c(6, 6, 5, 2),
  mgp = c(3.0, 1.0, 0),
  font = 2,
  font.lab = 2,
  font.axis = 2,
  font.main = 2,
  cex = 1.20,
  cex.lab = 1.45,
  cex.axis = 1.25,
  cex.main = 1.75
)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors",
  cex.colorLabels = 1.12,
  cex.main = 1.75
)

dev.off()

############################################################
## 13) Module-trait relationship
############################################################

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

row_labels <- colnames(MEs)
col_labels <- c("AML", "Healthy")
row_colors <- sub("^ME", "", row_labels)

mat_show <- moduleTraitCor
rownames(mat_show) <- row_labels
colnames(mat_show) <- col_labels

text_show <- matrix(
  paste0(
    sprintf("%.2f", moduleTraitCor),
    "\n(",
    formatC(moduleTraitP, format = "e", digits = 2),
    ")"
  ),
  nrow = nrow(moduleTraitCor),
  ncol = ncol(moduleTraitCor),
  byrow = FALSE
)

rownames(text_show) <- row_labels
colnames(text_show) <- col_labels

nR <- nrow(mat_show)

save_pdf_base("WGCNA_04_module_trait_heatmap_style.pdf", width = 9.6, height = 9.4)

layout(
  matrix(c(1, 2, 3, 4), nrow = 1),
  widths = c(1.18, 0.16, 2.35, 0.55)
)

# -------------------------------------------------
# panel 1: module names
# -------------------------------------------------
par(
  mar = c(4.2, 0, 4.8, 0),
  bg = COL_BG,
  font = 2,
  font.main = 2,
  font.axis = 2,
  font.lab = 2
)

plot.new()
plot.window(
  xlim = c(0, 1),
  ylim = c(0.5, nR + 0.5),
  xaxs = "i",
  yaxs = "i"
)

for (i in seq_len(nR)) {
  y <- nR - i + 1
  text(
    1.00,
    y,
    labels = row_labels[i],
    adj = 1,
    cex = 1.34,
    col = "black",
    font = 2
  )
}

# -------------------------------------------------
# panel 2: module color bar
# -------------------------------------------------
par(
  mar = c(4.2, 0, 4.8, 0),
  bg = COL_BG,
  font = 2,
  font.main = 2,
  font.axis = 2,
  font.lab = 2
)

plot.new()
plot.window(
  xlim = c(0, 1),
  ylim = c(0.5, nR + 0.5),
  xaxs = "i",
  yaxs = "i"
)

for (i in seq_len(nR)) {
  y <- nR - i + 1
  rect(
    0.10,
    y - 0.5,
    0.90,
    y + 0.5,
    col = row_colors[i],
    border = NA
  )
}

# -------------------------------------------------
# panel 3: heatmap body
# -------------------------------------------------
par(
  mar = c(4.2, 0.02, 4.8, 0.02),
  bg = COL_BG,
  font = 2,
  font.main = 2,
  font.axis = 2,
  font.lab = 2
)

plot.new()

x_centers <- c(1.00, 1.58)
cell_half_width <- 0.29

plot.window(
  xlim = c(0.74, 1.86),
  ylim = c(0.5, nR + 0.5),
  xaxs = "i",
  yaxs = "i"
)

for (i in seq_len(nR)) {
  for (j in seq_len(ncol(mat_show))) {
    
    y <- nR - i + 1
    z <- mat_show[i, j]
    
    idx <- round((z + 1) / 2 * 99) + 1
    idx <- max(1, min(100, idx))
    fill_col <- module_trait_palette[idx]
    
    rect(
      x_centers[j] - cell_half_width,
      y - 0.5,
      x_centers[j] + cell_half_width,
      y + 0.5,
      col = fill_col,
      border = "#D0D0D0",
      lwd = 0.8
    )
    
    text(
      x_centers[j],
      y,
      labels = text_show[i, j],
      cex = 1.25,
      col = "black",
      font = 2
    )
  }
}

axis(1, at = x_centers, labels = FALSE, tick = FALSE)

text(
  x = x_centers,
  y = 0.06,
  labels = col_labels,
  xpd = NA,
  srt = 0,
  font = 2,
  cex = 1.62
)

title(
  main = "Module–trait relationships",
  font.main = 2,
  cex.main = 2.05,
  line = 1
)

# -------------------------------------------------
# panel 4: legend
# -------------------------------------------------
par(
  mar = c(4.2, 0.02, 4.8, 0.40),
  bg = COL_BG,
  font = 2,
  font.main = 2,
  font.axis = 2,
  font.lab = 2
)

plot.new()
plot.window(
  xlim = c(0, 1.2),
  ylim = c(0, 1),
  xaxs = "i",
  yaxs = "i"
)

y0 <- 0.18
y1 <- 0.82
legend_y <- seq(y0, y1, length.out = 200)
legend_vals <- seq(-1, 1, length.out = 199)

for (k in seq_along(legend_vals)) {
  idx <- round((legend_vals[k] + 1) / 2 * 99) + 1
  idx <- max(1, min(100, idx))
  
  rect(
    0.18,
    legend_y[k],
    0.34,
    legend_y[k + 1],
    col = module_trait_palette[idx],
    border = NA
  )
}

rect(0.18, y0, 0.34, y1, border = "black", lwd = 0.9)

tick_vals <- c(-1, -0.5, 0, 0.5, 1)
tick_pos <- (tick_vals + 1) / 2 * (y1 - y0) + y0

segments(0.34, tick_pos, 0.42, tick_pos, col = "black", lwd = 0.8)

text(
  0.48,
  tick_pos,
  labels = c("-1", "-0.5", "0", "0.5", "1"),
  adj = 0,
  cex = 1.00,
  col = "black",
  font = 2
)

dev.off()

############################################################
## 14) Pick the most AML-positively associated module
############################################################

aml_cor <- moduleTraitCor[, "AML"]
aml_p <- moduleTraitP[, "AML"]

candidate_idx <- which(aml_cor > 0 & aml_p < 0.05)

if (length(candidate_idx) == 0) {
  warning("No AML-positive module with p < 0.05. Using the module with the highest AML correlation.")
  candidate_idx <- which(aml_cor > 0)
}

if (length(candidate_idx) == 0) {
  warning("No module has positive AML correlation. Using the module with the largest absolute AML correlation.")
  candidate_idx <- seq_along(aml_cor)
  best_idx <- candidate_idx[which.max(abs(aml_cor[candidate_idx]))]
} else {
  best_idx <- candidate_idx[which.max(aml_cor[candidate_idx])]
}

bestME <- rownames(moduleTraitCor)[best_idx]
bestModule <- sub("^ME", "", bestME)

cat(
  "Best AML-associated module:",
  bestME,
  " r =",
  aml_cor[best_idx],
  " p =",
  aml_p[best_idx],
  "\n"
)

############################################################
## 15) Extract genes in selected module & compute MM / GS
############################################################

moduleGenes <- net$colors == bestModule
moduleGeneNames <- colnames(datExpr)[moduleGenes]

cat("Number of genes in selected module:", length(moduleGeneNames), "\n")

MM <- as.data.frame(cor(datExpr, MEs, use = "p"))
MM_p <- as.data.frame(corPvalueStudent(as.matrix(MM), nrow(datExpr)))

GS <- as.data.frame(cor(datExpr, datTraits$AML, use = "p"))
GS_p <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))

colnames(GS) <- "GS_AML"
colnames(GS_p) <- "pGS_AML"

out <- data.frame(
  gene_symbol = colnames(datExpr),
  moduleColor = net$colors,
  GS,
  GS_p,
  MM,
  MM_p,
  stringsAsFactors = FALSE
)

out_module <- out[out$moduleColor == bestModule, ]

############################################################
## 16) Export module results
############################################################

write.xlsx(
  out_module,
  file = paste0("WGCNA_GSE30029_", bestModule, "_module_genes.xlsx"),
  overwrite = TRUE
)

write.table(
  out_module,
  file = paste0("WGCNA_GSE30029_", bestModule, "_module_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

module_gene_list <- data.frame(
  gene_symbol = out_module$gene_symbol,
  stringsAsFactors = FALSE
)

write.xlsx(
  module_gene_list,
  file = paste0("WGCNA_GSE30029_", bestModule, "_module_gene_list.xlsx"),
  overwrite = TRUE
)

write.table(
  module_gene_list$gene_symbol,
  file = paste0("WGCNA_GSE30029_", bestModule, "_module_gene_list.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
## 17) MM vs GS plot
############################################################

me_col <- match(bestME, colnames(MM))

x_mm <- abs(MM[moduleGenes, me_col])
y_gs <- abs(GS[moduleGenes, 1])

mm_gs_cor <- cor(x_mm, y_gs, use = "p")
mm_gs_p <- corPvalueStudent(mm_gs_cor, nSamples = length(x_mm))

plot_df <- data.frame(
  MM = as.numeric(x_mm),
  GS = as.numeric(y_gs)
)

if (mm_gs_p < 1e-200) {
  main_text <- paste0(
    "MM vs GS in ",
    bestModule,
    " module  cor = ",
    sprintf("%.2f", mm_gs_cor),
    ", P < 1e-200"
  )
} else {
  main_text <- paste0(
    "MM vs GS in ",
    bestModule,
    " module  cor = ",
    sprintf("%.2f", mm_gs_cor),
    ", P = ",
    format(mm_gs_p, scientific = TRUE, digits = 2)
  )
}

p_mm_gs <- ggplot(plot_df, aes(x = MM, y = GS)) +
  geom_point(
    shape = 16,
    size = 1.8,
    color = "black",
    alpha = 0.65
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.25,
    color = COL_AML,
    fill = "#E8B6B2",
    alpha = 0.28
  ) +
  scale_x_continuous(
    limits = c(0, max(plot_df$MM, na.rm = TRUE) * 1.04),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(plot_df$GS, na.rm = TRUE) * 1.04),
    expand = c(0, 0)
  ) +
  labs(
    title = main_text,
    x = paste0("Module Membership in ", bestME),
    y = "Gene Significance for AML"
  ) +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(face = "bold", color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 12)
    ),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, face = "bold"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.7, color = "black")
  )

save_ggplot_all(
  p_mm_gs,
  "WGCNA_05_MM_vs_GS_best_module_style",
  width = 10,
  height = 8,
  dpi = 600
)

############################################################
## 18) Extract MM and GS for selected module
############################################################

ME_target <- MEs[, bestME]

MM_target <- cor(datExpr, ME_target, use = "p")
pMM_target <- corPvalueStudent(MM_target, nrow(datExpr))

GS_target <- cor(datExpr, datTraits$AML, use = "p")
pGS_target <- corPvalueStudent(GS_target, nrow(datExpr))

############################################################
## 19) Build candidate table
############################################################

hubTable <- data.frame(
  gene_symbol = colnames(datExpr),
  moduleColor = net$colors,
  MM = as.numeric(MM_target),
  pMM = as.numeric(pMM_target),
  GS = as.numeric(GS_target),
  pGS = as.numeric(pGS_target),
  stringsAsFactors = FALSE
)

hubTable <- hubTable[hubTable$moduleColor == bestModule, ]
hubTable$absMM <- abs(hubTable$MM)
hubTable$absGS <- abs(hubTable$GS)

############################################################
## 20) hubGenes
############################################################

hubGenes <- hubTable[
  hubTable$absMM >= 0.75 &
    hubTable$absGS >= 0.4 &
    hubTable$pGS < 0.05,
]

hubGenes <- hubGenes[
  order(hubGenes$absMM, hubGenes$absGS, decreasing = TRUE),
]

cat("Number of hub genes identified:", nrow(hubGenes), "\n")
print(head(hubGenes, 10))

write.xlsx(
  hubGenes,
  file = paste0("WGCNA_GSE30029_", bestModule, "_hub_genes_2.xlsx"),
  overwrite = TRUE
)

write.table(
  hubGenes,
  file = paste0("WGCNA_GSE30029_", bestModule, "_hub_genes_2.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

hub_gene_list <- data.frame(
  gene_symbol = hubGenes$gene_symbol,
  stringsAsFactors = FALSE
)

write.xlsx(
  hub_gene_list,
  file = paste0("WGCNA_GSE30029_", bestModule, "_hub_gene_list_2.xlsx"),
  overwrite = TRUE
)

write.table(
  hub_gene_list$gene_symbol,
  file = paste0("WGCNA_GSE30029_", bestModule, "_hub_gene_list_2.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
## 21) core_gene
############################################################

core_gene <- hubTable[
  hubTable$absMM >= 0.80 &
    hubTable$absGS >= 0.40 &
    hubTable$pGS < 0.05,
]

core_gene <- core_gene[
  order(core_gene$absMM, core_gene$absGS, decreasing = TRUE),
]

cat("Number of core genes identified:", nrow(core_gene), "\n")
print(head(core_gene, 10))

write.xlsx(
  core_gene,
  file = paste0("WGCNA_GSE30029_", bestModule, "_core_gene.xlsx"),
  overwrite = TRUE
)

write.table(
  core_gene,
  file = paste0("WGCNA_GSE30029_", bestModule, "_core_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

core_gene_list <- data.frame(
  gene_symbol = core_gene$gene_symbol,
  stringsAsFactors = FALSE
)

write.xlsx(
  core_gene_list,
  file = paste0("WGCNA_GSE30029_", bestModule, "_core_gene_list.xlsx"),
  overwrite = TRUE
)

write.table(
  core_gene_list$gene_symbol,
  file = paste0("WGCNA_GSE30029_", bestModule, "_core_gene_list.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Downstream analysis finished for module:", bestModule, "\n")


# END ORIGINAL SCRIPT P05

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P06: 07_01_SC_annotation_数据清洗.R
# ==============================================================================
rm(list = ls())  # 清空环境变量

# 加载R包----
library(qs)
library(Seurat)
library(ggplot2)
library(data.table)
library(SeuratDisk)
library(decontX)
library(celda)
library(DoubletFinder)
library(harmony)
library(ggsci)
library(clustree)


############################################################-
# 0. 全局参数与路径 ----
############################################################-

AML_TXT_DIR <- "G:\\AML_SC_RNA\\GSE116256\\RAWdata\\新诊AML_未注释"
Normal_TXT_DIR <- "G:\\AML_SC_RNA\\GSE116256\\RAWdata\\健康对照_未注释"
OUT_QC_DIR  <- "G:\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\1.Sc_basic_analysis"

dir.create(OUT_QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ================= 图像辅助 =================
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

.fig_id <- 0L
fig_id <- function() {.fig_id <<- .fig_id + 1L; .fig_id}
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

############################################################-
# 1. 数据读取工具函数 ----
############################################################-
#Seurat读取数据----
##疾病样本单细胞对象创建----
setwd(AML_TXT_DIR)
#存放数据的地址，注意子文件夹的名字必须是barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz不能改动
#Seurat_R
#----treat
#--------P1
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------P2
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------P3
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#----normal
#--------N1
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------N2
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz
#--------N3
#------------barcodes.tsv.gz
#------------features.tsv.gz
#------------matrix.mtx.gz



read_txt_to_seurat <- function(file) {
  message("📥 Reading: ", basename(file))
  dat <- fread(file, data.table = FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  storage.mode(dat) <- "integer"
  
  seu <- CreateSeuratObject(
    counts = dat,
    project = "AML",
    min.cells = 3,
    min.features = 200
  )
  
  sid <- tools::file_path_sans_ext(basename(file))
  seu$type <- sid
  seu$orig.ident <- sid
  return(seu)
}

files <- list.files(AML_TXT_DIR, pattern = "\\.txt$", full.names = TRUE)
seurat_list <- lapply(files, read_txt_to_seurat)
names(seurat_list) <- tools::file_path_sans_ext(basename(files))


# 将所有 Seurat 对象合并成一个 Control 的 Seurat 对象
seurat_aml <- Reduce(function(x, y) merge(x, y), seurat_list)

# 查看合并后的 Seurat 对象
head(seurat_aml@meta.data)
dim(seurat_aml)


#进行数据质控
seurat_aml[["percent.mt"]] <- PercentageFeatureSet(seurat_aml, pattern = "^MT-")#线粒体基因
seurat_aml[["percent.rb"]] <- PercentageFeatureSet(seurat_aml, pattern = "^RP")#核糖体基因
seurat_aml[["percent.hb"]] <- PercentageFeatureSet(seurat_aml, pattern = "^HB[AB]")#血红蛋白基因

randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

p1 <- VlnPlot(
  seurat_aml, 
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "percent.hb"), 
  ncol = 5
) + scale_fill_manual(values = randomColors)

# 修改后的代码，确保提供完整的路径与文件名
pdf(file.path(OUT_QC_DIR, "Fig01_QC_violin_AML.pdf"), width = 20, height = 8)
print(p1)
dev.off()



#设置阈值
min_nCount <- 0
max_nCount <- 8000
min_nFeature <- 0
max_nFeature <- 2250
thr_mt <- 5#过滤高端
# thr_rb <- 50#过滤高端
# thr_hb <- 5#过滤高端

#筛选前
dim_before <- dim(seurat_aml)

#过滤
mask1 <- min_nCount <= seurat_aml$nCount_RNA & seurat_aml$nCount_RNA <= max_nCount
mask2 <- min_nFeature <= seurat_aml$nFeature_RNA &  seurat_aml$nFeature_RNA <= max_nFeature
mask3 <- seurat_aml$percent.mt <= thr_mt
# mask4 <- seurat_pnh$nFeature_RNA <= thr_rb
# mask5 <- seurat_pnh$percent.mt <= thr_hb
seurat_aml <- seurat_aml[, mask1 & mask2 & mask3]

#筛选后
dim_after <- dim(seurat_aml)

# 写入 txt
log_file <- file.path(OUT_QC_DIR, "Fig01_QC_pnh_log.txt")

cat(
  " 时间:\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "阈值:\n",
  "  nCount_RNA: ", min_nCount, " - ", max_nCount, "\n",
  "  nFeature_RNA: ", min_nFeature, " - ", max_nFeature, "\n",
  "  percent.mt <= ", thr_mt, "\n",
  if (exists("thr_rb")) {
    paste0("  percent.rb <= ", thr_rb, "\n")
  } else {
    "  percent.rb: 未筛选\n"
  },
  if (exists("thr_hb")) {
    paste0("  percent.hb <= ", thr_hb, "\n")
  } else {
    "  percent.hb: 未筛选\n"
  },
  "\n",
  "筛选前 dim: ", dim_before[1], " features × ", dim_before[2], " cells\n",
  "筛选后 dim: ", dim_after[1], " features × ", dim_after[2], " cells\n",
  "说明: dim 的第一个数字是基因/特征数，第二个是细胞数\n",
  "----------------------------------------\n",
  file = log_file, append = TRUE
)

qsave(seurat_aml, file.path(OUT_QC_DIR, "1.QC_AML_filtered.qs"))







##正常样本单细胞对象创建----
setwd(Normal_TXT_DIR)

read_txt_to_seurat <- function(file) {
  message("📥 Reading: ", basename(file))
  dat <- fread(file, data.table = FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  storage.mode(dat) <- "integer"
  
  seu <- CreateSeuratObject(
    counts = dat,
    project = "Normal",
    min.cells = 3,
    min.features = 200
  )
  
  sid <- tools::file_path_sans_ext(basename(file))
  seu$type <- sid
  seu$orig.ident <- sid
  return(seu)
}

files <- list.files(Normal_TXT_DIR, pattern = "\\.txt$", full.names = TRUE)
seurat_list <- lapply(files, read_txt_to_seurat)
names(seurat_list) <- tools::file_path_sans_ext(basename(files))


# 将所有 Seurat 对象合并成一个 Control 的 Seurat 对象
seurat_con <- Reduce(function(x, y) merge(x, y), seurat_list)

# 查看合并后的 Seurat 对象
head(seurat_con@meta.data)
dim(seurat_con)

#进行数据质控
seurat_con[["percent.mt"]] <- PercentageFeatureSet(seurat_con, pattern = "^MT-")#线粒体基因
seurat_con[["percent.rb"]] <- PercentageFeatureSet(seurat_con, pattern = "^RP")#核糖体基因
seurat_con[["percent.hb"]] <- PercentageFeatureSet(seurat_con, pattern = "^HB[AB]")#血红蛋白基因

randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(100,randomColor())

p2 <- VlnPlot(
  seurat_con,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "percent.hb"),
  ncol = 5
) + scale_fill_manual(values = randomColors)

pdf(file.path(OUT_QC_DIR, "Fig02_QC_violin_col.pdf"), width = 20, height = 8)
print(p2)
dev.off()


#设置阈值
min_nCount <- 0
max_nCount <- 10000
min_nFeature <- 0
max_nFeature <- 2500
thr_mt <- 5#过滤高端
# thr_rb <- 50#过滤高端
# thr_hb <- 5#过滤高端

#筛选前
dim_before <- dim(seurat_con)

#过滤
mask1 <- min_nCount <= seurat_con$nCount_RNA & seurat_con$nCount_RNA <= max_nCount
mask2 <- min_nFeature <= seurat_con$nFeature_RNA &  seurat_con$nFeature_RNA <= max_nFeature
mask3 <- seurat_con$percent.mt <= thr_mt
# mask4 <- seurat_con$nFeature_RNA <= thr_rb
# mask5 <- seurat_con$percent.mt <= thr_hb
seurat_con <- seurat_con[, mask1 & mask2 & mask3]

#筛选后
dim_after <- dim(seurat_con)

# 写入 txt
log_file <- file.path(OUT_QC_DIR, "Fig02_QC_col_log.txt")

cat(
  " 时间:\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "阈值:\n",
  "  nCount_RNA: ", min_nCount, " - ", max_nCount, "\n",
  "  nFeature_RNA: ", min_nFeature, " - ", max_nFeature, "\n",
  "  percent.mt <= ", thr_mt, "\n",
  if (exists("thr_rb")) {
    paste0("  percent.rb <= ", thr_rb, "\n")
  } else {
    "  percent.rb: 未筛选\n"
  },
  if (exists("thr_hb")) {
    paste0("  percent.hb <= ", thr_hb, "\n")
  } else {
    "  percent.hb: 未筛选\n"
  },
  "\n",
  "筛选前 dim: ", dim_before[1], " features × ", dim_before[2], " cells\n",
  "筛选后 dim: ", dim_after[1], " features × ", dim_after[2], " cells\n",
  "说明: dim 的第一个数字是基因/特征数，第二个是细胞数\n",
  "----------------------------------------\n",
  file = log_file, append = TRUE
)

qsave(seurat_con, file.path(OUT_QC_DIR, "2.QC_COL_filtered.qs"))

##合并疾病和正常的seurat样本----
seurat_obj <- merge(seurat_con, seurat_aml, add.cell.ids = c("Control", "AML"))

# 查看合并后的 Seurat 对象的信息
head(seurat_obj@meta.data)
dim(seurat_obj)


pdf(file.path(OUT_QC_DIR, "Fig03_Con_AML_VlnPlot_QC.pdf"))
VlnPlot(seurat_obj, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.hb"),ncol = 5)+
  scale_fill_manual(values =randomColors)
dev.off()

# V5版本的需要joinLayers
LayerData(seurat_obj, assay = "RNA", layer = "counts")
seurat_obj <- JoinLayers(seurat_obj)
dim(seurat_obj[["RNA"]]$counts )

#保存合并之后的变量
qsave(seurat_obj, file.path(OUT_QC_DIR, "3.Con_AML_Seurat_obj.qs"))

############################################################-
# 4. 标准化 / 归一化 / 高变基因 / 降维 ----
############################################################-

# Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)

plot1 <- VariableFeaturePlot(seurat_obj) + ggtitle(fig_title("Variable features (all)"))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + ggtitle(fig_title("Top variable features (top10)"))
save_pdf(file.path(OUT_QC_DIR, "Fig04_Variable_features.pdf"), plot2, w = 7, h = 6)

# Scale
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# PCA / UMAP
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

p_umap <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("UMAP (post-normalization)"))
save_pdf(file.path(OUT_QC_DIR, "Fig05_UMAP_post_normalization.pdf"), p_umap, w = 7.5, h = 5.5)

p_pca <- DimPlot(seurat_obj, reduction = "pca", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("PCA (post-normalization)"))
save_pdf(file.path(OUT_QC_DIR, "Fig06_PCA_post_normalization.pdf"), p_pca, w = 7.5, h = 5.5)


qsave(seurat_obj, file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))
############################################################-
# 5. 去除环境 RNA 污染 (decontX) ----
############################################################-

# seurat_obj_tmp <- qread(file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))

#2026.1.14质控后先去双再去除RNA污染
seurat_obj_tmp <- qread(file.path(OUT_QC_DIR, "7.Doublet_removed.qs"))

# decontX 接受 counts 矩阵；celda::decontX 返回 SingleCellExperiment
# Seurat v5：必须用 layer = "counts"
counts <- GetAssayData(
  object = seurat_obj_tmp,
  assay  = "RNA",
  layer  = "counts"
)
# 为简化与 Seurat 的衔接，这里直接运行并提取 contamination 向量
message(">>> Running decontX (ambient RNA decontamination) ...")
# decontX 接受 matrix / dgCMatrix
dxc <- celda::decontX(as.matrix(counts))  # 若内存较大，可考虑稀疏矩阵/子集
# 提取污染比例
contam <- dxc$contamination
# 写回 meta.data
seurat_obj_tmp$contamination <- contam

# 去污染前污染分布 (按 UMAP)
p_contam_before <- FeaturePlot(
  seurat_obj_tmp, features = "contamination", raster = FALSE
) + ggtitle(fig_title("Ambient contamination (before filtering)")) +
  theme_bw() + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(0, 1)) # 设置颜色映射范围

# 保存未过滤的污染图
save_pdf(file.path(OUT_QC_DIR, "Fig07_Contamination_before.pdf"), p_contam_before, w = 8, h = 6)

# 根据阈值过滤高污染细胞
# thr_contam <- 0.2 
# 
# sce_filt <- seurat_obj_tmp[, seurat_obj_tmp$contamination < thr_contam]


thr_contam <- quantile(seurat_obj_tmp$contamination, 0.95)

sce_filt <- subset(
  seurat_obj_tmp,
  subset =
    contamination <= thr_contam &
    nCount_RNA >= 500 & nFeature_RNA >= 200 & percent.mt <= 10
)


# 去污染后污染分布 (按 UMAP)
p_contam_after <- FeaturePlot(
  sce_filt, features = "contamination", raster = FALSE
) + ggtitle(fig_title("Ambient contamination (after filtering)")) +
  theme_bw() + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(0, 1)) # 设置颜色映射范围

# 保存过滤后的污染图
save_pdf(file.path(OUT_QC_DIR, "Fig08_Contamination_after.pdf"), p_contam_after, w = 8, h = 6)

dim(sce_filt)
qsave(sce_filt, file.path(OUT_QC_DIR, "5.Decontaminated.qs"))

# ############################################################-
# # 6. 细胞周期分析 ----
# ############################################################-
# 
# seurat_cc <- qread(file.path(OUT_QC_DIR, "5.Decontaminated.qs"))
# 
# # Seurat 自带 cc.genes (若版本不同，请根据 Seurat 文档获取最新基因集)
# s.genes   <- intersect(cc.genes$s.genes, rownames(seurat_cc))
# g2m.genes <- intersect(cc.genes$g2m.genes, rownames(seurat_cc))
# 
# seurat_cc <- CellCycleScoring(seurat_cc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # S/G2M violin
# p_cc_s  <- VlnPlot(seurat_cc, features = "S.Score", pt.size = 0, group.by = 'type') + NoLegend() +
#   ggtitle(fig_title("Cell cycle S.Score (by Type)"))
# p_cc_g2m <- VlnPlot(seurat_cc, features = "G2M.Score", pt.size = 0, group.by = 'type') + NoLegend() +
#   ggtitle(fig_title("Cell cycle G2M.Score (by Type)"))
# save_pdf(file.path(OUT_QC_DIR, "Fig09_CellCycle_S_and_G2M.pdf"), p_cc_s + p_cc_g2m, w = 10, h = 5)
# 
# # 去除细胞周期效应并重降维
# seurat_cc <- ScaleData(seurat_cc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_cc))
# p_umap_before_cc <- DimPlot(seurat_cc, reduction = "umap", group.by = "type") +
#   ggtitle(fig_title("UMAP before re-run (with CC regressed)"))
# seurat_cc_new <- RunPCA(seurat_cc, features = VariableFeatures(seurat_cc))
# seurat_cc_new <- RunUMAP(seurat_cc_new, dims = 1:15, reduction = "pca")
# p_umap_after_cc <- DimPlot(seurat_cc_new, reduction = "umap", group.by = "type") +
#   ggtitle(fig_title("UMAP after regressing cell cycle"))
# save_pdf(file.path(OUT_QC_DIR, "Fig10_CC_effect_before_after.pdf"), p_umap_before_cc + p_umap_after_cc, w = 10, h = 5)
# 
# qsave(seurat_cc_new, file.path(OUT_QC_DIR, "6.CC_regressed.qs"))

############################################################-
# 7. 双细胞分析 (DoubletFinder) ----
############################################################-
# 环境准备
library(scDblFinder)
library(SingleCellExperiment)
# BiocManager::install("scDblFinder")

# 读取对象
# seurat_df <- qread(file.path(OUT_QC_DIR, "6.CC_regressed.qs"))

# 2026.1.14质控后先去双再去除RNA污染
seurat_df <- qread(file.path(OUT_QC_DIR, "4.Normalized_scaled.qs"))


# 聚类用于后续可视化
seurat_df <- FindNeighbors(seurat_df, dims = 1:15)
seurat_df <- FindClusters(seurat_df, resolution = 0.5)

# 分组处理
sce_list <- SplitObject(seurat_df, split.by = "type")

# 收集失败分组
failed_groups <- c()

# 定义分析函数
safe_run_scdblfinder <- function(sce_object, verbose = TRUE) {
  tryCatch({
    log_msg <- function(...) if (verbose) cat(..., "\n")
    
    log_msg("🔄 转换为 SingleCellExperiment")
    sce <- as.SingleCellExperiment(sce_object)
    
    log_msg("🧪 运行 scDblFinder")
    sce <- scDblFinder(sce)
    
    log_msg("✅ 添加回 Seurat 对象中")
    sce_object$scDblFinder.class <- colData(sce)$scDblFinder.class
    sce_object$scDblFinder.score <- colData(sce)$scDblFinder.score
    
    return(sce_object)
  }, error = function(e) {
    message("❗ 捕获错误：", conditionMessage(e))
    return(NULL)
  })
}

# 循环处理每个样本
for (nm in names(sce_list)) {
  cat("\n=== 🔁 Running:", nm, "===\n")
  result <- safe_run_scdblfinder(sce_list[[nm]])
  if (is.null(result)) {
    failed_groups <- c(failed_groups, nm)
  } else {
    sce_list[[nm]] <- result
  }
}

cat("\n✅ 所有分组完成。\n")
if (length(failed_groups) > 0) {
  cat("❌ 以下分组失败或被跳过:\n")
  print(failed_groups)
}

# 合并结果
sce_all <- Reduce(function(x, y) merge(x, y), sce_list)
LayerData(sce_all, assay = "RNA", layer = "counts")
sce_all <- JoinLayers(sce_all)

sce_all <- RunPCA(sce_all, features = VariableFeatures(sce_all))
sce_all <- RunUMAP(sce_all, dims = 1:20)

# 设置颜色：固定 doublet 和 singlet 的配色
df_colors <- c("doublet" = "#F8766D", "singlet" = "#00BFC4")  # 对应 Seurat 默认颜色1和2

# 过滤前 UMAP
p_doublet_before <- DimPlot(sce_all, reduction = "umap", group.by = "scDblFinder.class") +
  scale_color_manual(values = df_colors) +
  ggtitle(fig_title("Doublet classification (before filtering)"))
save_pdf(file.path(OUT_QC_DIR, "Fig11_Doublet_before_filter.pdf"), p_doublet_before, w = 7.5, h = 5.5)

# 过滤后
sce_all_filtered <- subset(sce_all, subset = scDblFinder.class == "singlet")
p_doublet_after <- DimPlot(sce_all_filtered, reduction = "umap", group.by = "scDblFinder.class") +
  scale_color_manual(values = df_colors) +
  ggtitle(fig_title("Doublet classification (after filtering)"))
save_pdf(file.path(OUT_QC_DIR, "Fig12_Doublet_after_filter.pdf"), p_doublet_after, w = 7.5, h = 5.5)

dim(sce_all_filtered)
# 保存结果
qsave(sce_all_filtered, file.path(OUT_QC_DIR, "7.Doublet_removed.qs"))


############################################################-
# 8. 去批次效应 (Harmony) ----
############################################################-

# 并行(按需)
suppressWarnings({
  library(future)
  plan("multisession", workers = 6)  # Windows 推荐 multisession
})
options(future.globals.maxSize = 60 * 1024^3)

seurat_hm <- qread(file.path(OUT_QC_DIR, "5.Decontaminated.qs"))
seurat_hm <- RunHarmony(seurat_hm, group.by.vars = "orig.ident")

# 重新寻找高变基因（去批次后）
seurat_hm <- FindVariableFeatures(seurat_hm, selection.method = "vst", nfeatures = 2000)
top10_hm <- head(VariableFeatures(seurat_hm), 10)
plot_hm1 <- VariableFeaturePlot(seurat_hm) + ggtitle(fig_title("Variable features (after Harmony)"))
plot_hm2 <- LabelPoints(plot = plot_hm1, points = top10_hm, repel = TRUE) + ggtitle(fig_title("Top variable features (Harmony)"))
save_pdf(file.path(OUT_QC_DIR, "Fig13_Variable_features_after_Harmony.pdf"), plot_hm2, w = 7, h = 6)

# Harmony 空间的可视化
p_harmony <- DimPlot(seurat_hm, reduction = "harmony", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("Harmony embedding"))
save_pdf(file.path(OUT_QC_DIR, "Fig14_Harmony_embedding.pdf"), p_harmony, w = 7.5, h = 5.5)

# 在 Harmony 上做 UMAP
seurat_hm <- RunUMAP(seurat_hm, dims = 1:20, reduction = "harmony")
p_umap_hm <- DimPlot(seurat_hm, reduction = "umap", pt.size = 0.5) +
  theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5)) +
  ggtitle(fig_title("UMAP after Harmony"))
save_pdf(file.path(OUT_QC_DIR, "Fig15_UMAP_after_Harmony.pdf"), p_umap_hm, w = 7.5, h = 5.5)

# 按样本分面
save_pdf(file.path(OUT_QC_DIR, "Fig16_UMAP_split_by_orig.ident.pdf"), 
         DimPlot(seurat_hm, reduction = "umap", split.by = "orig.ident", pt.size = 0.5) +
           ggtitle(fig_title("UMAP split by sample (orig.ident)")),
         w = 12, h = 7.5)

qsave(seurat_hm, file.path(OUT_QC_DIR, "8.Harmony_corrected.qs"))

dim(seurat_hm)



############################################################-
# 9. 聚类与分辨率选择 ----
############################################################-
seurat_clu <- seurat_hm
seurat_clu <- qread(file.path(OUT_QC_DIR, "8.Harmony_corrected.qs"))

# 选择PC数：肘图 / 热图
save_pdf(file.path(OUT_QC_DIR, "Fig17_ElbowPlot.pdf"),
         ElbowPlot(seurat_clu, ndims = 50) + ggtitle(fig_title("Elbow plot")), w = 5, h = 4)

save_pdf(file.path(OUT_QC_DIR, "Fig18_DimHeatmap_PC1_20.pdf"),
         DimHeatmap(seurat_clu, dims = 1:20, cells = 1000, balanced = TRUE) +
           ggtitle(fig_title("PC1-20 heatmap")), w = 7.5, h = 11)

# 构图图谱并选择合适分辨率
seuratPC <- 30
seurat_clu <- FindNeighbors(seurat_clu, dims = 1:seuratPC, reduction = "pca")

for (res in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) {
  seurat_clu <- FindClusters(seurat_clu, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
# 查看层级关系
p_tree <- clustree(seurat_clu@meta.data, prefix = "RNA_snn_res.")
save_pdf(file.path(OUT_QC_DIR, "Fig19_Clustree_resolutions.pdf"), p_tree + ggtitle(fig_title("Clustree across resolutions")),
         w = 25, h = 20)

# 选择一个合适的 resolution 进行可视化
Idents(seurat_clu) <- seurat_clu@meta.data[["RNA_snn_res.0.5"]]
p_umap_final <- DimPlot(seurat_clu, reduction = "umap", pt.size = 0.1, label = TRUE) +
  ggtitle(fig_title("Final UMAP (resolution = 0.5)")) +
  theme_classic()
save_pdf(file.path(OUT_QC_DIR, "Fig20_UMAP_final_clusters_res0.5.pdf"), p_umap_final, w = 7.5, h = 5.5)

qsave(seurat_clu, file.path(OUT_QC_DIR, "9.Clustered_final.qs"))

############################################################-
# 10. 小结 ----
############################################################-
# 整个流程：
# 数据导入(多格式) → 质控 → 标准化/归一化/高变/降维 → 去环境RNA污染 → 细胞周期分析与回归 → 双细胞检测与过滤 →
# 去批次效应(Harmony) → 聚类与分辨率选择(Clustree) → 最终 UMAP 聚类图。
#
# 所有图像均按 Fig.xx 编号，图内元素英文；关键状态在 OUT_QC_DIR 下均有 qs 与 pdf/log 落盘。
message(">>> Pipeline finished. Outputs saved to: ", OUT_QC_DIR)


# 2026.1.14  GSE116256 16AML+5NORMAL  整体流程：数据导入(TXT) → 质控 → 标准化/归一化/高变/降维 → 双细胞检测与过滤 → 去环境RNA污染 →
# → 去批次效应(Harmony) → 聚类与分辨率选择(Clustree) → 最终 UMAP 聚类图。




# END ORIGINAL SCRIPT P06

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P07: 07_02_SC_annotation_聚类与注释.R
# ==============================================================================
rm(list = ls())  # 清空环境变量

# 加载R包----
library(qs)
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(openxlsx)
library(patchwork)
library(readr)
library(stringr)

# 
# Sys.which("make")
# install.packages("qs", type = "source")


# 根目录(按需修改)
QC_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\1.Sc_basic_analysis"
ANNO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)

# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}


setwd(ANNO_DIR)
#选择分辨率进行降维----
seurat_obj <- qread(file.path(QC_DIR, "9.Clustered_final.qs"))

head(seurat_obj@meta.data)

px=0.5
seurat_obj <- FindClusters(seurat_obj, resolution = px)

DefaultAssay(seurat_obj)
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# only.pos：只保留上调差异表达的基因
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 新增 'difference' 列，计算 pct.1 - pct.2
seurat_obj.markers$difference <- seurat_obj.markers$pct.1 - seurat_obj.markers$pct.2

# 按照每个簇内的 'difference' 降序排列
seurat_obj.markers_sorted <- seurat_obj.markers %>%
  group_by(cluster) %>%
  arrange(desc(difference)) %>%
  ungroup()


# 创建输出文件名
file_name <- sprintf("Tab01_resolution_%.1f_marker_gene.xlsx", px)
# 保存文件
write.xlsx(seurat_obj.markers_sorted, file = file_name)


# 提取每个细胞簇的前十个标记基因
top_genes_per_cluster <- lapply(unique(seurat_obj.markers_sorted$cluster), function(cluster) {
  top_genes <- seurat_obj.markers_sorted[seurat_obj.markers_sorted$cluster == cluster,]
  top_genes <- head(top_genes$gene, 10)  # 获取前10个基因
  paste(top_genes, collapse = ",")  # 将基因名称拼接成字符串
})
# 为 top_genes_per_cluster 添加簇的名称
names(top_genes_per_cluster) <- unique(seurat_obj.markers_sorted$cluster)

print(head(seurat_obj.markers_sorted, 3))

# 将每个簇的前十个标记基因连接成字符串并保存在txt文件中
txt_content <- sapply(names(top_genes_per_cluster), function(cluster) {
  paste0("cluster_", cluster, ": ", top_genes_per_cluster[[cluster]])
})
# 强制将 txt_content 转换为字符类型
txt_content <- as.character(txt_content)
# 构建自动命名的文件名
txt_filename <- sprintf("Tab02_resolution_%.1f_top10_marker_gene.txt", px)
# 写入文本
writeLines(txt_content, con = txt_filename)


#选择每个聚类前5各基因绘制热图
top5seurat_obj.markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = difference)
# 对于每个群集，挑选出5个差异表达最显著的基因  

col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))

# 创建格式化后的文件名
pdf_filename <- sprintf("Fig01_resolution_%.1f_heatmap_plot.pdf", px)

# 保存热图到 PDF
pdf(file = pdf_filename, width = 22, height = 16)

DoHeatmap(
  seurat_obj,
  features = top5seurat_obj.markers$gene,
  group.colors = col
) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(
    low = '#0099CC', mid = 'white', high = '#CC0033',
    name = 'Z-score'
  )

dev.off()



seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, reduction = "pca")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, reduction = "pca")

# 输出文件名
combined_pdf <- sprintf("Fig02_resolution_%.1f_umap_tsne_plot.pdf", px)
# 保存成一页 PDF
pdf(file = combined_pdf, width = 12, height = 5.5)
# 绘制 UMAP 和 tSNE
p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.1
) +
  ggtitle("UMAP") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right"
  )

p_tsne <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.1
) +
  ggtitle("tSNE") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right"
  )
# 拼图
print(p_umap + p_tsne)
dev.off()
#大数据用UMAP，小数据用TSNE
# UMAP在保留数据的全局结构的同时，更注重保留数据的局部结构，因此相对于TSNE更容易发现群簇的细节信息
# TSNE则在较低维度时能够展示出更好的聚类效果，但随着维度的增加，其聚类效果可能会受到影响


head(seurat_obj@meta.data)
# 构建 resolution 列名
res_col <- if (px == 0.5) {
  "RNA_snn_res.0.5"
} else {
  sprintf("RNA_snn_res.%.1f", px)
}


# 保留指定列
seurat_obj@meta.data <- seurat_obj@meta.data[, c("orig.ident", "type", res_col)]

# 构建保存文件名
rds_file <- sprintf("1.resolution_%.1f_seuratUMAP.qs", px)

# 保存对象
qsave(seurat_obj, rds_file)




#细胞注释----
seurat_obj <- qread("1.resolution_0.5_seuratUMAP.qs")#读取数据

file_path <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\ACT_Annotation results_top 1.txt"

# 正确读取含中文的 GB2312 编码文本
df <- read.delim(
  file(file_path, encoding = "GB18030"),
  stringsAsFactors = FALSE,
  header = TRUE
)
# 提取 cluster 编号（去掉"cluster_"前缀）
df$ClusterNum <- as.numeric(gsub("cluster_", "", df$Cluster))

# 清洗 Cell.Type，保留字母、数字、空格、逗号
df$Cell.Type.Clean <- str_trim(gsub("[^a-zA-Z0-9 ,]", "", df$Cell.Type))

# 分割 whg.overlap 基因为列表
df$GeneList <- strsplit(df$whg.overlap, ",")

df$GeneListTop3 <- lapply(df$GeneList, function(x) {
  if (is.null(x)) return(NULL)
  x <- lapply(x, function(gene) gsub(" ", "", gene))  # 对每个基因删除空格
  x <- unlist(x)  # lapply 会返回 list，要变回字符向量
  head(x, 3)
})

# 排序按 ClusterNum
df <- df[order(df$ClusterNum), ]
# 构造 marker_cluster <- c(...) 文字表示（带注释
marker_lines <- c("marker_cluster <- c(")

# 初始化空向量收集 NA 报告
na_messages <- c()

for (i in seq_len(nrow(df))) {
  cluster_id <- df$ClusterNum[i]
  
  #处理基因列
  genes <- df$GeneListTop3[[i]]
  if (is.null(genes) || all(is.na(genes))) {
    genes <- c("NA", "NA", "NA")
    na_messages <- c(na_messages, sprintf("第 %s 簇的 whg.overlap 是 NA", cluster_id))
  } else {
    genes <- trimws(genes)
    genes[genes == ""] <- "NA"
    genes[is.na(genes)] <- "NA"
  }
  
  # 处理 Cell.Type 
  celltype_clean <- df$Cell.Type.Clean[i]
  if (is.na(celltype_clean) || celltype_clean == "") {
    celltype_clean <- "NA"
    na_messages <- c(na_messages, sprintf("第 %s 簇的 Cell.Type 是 NA", cluster_id))
  }
  
  # 构造输出行
  gene_line <- paste0('  "', paste(genes, collapse = '", "'), '"')
  # 只有不是最后一行，才添加逗号
  if (i < nrow(df)) {
    gene_line <- paste0(gene_line, ",")
  }
  comment <- sprintf(" #%s. %s", cluster_id, celltype_clean)
  marker_lines <- c(marker_lines, paste0(gene_line, comment))
}

marker_lines <- c(marker_lines, ")")

# 写入 marker_cluster 向量格式代码
writeLines(marker_lines, "0.5_marker_cluster_output.R")
cat(paste(marker_lines, collapse = "\n"))

# 打印 NA 提示
if (length(na_messages) > 0) {
  cat("\n# 以下是缺失信息提示：\n")
  for (msg in na_messages) cat(msg, "\n")
}


# FCER1A/CLEC10A/CLEC9A/IL3RA/GZMB/LST1/TYROBP/FCN1
# p6 <- FeaturePlot(seurat_obj, features = c("FCER1A","CLEC10A","CLEC9A","IL3RA","GZMB","LST1","TYROBP","FCN1"))#可以选择不同的基因（为某种细胞类型特征基因），进而初步判断分簇是否与预计相同
# p6
# ggsave("gene_feature.pdf", plot = p6, device = "pdf", width = 15, height = 10)


#文本文件打开生成的0.5_marker_cluster_output.R
#将NA的基因、NA的簇填补完整
#保存

source("1_marker_cluster_output.R")
marker_cluster <- unique(marker_cluster)

# 可视化表达情况
dotplot <- DotPlot(seurat_obj, features = marker_cluster, assay = 'RNA') + coord_flip()
# 构建文件名
dotplot_file <- sprintf("Fig03_resolution_%.1f_dotplot.pdf", px)
# 保存 DotPlot 图
ggsave(dotplot_file, plot = dotplot, width = 15, height = 12)


heatmap <- DoHeatmap(seurat_obj, features = marker_cluster) + NoLegend()
# 构建文件名
heatmap_file <- sprintf("Fig04_resolution_%.1f_marker_heatmap.pdf", px)
# 保存热图
ggsave(heatmap_file, plot = heatmap, width = 15, height = 12)


# 根据 marker 基因手动注释细胞类型
# 提取 cluster 编号
df$ClusterNum <- as.numeric(gsub("cluster_", "", df$Cluster))

# 清洗 Cell.Type（保留字母、数字、空格、逗号）
df$Cell.Type.Clean <- str_trim(gsub("[^a-zA-Z0-9 ,]", "", df$Cell.Type))

# 替换空字符串为 "NA"（真正 NA 会被保留）
df$Cell.Type.Clean[df$Cell.Type.Clean == ""] <- "NA"
df$Cell.Type.Clean[is.na(df$Cell.Type.Clean)] <- "NA"

# 排序
df <- df[order(df$ClusterNum), ]

# 构造 new.cluster.ids 向量文本
cluster_lines <- c("new.cluster.ids <- c(")

# 遍历每一行
for (i in seq_len(nrow(df))) {
  cell_type <- df$Cell.Type.Clean[i]
  cluster_id <- df$ClusterNum[i]
  
  # 构造行文本
  line <- paste0('  "', cell_type, '"')
  
  # 如果不是最后一行，则加上逗号
  if (i < nrow(df)) {
    line <- paste0(line, ",")
  }
  
  # 添加注释
  comment <- sprintf(" #%d", cluster_id)
  cluster_lines <- c(cluster_lines, paste0(line, comment))
}
# 添加结束括号
cluster_lines <- c(cluster_lines, ")")
# 写入文件
writeLines(cluster_lines, "0.5_new_cluster_ids_output.R")
# 打印结果
cat(paste(cluster_lines, collapse = "\n"))


source("1_new_cluster_ids_output.R")
# 更新 Seurat 对象的聚类标识符
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# 拼出列名，如 RNA_snn_res.0.1
cluster_col <- paste0("RNA_snn_res.", px)
# 检查列是否存在
if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("找不到聚类列:", cluster_col))
}
# 将其转为因子后再转为数字，确保整型下标
cluster_ids <- as.integer(as.character(seurat_obj@meta.data[[cluster_col]]))
# 防御性检查：cluster ID 超过定义的细胞类型数量时赋值为 "Unknown"
seurat_obj$celltype <- ifelse(
  !is.na(cluster_ids) & cluster_ids < length(new.cluster.ids),
  new.cluster.ids[cluster_ids + 1],  # R索引从1开始，cluster从0开始
  "Unknown"
)
# 查看结果
head(seurat_obj@meta.data)


#绘制注释后的细胞聚类图----
# 加载所需包
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ggplot2)
library(ggrepel)
library(ggunchull)
library(dplyr)
library(stringr)
library(viridis)
# -- 合并 tSNE 坐标 --
meta <- cbind(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions$tsne@cell.embeddings))
colnames(meta)[(ncol(meta)-1):ncol(meta)] <- c("tSNE_1", "tSNE_2")

# -- 替换 NA 并准备标签 --
meta$cell_type <- as.character(meta$celltype)
meta$cell_type[is.na(meta$celltype)] <- "NA"
meta$main_type_wrapped <- str_wrap(meta$cell_type, width = 20)

# -- 获取标签位置 --
main_type_med <- meta %>%
  group_by(main_type_wrapped) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2),
    .groups = "drop"
  )

# 分配颜色
cell_types <- unique(meta$cell_type)
non_na_types <- setdiff(cell_types, "NA")
cols <- viridis(length(non_na_types), option = "D")
names(cols) <- non_na_types
if ("NA" %in% cell_types) cols <- c(cols, "NA" = "grey80")

x_limits <- c(min(meta$tSNE_1) - 1, max(meta$tSNE_1) + 1)
y_limits <- c(min(meta$tSNE_2) - 1, max(meta$tSNE_2) + 1)

# 绘图
annotated_tSNE <- ggplot(meta, aes(x = tSNE_1, y = tSNE_2)) +
  stat_unchull(aes(fill = cell_type, color = cell_type),
               alpha = 0.15, linewidth = 0.8, linetype = "dashed", delta = 0.3,
               show.legend = FALSE) +
  geom_point(aes(color = cell_type), size = 0.4, show.legend = FALSE) +
  geom_text_repel(data = main_type_med,
                  aes(label = main_type_wrapped),
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.size = 0.3,
                  max.overlaps = 100,
                  color = "black",
                  show.legend = FALSE) +
  # 坐标轴设置（显示标签，不显示刻度）
  scale_x_continuous(name = "tSNE 1", breaks = NULL, limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(name = "tSNE 2", breaks = NULL, limits = y_limits, expand = c(0, 0)) +
  # 主题美化
  theme_minimal(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),  # 放大轴标题
    axis.text = element_blank(),                          # 不显示刻度
    axis.ticks = element_blank(),                         # 不显示刻度线
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 0.6,
                             arrow = arrow(type = "closed", length = unit(0.25, "cm")))  # 添加箭头
  ) +
  ggtitle("")  # 无标题

# 保存图像为 PDF
umap_file <- sprintf("Fig05_resolution_%.1f_celltype_tSNE.pdf", px)
ggsave(umap_file, plot = annotated_tSNE, width = 10, height = 8)




# -- 合并 UMAP 坐标 --
meta <- cbind(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions$umap@cell.embeddings))
colnames(meta)[(ncol(meta)-1):ncol(meta)] <- c("UMAP_1", "UMAP_2")

# -- 替换 NA 并准备标签 --
meta$cell_type <- as.character(meta$celltype)
meta$cell_type[is.na(meta$celltype)] <- "NA"
meta$main_type_wrapped <- str_wrap(meta$cell_type, width = 20)

# -- 获取标签位置 --
main_type_med <- meta %>%
  group_by(main_type_wrapped) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# 分配颜色
cell_types <- unique(meta$cell_type)
non_na_types <- setdiff(cell_types, "NA")
cols <- viridis(length(non_na_types), option = "D")
names(cols) <- non_na_types
if ("NA" %in% cell_types) cols <- c(cols, "NA" = "grey80")

x_limits <- c(min(meta$UMAP_1) - 1, max(meta$UMAP_1) + 1)
y_limits <- c(min(meta$UMAP_2) - 1, max(meta$UMAP_2) + 1)

# 绘图
annotated_umap <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_unchull(aes(fill = cell_type, color = cell_type),
               alpha = 0.15, linewidth = 0.8, linetype = "dashed", delta = 0.3,
               show.legend = FALSE) +
  geom_point(aes(color = cell_type), size = 0.4, show.legend = FALSE) +
  geom_text_repel(data = main_type_med,
                  aes(label = main_type_wrapped),
                  fontface = "bold",
                  size = 5,
                  box.padding = 0.5,
                  segment.size = 0.3,
                  max.overlaps = 100,
                  color = "black",
                  show.legend = FALSE) +
  # 坐标轴设置（显示标签，不显示刻度）
  scale_x_continuous(name = "UMAP 1", breaks = NULL, limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(name = "UMAP 2", breaks = NULL, limits = y_limits, expand = c(0, 0)) +
  # 主题美化
  theme_minimal(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),  # 放大轴标题
    axis.text = element_blank(),                          # 不显示刻度
    axis.ticks = element_blank(),                         # 不显示刻度线
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 0.6,
                             arrow = arrow(type = "closed", length = unit(0.25, "cm")))  # 添加箭头
  ) +
  ggtitle("")  # 无标题

# 保存图像为 PDF
umap_file <- sprintf("Fig05_resolution_%.1f_celltype_umap.pdf", px)
ggsave(umap_file, plot = annotated_umap, width = 10, height = 8)


seurat_obj$type <- ifelse(
  grepl("^AML", seurat_obj$orig.ident),
  "AML",
  "Healthy"
)
table(seurat_obj$type)
table(seurat_obj$celltype)
head(seurat_obj@meta.data[, c("orig.ident", "type", "celltype")])


# 保存seurat变量
qsave(seurat_obj,"2.resolution_0.5_annotation.qs")


#install.packages("wordcloud")

# #细胞聚类圈图----
# library(Seurat)
# library(wordcloud)
# library(circlize)
# library(RColorBrewer)
# source("D:\\plot1cell\\plot1cell\\plot_cell_fraction.R")
# source("D:\\plot1cell\\plot1cell\\plot_circlize.R")
# source("D:\\plot1cell\\plot1cell\\plot_dot.R")
# source("D:\\plot1cell\\plot1cell\\plot_extra.R")
# source("D:\\plot1cell\\plot1cell\\utils.R")
# 
# # 通过网盘分享的文件：plot1cell.zip
# # 链接: https://pan.baidu.com/s/1rHuapWz0CaHaBDwQWASpiw 提取码: 2x86
# # --来自百度网盘超级会员v4的分享
# 
# seurat_object<-qread("2.resolution_0.5_annotation.qs")
# seurat_object$celltype<-seurat_object@active.ident
# 
# colnames(seurat_object@meta.data)
# circ_data <- prepare_circlize_data(seurat_object, scale = 0.8 )
# library(colorspace)
# # 生成随机颜色
# cluster_colors <- rainbow_hcl(length(levels(seurat_object)))
# group_colors <- rainbow_hcl(length(names(table(seurat_object$celltype))))
# rep_colors <- rainbow_hcl(length(names(table(seurat_object$orig.ident))))
# # #绘图
# # png(filename =  '01.circlize_plot.png', width = 8, height = 8,units = 'in', res = 300)
# # plot_circlize(circ_data, do.label = T, pt.size = 0.01,
# #               col.use = cluster_colors, bg.color = 'white',
# #               kde2d.n = 500, repel = TRUE, label.cex = 1,
# #               y_offset = 0,x_offset = -0.03)  # 添加 y 偏移量，正值表示上移
# # # 添加细胞群注释信息
# # add_track(circ_data, group = "celltype", colors = group_colors, track_num = 4)
# # add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num =5)
# # dev.off()
# library(circlize)
# pdf('Fig06_circlize_plot.pdf', width = 8, height = 8)
# plot_circlize(circ_data, do.label = T, pt.size = 0.01,
#               col.use = cluster_colors, bg.color = 'white',
#               kde2d.n = 500, repel = TRUE, label.cex = 1,
#               y_offset = 0,x_offset = 0)  # 添加 y 偏移量，正值表示上移
# add_track(circ_data, group = "celltype", colors = group_colors, track_num = 4)
# add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num =5)
# dev.off()
# 


rm(list = ls())

# =========================================================
# 细胞聚类圈图（最终版）
# 要求：
# 1）圈边只保留2条轨道：cluster -> sample_group
# 2）保留 Cluster 标签
# 3）图注只放一列
# 4）图注上半部分为 celltype（带全称+颜色块）
# 5）图注下半部分为 sample_group
# =========================================================

# ================================
# 0. 加载R包
# ================================
library(qs)
library(Seurat)
library(circlize)
library(colorspace)
library(RColorBrewer)
library(grid)
#install.packages("wordcloud")  # 如果没装
library(wordcloud)

# plot1cell函数
source("D:\\plot1cell\\plot1cell\\plot_cell_fraction.R")
source("D:\\plot1cell\\plot1cell\\plot_circlize.R")
source("D:\\plot1cell\\plot1cell\\plot_dot.R")
source("D:\\plot1cell\\plot1cell\\plot_extra.R")
source("D:\\plot1cell\\plot1cell\\utils.R")

# ================================
# 1. 读取对象
# ================================
seurat_object <- qread("2.resolution_0.5_annotation.qs")

# 使用 active.ident 作为 celltype
seurat_object$celltype <- as.character(seurat_object@active.ident)

# AML / NORMAL 分组
seurat_object$sample_group <- ifelse(
  grepl("^AML", seurat_object$orig.ident, ignore.case = TRUE),
  "AML",
  "Healthy"
)

seurat_object$sample_group <- factor(
  seurat_object$sample_group,
  levels = c("Healthy", "AML")
)

# 保留 celltype 顺序
if (is.factor(seurat_object@active.ident)) {
  seurat_object$celltype <- factor(
    seurat_object$celltype,
    levels = levels(seurat_object@active.ident)
  )
} else {
  seurat_object$celltype <- factor(
    seurat_object$celltype,
    levels = unique(seurat_object$celltype)
  )
}

Idents(seurat_object) <- seurat_object$celltype

# ================================
# 2. 准备circlize数据
# ================================
circ_data <- prepare_circlize_data(seurat_object, scale = 0.8)

# ================================
# 3. 颜色设置
# ================================
# celltype颜色（主图点颜色 + 图注颜色）
n_celltype <- length(levels(seurat_object$celltype))
if (n_celltype <= 12) {
  celltype_colors <- brewer.pal(max(3, n_celltype), "Paired")[seq_len(n_celltype)]
} else {
  celltype_colors <- qualitative_hcl(n_celltype, palette = "Dark 3")
}
names(celltype_colors) <- levels(seurat_object$celltype)

# sample_group颜色
sample_group_colors <- c(
  "Healthy" = "#6A9BCB",
  "AML"    = "#CB5148"
)

# ================================
# 4. celltype全称
# 注意：名称要与你对象中的celltype完全一致
# ================================
celltype_fullname <- c(
  "HSC-like"      = "HSC-like (hematopoietic stem cell-like cell)",
  "LymphoidProg"  = "LymphoidProg (lymphoid progenitor cell)",
  "Progenitor"    = "Progenitor (progenitor cell)",
  "CyclingProg"   = "CyclingProg (cycling progenitor cell)",
  "GMP-like"      = "GMP-like (granulocyte-monocyte progenitor-like cell)",
  "MyeloidProg"   = "MyeloidProg (myeloid progenitor cell)",
  "GranuleProg"   = "GranuleProg (granulocytic progenitor cell)",
  "Promyelo"      = "Promyelocyte",
  "cMono"         = "cMono (classical monocyte)",
  "intMono"       = "intMono (intermediate monocyte)",
  "ncMono"        = "ncMono (non-classical monocyte)",
  "Macrophage"    = "Macrophage",
  "cDC"           = "cDC (conventional dendritic cell)",
  "moDC"          = "moDC (monocyte-derived dendritic cell)",
  "pDC"           = "pDC (plasmacytoid dendritic cell)",
  "NaiveT"        = "NaiveT (naive T cell)",
  "NK_T"          = "NK_T (NK/T cell)",
  "Bcell"         = "Bcell (B cell)",
  "Plasma"        = "Plasma (plasma cell)",
  "Mast"          = "Mast (mast cell)",
  "Erythrocyte"   = "Erythrocyte"
)

# 若有未写入全称的celltype，自动补原名
miss_ct <- setdiff(levels(seurat_object$celltype), names(celltype_fullname))
if (length(miss_ct) > 0) {
  tmp <- setNames(miss_ct, miss_ct)
  celltype_fullname <- c(celltype_fullname, tmp)
}
celltype_fullname <- celltype_fullname[levels(seurat_object$celltype)]

# ================================
# 5. 绘图函数
# ================================
draw_circlize_final <- function() {
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(3.0, 1.6))
  
  # 左侧：主体圈图
  par(mar = c(1, 1, 1, 1))
  circos.clear()
  par(font = 2, cex = 1.1)
  
  plot_circlize(
    circ_data,
    do.label = TRUE,
    pt.size = 0.3,
    col.use = celltype_colors,
    bg.color = "white",
    kde2d.n = 500,
    repel = TRUE,
    label.cex = 1.28,
    y_offset = 0,
    x_offset = 0
  )
  
  # sample_group 更靠外
  add_track(
    circ_data,
    group = "sample_group",
    colors = sample_group_colors,
    track_num = 2
  )
  
  # 右侧：单列图注
  par(mar = c(1, 1, 1, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  
  y0 <- 0.90
  text(
    x = 0.02, y = y0,
    labels = "Cell type",
    adj = c(0, 1),
    cex = 1.15,
    font = 2
  )
  
  yy <- y0 - 0.04
  for (nm in levels(seurat_object$celltype)) {
    rect(
      0.03, yy - 0.012,
      0.07, yy + 0.012,
      col = celltype_colors[nm],
      border = NA
    )
    text(
      0.10, yy,
      labels = celltype_fullname[nm],
      adj = c(0, 0.5),
      cex = 0.74,
      font = 2
    )
    yy <- yy - 0.032
  }
  
  yy <- yy - 0.03
  text(
    x = 0.02, y = yy,
    labels = "Sample group",
    adj = c(0, 1),
    cex = 1.15,
    font = 2
  )
  
  yy <- yy - 0.04
  for (nm in names(sample_group_colors)) {
    rect(
      0.03, yy - 0.012,
      0.07, yy + 0.012,
      col = sample_group_colors[nm],
      border = NA
    )
    text(
      0.10, yy,
      labels = nm,
      adj = c(0, 0.5),
      cex = 0.95,
      font = 2
    )
    yy <- yy - 0.045
  }
}

pdf("Fig06_circlize_plot_final_v3.pdf", width = 14.5, height = 10)
draw_circlize_final()
dev.off()

png(
  filename = "Fig06_circlize_plot_final_v3.png",
  width = 14.5,
  height = 10,
  units = "in",
  res = 300
)
draw_circlize_final()
dev.off()

cat("最终修订版圈图已输出：\n")
cat(" - Fig06_circlize_plot_final_v2.pdf\n")
cat(" - Fig06_circlize_plot_final_v2.png\n")




# seurat_obj$celltype <- gsub(
#   "Naive/Memory T cells",
#   "Naive or Memory T cells",
#   seurat_obj$celltype
# )
# Idents(seurat_obj) <- seurat_obj$celltype
# levels(seurat_obj)

#单独保存各细胞簇----
# 获取meta数据的前几行
head(seurat_obj@meta.data)

# 获取细胞类型的所有levels
cell_types <- levels(seurat_obj)

# 遍历每个细胞类型，进行子集选择和保存
for (cell_type in cell_types) {
  # 创建子集对象
  cell_type_seurat <- subset(seurat_obj, subset = celltype == cell_type)
  
  # 创建保存文件名
  file_name <- paste0("3.", gsub(" ", "_", cell_type), "_seurat.qS")
  
  # 保存为新的 .rds 文件
  qsave(cell_type_seurat, file_name)
}
















# END ORIGINAL SCRIPT P07

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P08: 07_03_SC_annotation_名称修改与marker气泡图绘制.R
# ==============================================================================


rm(list = ls())  # 清空环境变量

# 加载R包----
library(qs)
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(openxlsx)
library(patchwork)
library(readr)
library(stringr)

# 
# Sys.which("make")
# install.packages("qs", type = "source")


# 根目录(按需修改)
#QC_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\1.Sc_basic_analysis"
ANNO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)

# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

setwd(ANNO_DIR)

px=0.5


seurat_obj <- qread("1.resolution_0.5_seuratUMAP.qs")#读取数据

head(seurat_obj@meta.data)


# ================================
# 承接上文：读取对象后继续
# ================================

# 查看聚类列是否存在
colnames(seurat_obj@meta.data)

# 这里你的聚类列名
cluster_col <- "RNA_snn_res.0.5"

if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste0("在 meta.data 中未找到聚类列：", cluster_col))
}

# 转成字符，避免后面映射出问题
seurat_obj[[cluster_col]] <- as.character(seurat_obj[[cluster_col]][, 1])

# ================================
# 1. cluster -> 注释名称映射
# ================================
celltype_map <- c(
  "0"  = "Progenitor",
  "1"  = "cMono",
  "2"  = "Progenitor",
  "3"  = "ncMono",
  "4"  = "GMP-like",
  "5"  = "NaiveT",
  "6"  = "intMono",
  "7"  = "MyeloidProg",
  "8"  = "GranuleProg",
  "9"  = "LymphoidProg",
  "10" = "cDC",
  "11" = "NK_T",
  "12" = "Erythrocyte",
  "13" = "CyclingProg",
  "14" = "Promyelo",
  "15" = "Macrophage",
  "16" = "HSC-like",
  "17" = "moDC",
  "18" = "Mast",
  "19" = "Plasma",
  "20" = "Bcell",
  "21" = "pDC"
)

# 检查对象中的cluster编号
cluster_found <- sort(unique(seurat_obj[[cluster_col]][, 1]))
cat("当前对象中的 cluster 为：\n")
print(cluster_found)

# 检查是否有未定义注释的cluster
cluster_not_mapped <- setdiff(cluster_found, names(celltype_map))
if (length(cluster_not_mapped) > 0) {
  warning("以下 cluster 在 celltype_map 中未定义：",
          paste(cluster_not_mapped, collapse = ", "))
}

# 新建注释列
seurat_obj$celltype_cluster <- unname(celltype_map[seurat_obj[[cluster_col]][, 1]])

# 若有NA，提醒
if (any(is.na(seurat_obj$celltype_cluster))) {
  warning("celltype_cluster 中存在 NA，请检查 cluster 编号与映射关系。")
}

# ================================
# 2. 设置你希望的细胞顺序（按谱系）
# ================================
celltype_order <- c(
  # Stem / Progenitor
  "HSC-like",
  "LymphoidProg",
  "Progenitor",
  "CyclingProg",
  "GMP-like",
  "MyeloidProg",
  "GranuleProg",
  "Promyelo",

  
  # Monocyte / Macrophage
  "cMono",
  "intMono",
  "ncMono",
  "Macrophage",
  
  # DC
  "cDC",
  "moDC",
  "pDC",
  
  # Lymphoid
  "NaiveT",
  "NK_T",
  "Bcell",
  "Plasma",
  
  # Other
  "Mast",
  "Erythrocyte"
)

# 只保留对象里实际存在的类型，避免空level报错
celltype_order_use <- celltype_order[celltype_order %in% unique(seurat_obj$celltype_cluster)]

seurat_obj$celltype_cluster <- factor(
  seurat_obj$celltype_cluster,
  levels = celltype_order_use
)

# 设置身份类
Idents(seurat_obj) <- seurat_obj$celltype_cluster

cat("最终细胞注释顺序：\n")
print(levels(seurat_obj$celltype_cluster))


# ================================
# 3. 定义 marker（按模块排序）
# ================================
marker_genes <- c(
  # ===== Stem / Progenitor =====
  "FAM30A", "MECOM", "CD34",
  "CD3D",
  "NRIP1", "MSI2",
  "PCLAF", "TOP2A", "PCNA",
  "POU4F1", "MPO", "RUNX1T1",
  "CD36", "HLA-DRB5",
  "CPA3", 
  "CXCL8", "ITM2A",
  
  # ===== Monocyte / Macrophage =====
  "MARCKS", "S100A8", "VCAN",
  "MNDA", "IFI30", "MS4A6A",
  "S100A9",
  "CD163", "TNFAIP2",
  
  # ===== DC =====
  "NAPSB", "JAML",
  "HLA-DQB1", "HLA-DQA1",
  "UGCG", "IRF4", "CCDC50",
  
  # ===== Lymphoid =====
  "IL7R", "ETS1", "LTB",
  "CCL5", "GZMA", "GNLY",
  "MS4A1", "BANK1", "BIRC3",
  "MZB1", "SEC11C", "SLAMF7",
  
  # ===== Other =====
  "TPSB2", "TPSD1", "CD69",
  "HBB", "HBD", "CA1"
)

marker_genes <- unique(marker_genes)

# 只保留对象中实际存在的基因
gene_use <- marker_genes[marker_genes %in% rownames(seurat_obj)]
gene_missing <- setdiff(marker_genes, rownames(seurat_obj))

cat("DotPlot可用marker数：", length(gene_use), "\n")
if (length(gene_missing) > 0) {
  cat("以下marker在对象中未找到，将自动跳过：\n")
  print(gene_missing)
}

if (length(gene_use) == 0) {
  stop("一个可用marker都没有，请检查基因名格式（大小写/物种/是否为symbol）。")
}



# ================================
# 6. 画注释后的 DotPlot
# ================================
p_dot <- DotPlot(
  object = seurat_obj,
  features = gene_use,
  group.by = "celltype_cluster",
  assay = DefaultAssay(seurat_obj)
) +
  coord_flip() +
  scale_color_gradientn(
    colours = c("#E0E0E0", "#6BAED6", "#08519C")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9, face = "italic"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  ) +
  labs(
    title = "Marker gene expression across annotated cell populations",
    x = "Marker genes",        
    y = "Cell populations"     
  )


save_pdf(
  file = paste0("Fig", px, "_marker_DotPlot_annotated.pdf"),
  expr = p_dot,
  w = 10,
  h = 16
)

ggsave(
  filename = paste0("Fig", px, "_marker_DotPlot_annotated.png"),
  plot = p_dot,
  width = 10,
  height = 16,
  dpi = 300
)
print(p_dot)







# END ORIGINAL SCRIPT P08

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P09: 07_SC_annotation/1_marker_cluster_output.R
# ==============================================================================
marker_cluster <- c(
  "NRIP1", "AREG","MSI2", #0. Progenitor-like cell (Progenitor)
  "MARCKS", "S100A8", "VCAN", #1. Classical Monocyte (cMono)
  "AREG", "SPINK2", "MSI2", #2. Progenitor-like cell  (Progenitor)
  "MNDA", "IFI30", "MS4A6A", #3. Non-classical Monocyte (ncMono)
  "POU4F1", "MPO", "RUNX1T1","CD34",#4. Granulocyte monocyte progenitor cell(GMP-like)
  "IL7R", "ETS1", "LTB", #5. Naive T cell  (NaiveT)
  "MNDA", "S100A8", "S100A9", #6.  Intermediate Monocyte (intMono)
  "CD36", "SPINK2", "HLA-DRB5", #7.  Myeloid Progenitor  (MyeloidProg)
  "SPINK2", "CPA3","PRSS57", #8. Granule Progenitor （GranuleProg）  
  "CD3D", "IFITM2", "PAN3", #9. Lymphoid Progenitor  （LymphoidProg） 
  "NAPSB", "HLA-DRB5", "JAML", #10. Conventional Dendritic cell(cDC)
  "CCL5", "GZMA", "GNLY", #11. Nature Killer or T cell(NK_T)
  "HBB", "CA1", "HBD", #12. Erythrocyte
  "PCLAF", "TOP2A", "PCNA", #13. Cycling Progenitor  （CyclingProg）
  "CXCL8", "ITM2A", #14. Promyelocyte-like cell （Promyelo）
  "CD163", "VCAN", "TNFAIP2", #15. Macrophage
  "FAM30A", "CD34", "MECOM", #16. Hematopoietic stem cell(HSC-like)
  "HLA-DQB1", "ATP2B1", "HLA-DQA1", #17. Monocyte Derived Dendritic Cell(moDC)
  "TPSB2", "CD69", "TPSD1", #18. Mast cell
  "MZB1", "SEC11C", "SLAMF7", #19. Plasma cell
  "MS4A1", "BANK1", "BIRC3",#20. B cell
  "UGCG", "IRF4", "CCDC50" #21. Plasmacytoid dendritic cell(pDC)
)

# END ORIGINAL SCRIPT P09

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P10: 07_SC_annotation/1_new_cluster_ids_output.R
# ==============================================================================
new.cluster.ids <- c(
  "Progenitor",
  "cMono",
  "Progenitor",
  "ncMono",
  "GMP-like",
  "NaiveT",
  "intMono",
  "MyeloidProg",
  "GranuleProg",
  "LymphoidProg",
  "cDC",
  "NK_T",
  "Erythrocyte",
  "CyclingProg",
  "Promyelo",
  "Macrophage",
  "HSC-like",
  "moDC",
  "Mast",
  "Plasma",
  "Bcell",
  "pDC"
)

# END ORIGINAL SCRIPT P10

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P11: 07_SC_annotation/plot_annotation_pub.R
# ==============================================================================
############################################################
## 01_plot_annotation_pub.R
## Publication-style annotation UMAP/tSNE plots
## Revised version: visible clusters + short axes
##
## Main style:
##   1) 21 short cell-type names
##   2) lineage-aware unified colors from 00_pub_theme_colors.R
##   3) visible but still small cell points
##   4) bold cell labels, moderately enlarged
##   5) no rectangular border / no hull boundary
##   6) short UMAP1/UMAP2 or TSNE1/TSNE2 axes at lower-left corner
############################################################

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
})

source("00_pub_theme_colors.R")

############################################################
## 1. Parameters
############################################################

OUT_DIR <- file.path(ANNO_DIR, "pub_figures_annotation")
safe_mkdir(OUT_DIR)

CELLTYPE_COL <- "celltype"
CELLTYPE_SHORT_COL <- "celltype_short"
GROUP_COL_CANDIDATES <- c("type", "group", "orig.ident", "sample", "sample_id")

## Key revised parameters:
## Previous points were too small/faint. These values keep points fine but visible.
POINT_SIZE_CELLTYPE <- 0.45
POINT_ALPHA_CELLTYPE <- 0.95
LABEL_SIZE_CELLTYPE <- 4.1
LABEL_FONTFACE <- "bold"

POINT_SIZE_GROUP <- 0.28
POINT_ALPHA_GROUP <- 0.80

FIG_W_SINGLE <- 6.4
FIG_H_SINGLE <- 5.8
FIG_W_COMBO <- 11.8
FIG_H_COMBO <- 5.8

############################################################
## 2. Load annotated Seurat object
############################################################

if (!file.exists(SEURAT_ANNO_FILE)) {
  stop("找不到注释后 Seurat 对象：", SEURAT_ANNO_FILE,
       "\n请检查 00_pub_theme_colors.R 中 PROJECT_DIR / ANNO_DIR / SEURAT_ANNO_FILE 设置。")
}

seurat_obj <- qread(SEURAT_ANNO_FILE)

if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  stop("meta.data 中找不到列：", CELLTYPE_COL,
       "\n当前 meta.data 列名包括：", paste(colnames(seurat_obj@meta.data), collapse = ", "))
}

seurat_obj <- add_celltype_short(
  seurat_obj,
  celltype_col = CELLTYPE_COL,
  short_col = CELLTYPE_SHORT_COL
)

seurat_obj@meta.data[[CELLTYPE_SHORT_COL]] <- factor(
  as.character(seurat_obj@meta.data[[CELLTYPE_SHORT_COL]]),
  levels = celltype_levels_short
)

celltypes_use <- levels(droplevels(seurat_obj@meta.data[[CELLTYPE_SHORT_COL]]))
cell_cols_use <- get_celltype_colors(celltypes_use, use_short = TRUE)

write.csv(
  data.frame(celltype = names(cell_cols_use), color = unname(cell_cols_use)),
  file.path(OUT_DIR, "Tab_celltype_colors_used.csv"),
  row.names = FALSE
)

message("Cell types used: ", paste(celltypes_use, collapse = ", "))
message("Output directory: ", OUT_DIR)

############################################################
## 3. Helper functions
############################################################

reduction_exists <- function(obj, reduction) {
  reduction %in% names(obj@reductions)
}

get_reduction_df <- function(obj, reduction, group_col = NULL) {
  if (!reduction_exists(obj, reduction)) {
    stop("Seurat object does not contain reduction: ", reduction)
  }
  emb <- as.data.frame(Embeddings(obj, reduction = reduction)[, 1:2])
  colnames(emb) <- c("dim1", "dim2")
  emb$cell <- rownames(emb)
  emb$celltype_short <- obj@meta.data[emb$cell, CELLTYPE_SHORT_COL, drop = TRUE]
  emb$celltype_short <- factor(as.character(emb$celltype_short), levels = celltype_levels_short)
  if (!is.null(group_col) && group_col %in% colnames(obj@meta.data)) {
    emb$group <- obj@meta.data[emb$cell, group_col, drop = TRUE]
  } else {
    emb$group <- NA_character_
  }
  emb
}

get_label_df <- function(plot_df) {
  plot_df %>%
    dplyr::filter(!is.na(celltype_short)) %>%
    dplyr::group_by(celltype_short) %>%
    dplyr::summarise(
      dim1 = median(dim1, na.rm = TRUE),
      dim2 = median(dim2, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n > 0)
}

add_short_axes <- function(p,
                           plot_df,
                           x_lab = "UMAP1",
                           y_lab = "UMAP2",
                           axis_frac = 0.125,
                           offset_frac = -0.025,
                           text_size = 3.7,
                           line_width = 0.55) {
  xr <- range(plot_df$dim1, na.rm = TRUE)
  yr <- range(plot_df$dim2, na.rm = TRUE)
  x_span <- diff(xr)
  y_span <- diff(yr)

  x0 <- xr[1] + offset_frac * x_span
  y0 <- yr[1] + offset_frac * y_span
  x1 <- x0 + axis_frac * x_span
  y1 <- y0 + axis_frac * y_span

  p +
    annotate("segment", x = x0, xend = x1, y = y0, yend = y0,
             linewidth = line_width, color = "black",
             arrow = arrow(length = unit(0.085, "inches"), type = "closed")) +
    annotate("segment", x = x0, xend = x0, y = y0, yend = y1,
             linewidth = line_width, color = "black",
             arrow = arrow(length = unit(0.085, "inches"), type = "closed")) +
    annotate("text", x = (x0 + x1) / 2, y = y0 - 0.036 * y_span,
             label = x_lab, size = text_size, fontface = "bold", color = "black") +
    annotate("text", x = x0 - 0.036 * x_span, y = (y0 + y1) / 2,
             label = y_lab, size = text_size, fontface = "bold", color = "black", angle = 90)
}

plot_celltype_reduction_shortaxis <- function(obj,
                                               reduction = "umap",
                                               title_text = NULL,
                                               point_size = POINT_SIZE_CELLTYPE,
                                               point_alpha = POINT_ALPHA_CELLTYPE,
                                               label_size = LABEL_SIZE_CELLTYPE,
                                               show_legend = FALSE) {
  df <- get_reduction_df(obj, reduction = reduction, group_col = NULL)
  df <- df %>% dplyr::filter(!is.na(celltype_short))
  lab_df <- get_label_df(df)
  cols <- get_celltype_colors(levels(droplevels(df$celltype_short)), use_short = TRUE)

  red_upper <- toupper(reduction)
  if (is.null(title_text)) title_text <- paste0("Cell annotation")

  p <- ggplot(df, aes(x = dim1, y = dim2, color = celltype_short)) +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0, shape = 16) +
    geom_text_repel(
      data = lab_df,
      aes(x = dim1, y = dim2, label = celltype_short),
      inherit.aes = FALSE,
      size = label_size,
      fontface = LABEL_FONTFACE,
      color = "black",
      box.padding = 0.25,
      point.padding = 0.18,
      min.segment.length = Inf,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_color_manual(
      values = cols,
      drop = FALSE,
      name = "Cell type"
    ) +
    guides(
      color = guide_legend(
        title = "Cell type",
        ncol = 1,
        override.aes = list(
          size = 3.0,
          alpha = 1
        )
      )
    ) +
    coord_equal(clip = "off") +
    labs(title = NULL, color = "Cell type") +
    theme_dim_shortaxis(base_size = 12) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 11, face = "bold", color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      legend.key.size = grid::unit(0.32, "cm"),
      legend.spacing.y = grid::unit(0.00, "cm"),
      legend.box.spacing = grid::unit(0.12, "cm"),
      plot.title = element_blank(),
      plot.margin = margin(8, 8, 12, 14)
    )

  add_short_axes(
    p,
    plot_df = df,
    x_lab = paste0(red_upper, "1"),
    y_lab = paste0(red_upper, "2")
  )
}

plot_group_reduction_shortaxis <- function(obj,
                                            reduction = "umap",
                                            group_col,
                                            title_text = NULL,
                                            point_size = POINT_SIZE_GROUP,
                                            point_alpha = POINT_ALPHA_GROUP) {
  df <- get_reduction_df(obj, reduction = reduction, group_col = group_col)
  df$group <- factor(as.character(df$group))
  groups_use <- levels(droplevels(df$group))
  cols <- get_named_colors(groups_use, palette_named = group_colors)

  red_upper <- toupper(reduction)
  if (is.null(title_text)) title_text <- paste0(group_col)

  p <- ggplot(df, aes(x = dim1, y = dim2, color = group)) +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0, shape = 16) +
    scale_color_manual(values = cols, drop = FALSE, name = group_col) +
    guides(
      color = guide_legend(
        title = group_col,
        override.aes = list(size = 3.0, alpha = 1)
      )
    ) +
    coord_equal(clip = "off") +
    labs(title = NULL, color = group_col) +
    theme_dim_shortaxis(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      legend.key.size = grid::unit(0.32, "cm"),
      legend.spacing.y = grid::unit(0.00, "cm"),
      legend.box.spacing = grid::unit(0.12, "cm"),
      plot.title = element_blank(),
      plot.margin = margin(8, 8, 12, 14)
    )

  add_short_axes(
    p,
    plot_df = df,
    x_lab = paste0(red_upper, "1"),
    y_lab = paste0(red_upper, "2")
  )
}

############################################################
## 4. Export UMAP/tSNE annotation plots
############################################################

if (reduction_exists(seurat_obj, "umap")) {
  p_umap_noleg <- plot_celltype_reduction_shortaxis(
    seurat_obj,
    reduction = "umap",
    title_text = "Cell annotation",
    show_legend = FALSE
  )
  save_plot_pub(p_umap_noleg, file.path(OUT_DIR, "Fig01_celltype_UMAP_shortaxis_visible_noLegend"),
                width = FIG_W_SINGLE, height = FIG_H_SINGLE)

  p_umap <- plot_celltype_reduction_shortaxis(
    seurat_obj,
    reduction = "umap",
    title_text = "Cell annotation",
    show_legend = TRUE
  )
  save_plot_pub(p_umap, file.path(OUT_DIR, "Fig01_celltype_UMAP_shortaxis_visible"),
                width = FIG_W_SINGLE + 0.7, height = FIG_H_SINGLE)
}

if (reduction_exists(seurat_obj, "tsne")) {
  p_tsne_noleg <- plot_celltype_reduction_shortaxis(
    seurat_obj,
    reduction = "tsne",
    title_text = "Cell annotation",
    show_legend = FALSE
  )
  save_plot_pub(p_tsne_noleg, file.path(OUT_DIR, "Fig02_celltype_tSNE_shortaxis_visible_noLegend"),
                width = FIG_W_SINGLE, height = FIG_H_SINGLE)

  p_tsne <- plot_celltype_reduction_shortaxis(
    seurat_obj,
    reduction = "tsne",
    title_text = "Cell annotation",
    show_legend = TRUE
  )
  save_plot_pub(p_tsne, file.path(OUT_DIR, "Fig02_celltype_tSNE_shortaxis_visible"),
                width = FIG_W_SINGLE + 0.7, height = FIG_H_SINGLE)
}

if (exists("p_umap_noleg") && exists("p_tsne_noleg")) {
  p_combo <- p_umap_noleg + p_tsne_noleg + plot_layout(ncol = 2)
  save_plot_pub(p_combo, file.path(OUT_DIR, "Fig03_celltype_UMAP_tSNE_shortaxis_visible"),
                width = FIG_W_COMBO, height = FIG_H_COMBO)
}

############################################################
## 5. Export group UMAP if group/type column exists
############################################################

group_col <- GROUP_COL_CANDIDATES[GROUP_COL_CANDIDATES %in% colnames(seurat_obj@meta.data)][1]

if (!is.na(group_col) && reduction_exists(seurat_obj, "umap")) {
  p_group_umap <- plot_group_reduction_shortaxis(
    seurat_obj,
    reduction = "umap",
    group_col = group_col,
    title_text = paste0("Sample group: ", group_col)
  )
  save_plot_pub(p_group_umap, file.path(OUT_DIR, paste0("Fig04_group_UMAP_", group_col, "_shortaxis_visible")),
                width = FIG_W_SINGLE + 0.7, height = FIG_H_SINGLE)
}

############################################################
## 6. Export metadata summary for checking
############################################################

meta_summary <- seurat_obj@meta.data %>%
  dplyr::mutate(celltype_short = as.character(.data[[CELLTYPE_SHORT_COL]])) %>%
  dplyr::count(celltype_short, name = "n_cells") %>%
  dplyr::arrange(factor(celltype_short, levels = celltype_levels_short))

write.csv(meta_summary, file.path(OUT_DIR, "Tab_celltype_cell_counts.csv"), row.names = FALSE)

message("Annotation plots exported to: ", OUT_DIR)


# END ORIGINAL SCRIPT P11

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P12: 08_SC_cell_ratio.R
# ==============================================================================
rm(list = ls())

# =========================================================
# 单细胞细胞比例分析（仅按 AML / NORMAL 分组）
# 功能：
# 1）读取注释后的 Seurat 对象
# 2）使用 celltype 作为细胞类型分组
# 3）新增 sample_group：AML / NORMAL
# 4）仅基于 AML / NORMAL 两组做细胞比例分析
# 5）输出统计表和图形
# =========================================================

# ================================
# 0. 加载R包
# ================================
library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(scales)
library(RColorBrewer)
library(plotrix)

# ================================
# 1. 路径设置
# ================================
ANNO_DIR  <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
RATIO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\3.Sc_cell_ratio\\new"

if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(RATIO_DIR)) dir.create(RATIO_DIR, recursive = TRUE)

setwd(ANNO_DIR)

# ================================
# 2. 读取Seurat对象
# ================================
seurat_obj <- qread("2.resolution_0.5_annotation.qs")

# 如果 meta.data 中没有 celltype，则使用 active.ident
if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$celltype <- as.character(Idents(seurat_obj))
}

seurat_obj$celltype <- as.character(seurat_obj$celltype)

# 保留 celltype 顺序
if (is.factor(seurat_obj@meta.data$celltype)) {
  celltype_levels <- levels(seurat_obj@meta.data$celltype)
} else if (is.factor(Idents(seurat_obj))) {
  celltype_levels <- levels(Idents(seurat_obj))
} else {
  celltype_levels <- unique(seurat_obj$celltype)
}
seurat_obj$celltype <- factor(seurat_obj$celltype, levels = celltype_levels)

# ================================
# 3. 新增 AML / NORMAL 分组
# ================================
# 规则：
# 样本名以 AML 开头 -> AML
# 其他样本（BM、Control等） -> NORMAL
seurat_obj$sample_group <- ifelse(
  grepl("^AML", seurat_obj$orig.ident, ignore.case = TRUE),
  "AML",
  "Healthy"
)

seurat_obj$sample_group <- factor(
  seurat_obj$sample_group,
  levels = c("Healthy", "AML")
)

cat("细胞类型统计：\n")
print(table(seurat_obj$celltype))

cat("\nAML / NORMAL 分组统计：\n")
print(table(seurat_obj$sample_group))

cat("\norig.ident 与 sample_group 对应关系：\n")
print(table(seurat_obj$orig.ident, seurat_obj$sample_group))

# ================================
# 4. 提取元数据
# ================================
meta <- seurat_obj@meta.data %>%
  dplyr::select(celltype, sample_group) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_levels),
    sample_group = factor(sample_group, levels = c("Healthy", "AML"))
  )

setwd(RATIO_DIR)

write.csv(meta, "0.cell_metadata_AML_Normal.csv", row.names = FALSE)

# ================================
# 5. 统计 AML / NORMAL 中各细胞类型数量
# ================================
group_celltype_count <- meta %>%
  count(sample_group, celltype, name = "CellNumber")

# 计算比例
group_celltype_ratio <- group_celltype_count %>%
  group_by(sample_group) %>%
  mutate(
    TotalCells = sum(CellNumber),
    Proportion = CellNumber / TotalCells
  ) %>%
  ungroup()

write.csv(group_celltype_count, "1.group_celltype_count.csv", row.names = FALSE)
write.csv(group_celltype_ratio, "2.group_celltype_ratio.csv", row.names = FALSE)

cat("\nAML / NORMAL 各细胞类型数量：\n")
print(group_celltype_count)

cat("\nAML / NORMAL 各细胞类型比例：\n")
print(group_celltype_ratio)

# ================================
# 6. 生成宽表（便于看结果）
# ================================
count_wide <- group_celltype_count %>%
  tidyr::pivot_wider(
    names_from = sample_group,
    values_from = CellNumber,
    values_fill = 0
  )

ratio_wide <- group_celltype_ratio %>%
  dplyr::select(sample_group, celltype, Proportion) %>%
  tidyr::pivot_wider(
    names_from = sample_group,
    values_from = Proportion,
    values_fill = 0
  )

write.csv(count_wide, "3.group_celltype_count_wide.csv", row.names = FALSE)
write.csv(ratio_wide, "4.group_celltype_ratio_wide.csv", row.names = FALSE)

# ================================
# 7. 配色
# ================================
n_celltype <- length(celltype_levels)

if (n_celltype <= 12) {
  celltype_colors <- brewer.pal(max(3, n_celltype), "Paired")[seq_len(n_celltype)]
} else {
  celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltype)
}
names(celltype_colors) <- celltype_levels

group_colors <- c(
  "Healthy" = "#6A9BCB",
  "AML"    = "#CB5148"
)

# ================================
# 8. 图1：AML / NORMAL 细胞组成堆叠比例图
# ================================
p1 <- ggplot(group_celltype_ratio, aes(x = sample_group, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Cell-type composition in AML and Healthy"
  )

ggsave("5.AML_NORMAL_stacked_ratio.pdf", p1, width = 8, height = 6)
ggsave("5.AML_NORMAL_stacked_ratio.png", p1, width = 8, height = 6, dpi = 300)

# ================================
# 9. 图2：AML / NORMAL 细胞绝对数量堆叠图
# ================================
p2 <- ggplot(group_celltype_count, aes(x = sample_group, y = CellNumber, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Group",
    y = "Cell number",
    fill = "Cell type",
    title = "Absolute cell counts in AML and Healthy"
  )

ggsave("6.AML_NORMAL_stacked_count.pdf", p2, width = 8, height = 6)
ggsave("6.AML_NORMAL_stacked_count.png", p2, width = 8, height = 6, dpi = 300)

# ================================
# 10. 图3：各 celltype 在 AML / NORMAL 中的比例对比柱状图
# ================================
p3 <- ggplot(group_celltype_ratio, aes(x = celltype, y = Proportion, fill = sample_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Cell type",
    y = "Cell proportion",
    fill = "Group",
    title = "Cell-type proportions in AML vs Healthy"
  )

ggsave("7.AML_NORMAL_celltype_ratio_barplot.pdf", p3, width = 12, height = 7)
ggsave("7.AML_NORMAL_celltype_ratio_barplot.png", p3, width = 12, height = 7, dpi = 300)

# ================================
# 11. 图4：各 celltype 在 AML / NORMAL 中的绝对数量对比柱状图
# ================================
p4 <- ggplot(group_celltype_count, aes(x = celltype, y = CellNumber, fill = sample_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Cell type",
    y = "Cell number",
    fill = "Group",
    title = "Cell-type counts in AML vs Healthy"
  )

ggsave("8.AML_NORMAL_celltype_count_barplot.pdf", p4, width = 12, height = 7)
ggsave("8.AML_NORMAL_celltype_count_barplot.png", p4, width = 12, height = 7, dpi = 300)

# ================================
# 12. 甜甜圈图：AML 和 NORMAL 分开画
# ================================
doughnut <- function(x, labels = names(x), edges = 200, outer.radius = 0.8,
                     inner.radius = 0.6, clockwise = FALSE,
                     init.angle = if (clockwise) 90 else 0, density = NULL,
                     angle = 45, col = NULL, border = FALSE, lty = NULL,
                     main = NULL, ...) {
  if (!is.numeric(x) || any(is.na(x) | x < 0)) stop("'x' values must be positive.")
  if (is.null(labels)) labels <- as.character(seq_along(x)) else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x) / sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) xlim <- (pin[1L] / pin[2L]) * xlim else ylim <- (pin[2L] / pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) col <- if (is.null(density)) palette() else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise) -2 * pi else 2 * pi
  t2xy <- function(t, radius) {
    t2p <- twopi * t + init.angle * pi / 180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n), outer.radius)
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i],
            border = border[i], col = col[i], lty = lty[i])
    Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
      text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
           xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0), ...)
    }
    Pin <- t2xy(seq.int(0, 1, length.out = n * nx), inner.radius)
    polygon(Pin$x, Pin$y, density = density[i], angle = angle[i],
            border = border[i], col = "white", lty = lty[i])
  }
  title(main = main, ...)
  invisible(NULL)
}

for (grp in levels(meta$sample_group)) {
  df_grp <- group_celltype_count %>%
    filter(sample_group == grp) %>%
    mutate(celltype = factor(celltype, levels = celltype_levels)) %>%
    arrange(celltype)
  
  labs_grp <- paste0(
    df_grp$celltype, " (",
    round(df_grp$CellNumber / sum(df_grp$CellNumber) * 100, 2), "%)"
  )
  
  pdf(paste0("9.doughnut_", grp, ".pdf"), width = 12, height = 9)
  doughnut(
    df_grp$CellNumber,
    labels = labs_grp,
    init.angle = 0,
    col = celltype_colors[as.character(df_grp$celltype)],
    border = "white",
    inner.radius = 0.45,
    cex = 0.85,
    main = paste0(grp, " cell-type composition")
  )
  dev.off()
}

# ================================
# 13. 饼图：AML 和 NORMAL 分开画
# ================================
for (grp in levels(meta$sample_group)) {
  df_grp <- group_celltype_count %>%
    filter(sample_group == grp) %>%
    mutate(celltype = factor(celltype, levels = celltype_levels)) %>%
    arrange(celltype)
  
  labs_grp <- paste0(
    df_grp$celltype, " (",
    round(df_grp$CellNumber / sum(df_grp$CellNumber) * 100, 2), "%)"
  )
  
  pdf(paste0("10.pie_", grp, ".pdf"), width = 12, height = 9)
  pie(
    df_grp$CellNumber,
    labels = labs_grp,
    radius = 1,
    clockwise = FALSE,
    main = paste0(grp, " cell-type composition"),
    col = celltype_colors[as.character(df_grp$celltype)]
  )
  dev.off()
}

cat("AML / NORMAL 细胞比例分析完成。\n")
cat("输出文件包括：\n")
cat("0.cell_metadata_AML_NORMAL.csv\n")
cat("1.group_celltype_count.csv\n")
cat("2.group_celltype_ratio.csv\n")
cat("3.group_celltype_count_wide.csv\n")
cat("4.group_celltype_ratio_wide.csv\n")
cat("5.AML_NORMAL_stacked_ratio.pdf/png\n")
cat("6.AML_NORMAL_stacked_count.pdf/png\n")
cat("7.AML_NORMAL_celltype_ratio_barplot.pdf/png\n")
cat("8.AML_NORMAL_celltype_count_barplot.pdf/png\n")
cat("9.doughnut_NORMAL.pdf\n")
cat("9.doughnut_AML.pdf\n")
cat("10.pie_NORMAL.pdf\n")
cat("10.pie_AML.pdf\n")





# ================================
# 新增结果：X轴为细胞类型，Y轴为细胞比率，
#          每个细胞类型下仅比较 AML 和 NORMAL
# ================================

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 如果还没有 sample_group，则先生成
if (!"sample_group" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$sample_group <- ifelse(
    grepl("^AML", seurat_obj$orig.ident, ignore.case = TRUE),
    "AML",
    "Healthy"
  )
}

seurat_obj$sample_group <- factor(
  seurat_obj$sample_group,
  levels = c("Healthy", "AML")
)

# 如果还没有 celltype，则用 active.ident
if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$celltype <- as.character(Idents(seurat_obj))
}

# 保持 celltype 顺序
if (is.factor(seurat_obj$celltype)) {
  celltype_levels <- levels(seurat_obj$celltype)
} else {
  celltype_levels <- unique(as.character(seurat_obj$celltype))
}
seurat_obj$celltype <- factor(seurat_obj$celltype, levels = celltype_levels)

# 提取元数据
plot_meta <- seurat_obj@meta.data %>%
  dplyr::select(celltype, sample_group)

# 统计 AML / NORMAL 中各细胞类型比例
ratio_df <- plot_meta %>%
  count(sample_group, celltype, name = "CellNumber") %>%
  group_by(sample_group) %>%
  mutate(
    TotalCells = sum(CellNumber),
    Ratio = CellNumber / TotalCells
  ) %>%
  ungroup()

# 查看结果
print(ratio_df)

# 颜色
group_colors <- c(
  "Healthy" = "#6A9BCB",
  "AML"    = "#CB5148"
)

# 绘图：X轴=celltype，Y轴=ratio，分组=AML/NORMAL
p_ratio_compare <- ggplot(
  ratio_df,
  aes(x = celltype, y = Ratio, fill = sample_group)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(
    x = "Cell type",
    y = "Cell ratio",
    fill = "Source",
    title = "Cell ratio of each cell type in AML and Healthy"
  )

# 保存
ggsave(
  filename = "新增_细胞类型比例_AML_vs_NORMAL.pdf",
  plot = p_ratio_compare,
  width = 12,
  height = 7
)

ggsave(
  filename = "新增_细胞类型比例_AML_vs_NORMAL.png",
  plot = p_ratio_compare,
  width = 12,
  height = 7,
  dpi = 300
)

print(p_ratio_compare)





p_ratio_facet <- ggplot(
  ratio_df,
  aes(x = Ratio, y = sample_group, fill = sample_group)
) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  facet_wrap(~ celltype, scales = "free_x") +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Cell ratio",
    y = "Group",
    title = "AML vs Healthy ratio across cell types"
  )

ggsave(
  filename = "新增_分面细胞比例_AML_vs_NORMAL.pdf",
  plot = p_ratio_facet,
  width = 14,
  height = 10
)

ggsave(
  filename = "新增_分面细胞比例_AML_vs_NORMAL.png",
  plot = p_ratio_facet,
  width = 14,
  height = 10,
  dpi = 300
)


library(ggplot2)
library(dplyr)
library(scales)

# ================================
# 1. 准备数据
# ================================

meta <- seurat_obj@meta.data

# 确保有 celltype
meta$celltype <- as.character(Idents(seurat_obj))

# AML / NORMAL 分组
meta$sample_group <- ifelse(
  grepl("^AML", meta$orig.ident, ignore.case = TRUE),
  "AML",
  "Healthy"
)

meta$sample_group <- factor(meta$sample_group, levels = c("Healthy", "AML"))

# ================================
# 2. 计算比例（关键）
# ================================
# 每个 celltype 内，分别在 AML / NORMAL 中，
# 再按 orig.ident 计算比例

plot_df <- meta %>%
  count(celltype, sample_group, orig.ident, name = "CellNumber") %>%
  group_by(celltype, sample_group) %>%
  mutate(
    Total = sum(CellNumber),
    Ratio = CellNumber / Total
  ) %>%
  ungroup()

# ================================
# 3. 作图
# ================================

p_stack_compare <- ggplot(
  plot_df,
  aes(x = celltype, y = Ratio, fill = orig.ident)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.2
  ) +
  
  # 分组：AML vs NORMAL
  facet_grid(~ sample_group) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  
  labs(
    x = "Cell type",
    y = "Cell proportion",
    fill = "Samples",
    title = "Cell composition in AML vs Healthy"
  )

# ================================
# 4. 保存
# ================================
ggsave(
  "Fig_celltype_stack_AML_vs_NORMAL.pdf",
  p_stack_compare,
  width = 14,
  height = 7
)

ggsave(
  "Fig_celltype_stack_AML_vs_NORMAL.png",
  p_stack_compare,
  width = 14,
  height = 7,
  dpi = 300
)

print(p_stack_compare)




# ================================
# AML / NORMAL 在每种细胞类型中的占比
# 横轴：celltype
# 纵轴：比例
# 每个 celltype 只有1根柱子
# 柱子内部：AML / NORMAL 堆叠
# ================================

library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 读取对象
setwd(ANNO_DIR)
seurat_obj <- qread("2.resolution_0.5_annotation.qs")


setwd(RATIO_DIR)
# 使用 celltype
if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$celltype <- as.character(Idents(seurat_obj))
}

# 新增 AML / NORMAL 分组
seurat_obj$sample_group <- ifelse(
  grepl("^AML", seurat_obj$orig.ident, ignore.case = TRUE),
  "AML",
  "Healthy"
)

seurat_obj$sample_group <- factor(
  seurat_obj$sample_group,
  levels = c("Healthy", "AML")
)

# 保持 celltype 顺序
if (is.factor(seurat_obj$celltype)) {
  celltype_levels <- levels(seurat_obj$celltype)
} else {
  celltype_levels <- unique(as.character(seurat_obj$celltype))
}
seurat_obj$celltype <- factor(seurat_obj$celltype, levels = celltype_levels)

# 提取元数据
plot_df <- seurat_obj@meta.data %>%
  dplyr::select(celltype, sample_group) %>%
  count(celltype, sample_group, name = "CellNumber") %>%
  group_by(celltype) %>%
  mutate(
    Total = sum(CellNumber),
    Ratio = CellNumber / Total
  ) %>%
  ungroup()

# 查看结果
print(plot_df)

# 颜色
group_colors <- c(
  "Healthy" = "#6A9BCB",
  "AML"    = "#CB5148"
)

# 绘图：每个celltype只有1根柱子，内部堆叠AML/NORMAL
p <- ggplot(plot_df, aes(x = celltype, y = Ratio, fill = sample_group)) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(
    x = "Cell type",
    y = "Cell proportion",
    fill = "Source",
    title = "Proportion of AML and Healthy cells in each cell type"
    
    )

# 保存
ggsave("AML_NORMAL_ratio_in_each_celltype_1.pdf", p, width = 12, height = 7)
ggsave("AML_NORMAL_ratio_in_each_celltype_1.png", p, width = 12, height = 7, dpi = 300)

print(p)





rm(list = ls())

# =========================================================
# 单细胞细胞比例分析：精简版 + Ro/E 热图 + 堆叠趋势图
# 功能：
# 1）读取注释后的 Seurat 对象
# 2）使用 celltype 作为细胞类型分组
# 3）生成 sample_group：Healthy / AML
# 4）输出细胞数量、细胞比例、Ro/E 结果表
# 5）绘图精简为 3 类：
#    Fig1：总体细胞组成堆叠比例图
#    Fig2：Ro/E 热图
#    Fig3：细胞组成堆叠趋势图
# =========================================================

# ================================
# 0. 加载 R 包
# ================================
suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
})

# ================================
# 1. 路径设置
# ================================
ANNO_DIR  <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
RATIO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\3.Sc_cell_ratio\\final"

if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(RATIO_DIR)) dir.create(RATIO_DIR, recursive = TRUE)

# 注释后 Seurat 对象文件名
SEURAT_FILE <- "2.resolution_0.5_annotation.qs"

# ================================
# 2. 分组参数
# ================================
# 默认：根据 orig.ident 生成 sample_group，用于 AML vs Healthy 比较
# 如果你后续想使用 meta.data 中已有的 group 列，例如 HCC_T / NCHCC_T / HCC_N / NCHCC_N，
# 可将 GROUP_COL 改为 "group"，并相应修改 GROUP_LEVELS 与 GROUP_LABELS。

GROUP_COL <- "sample_group"

GROUP_LEVELS <- c("Healthy", "AML")
GROUP_LABELS <- c("Healthy", "AML")

# 四组示例，暂时不启用：
# GROUP_COL <- "group"
# GROUP_LEVELS <- c("tumor_R", "tumor_NR", "paratumor_R", "paratumor_NR")
# GROUP_LABELS <- c("NCHCC_T", "HCC_T", "NCHCC_N", "HCC_N")

# 趋势图中每个色块内显示百分比的最低阈值
TREND_LABEL_MIN_PERCENT <- 3

# ================================
# 3. 读取 Seurat 对象
# ================================
setwd(ANNO_DIR)
seurat_obj <- qread(SEURAT_FILE)

# 如果 meta.data 中没有 celltype，则使用 active.ident
if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$celltype <- as.character(Idents(seurat_obj))
}

# 保留 celltype 顺序
if (is.factor(seurat_obj@meta.data$celltype)) {
  celltype_levels <- levels(seurat_obj@meta.data$celltype)
} else if (is.factor(Idents(seurat_obj))) {
  celltype_levels <- levels(Idents(seurat_obj))
} else {
  celltype_levels <- unique(as.character(seurat_obj@meta.data$celltype))
}

seurat_obj$celltype <- factor(
  as.character(seurat_obj$celltype),
  levels = celltype_levels
)

# 默认新增 AML / Healthy 分组
# 规则：样本名以 AML 开头 -> AML；其他样本 -> Healthy
seurat_obj$sample_group <- ifelse(
  grepl("^AML", seurat_obj$orig.ident, ignore.case = TRUE),
  "AML",
  "Healthy"
)

seurat_obj$sample_group <- factor(
  seurat_obj$sample_group,
  levels = c("Healthy", "AML")
)

# 检查分组列
if (!GROUP_COL %in% colnames(seurat_obj@meta.data)) {
  stop("Metadata does not contain GROUP_COL: ", GROUP_COL)
}

# ================================
# 4. 提取元数据并统一分组标签
# ================================
meta <- seurat_obj@meta.data %>%
  dplyr::mutate(
    celltype = factor(as.character(celltype), levels = celltype_levels),
    group_raw = as.character(.data[[GROUP_COL]]),
    group = factor(group_raw, levels = GROUP_LEVELS, labels = GROUP_LABELS)
  ) %>%
  dplyr::select(orig.ident, celltype, group_raw, group) %>%
  dplyr::filter(!is.na(celltype), !is.na(group))

# 如果分组映射失败，提示原始值，方便修改 GROUP_LEVELS
unmatched_groups <- setdiff(unique(as.character(seurat_obj@meta.data[[GROUP_COL]])), GROUP_LEVELS)
if (length(unmatched_groups) > 0) {
  message("[Warning] 以下分组值未被 GROUP_LEVELS 使用：")
  print(unmatched_groups)
}

cat("\n细胞类型统计：\n")
print(table(meta$celltype))

cat("\n分组统计：\n")
print(table(meta$group))

cat("\norig.ident 与 group 对应关系：\n")
print(table(meta$orig.ident, meta$group))

setwd(RATIO_DIR)
write.csv(meta, "0.cell_metadata_for_ratio.csv", row.names = FALSE)

# ================================
# 5. 统计细胞数量与比例
# ================================
group_celltype_count <- meta %>%
  dplyr::count(group, celltype, name = "CellNumber") %>%
  tidyr::complete(
    group = factor(GROUP_LABELS, levels = GROUP_LABELS),
    celltype = factor(celltype_levels, levels = celltype_levels),
    fill = list(CellNumber = 0)
  ) %>%
  dplyr::arrange(group, celltype)

group_celltype_ratio <- group_celltype_count %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    TotalCells = sum(CellNumber),
    Proportion = ifelse(TotalCells > 0, CellNumber / TotalCells, NA_real_),
    Percent = Proportion * 100
  ) %>%
  dplyr::ungroup()

count_wide <- group_celltype_count %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = CellNumber,
    values_fill = 0
  )

ratio_wide <- group_celltype_ratio %>%
  dplyr::select(group, celltype, Proportion, Percent) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(Proportion, Percent),
    values_fill = 0
  )

write.csv(group_celltype_count, "1.group_celltype_count.csv", row.names = FALSE)
write.csv(group_celltype_ratio, "2.group_celltype_ratio.csv", row.names = FALSE)
write.csv(count_wide, "3.group_celltype_count_wide.csv", row.names = FALSE)
write.csv(ratio_wide, "4.group_celltype_ratio_wide.csv", row.names = FALSE)

# ================================
# 6. 配色
# ================================
n_celltype <- length(celltype_levels)

if (n_celltype <= 12) {
  celltype_colors <- brewer.pal(max(3, n_celltype), "Paired")[seq_len(n_celltype)]
} else {
  celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltype)
}
names(celltype_colors) <- celltype_levels

group_colors <- c(
  "Healthy" = "#6A9BCB",
  "AML"     = "#CB5148"
)

# 如果 GROUP_LABELS 不是 Healthy / AML，则自动生成分组颜色
if (!all(GROUP_LABELS %in% names(group_colors))) {
  group_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(GROUP_LABELS)),
    GROUP_LABELS
  )
}

# ================================
# 7. Fig1：总体细胞组成堆叠比例图
#    原脚本中的堆叠比例图、堆叠数量图、分组柱状图、饼图、甜甜圈图等，
#    这里统一精简为这一种基础细胞比例展示方式。
# ================================
p_stacked <- ggplot(
  group_celltype_ratio,
  aes(x = group, y = Proportion, fill = celltype)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors, breaks = celltype_levels) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Cell-type composition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("5.Fig1_celltype_stacked_proportion.pdf", p_stacked, width = 7, height = 5.5)
ggsave("5.Fig1_celltype_stacked_proportion.png", p_stacked, width = 7, height = 5.5, dpi = 300)

# ================================
# 8. Fig2：Ro/E 热图
#    Ro/E = observed / expected
#    observed：实际细胞数
#    expected：基于列联表独立性假设得到的期望细胞数
#    为了让 Ro/E = 1 对应颜色中心，热图颜色使用 log2(Ro/E)，显示内容为符号。
# ================================
tab <- as.matrix(table(meta$celltype, meta$group))
chisq_res <- suppressWarnings(chisq.test(tab))
expected <- chisq_res$expected
roe_values <- tab / expected

map_roe_to_symbol <- function(x) {
  ifelse(
    x > 1, "+++",
    ifelse(
      x > 0.8 & x <= 1, "++",
      ifelse(
        x > 0.2 & x <= 0.8, "+",
        ifelse(x > 0 & x <= 0.2, "+/-", "")
      )
    )
  )
}

symbol_matrix <- apply(roe_values, c(1, 2), map_roe_to_symbol)

roe_table <- as.data.frame(as.table(roe_values)) %>%
  dplyr::rename(celltype = Var1, group = Var2, RoE = Freq)

expected_table <- as.data.frame(as.table(expected)) %>%
  dplyr::rename(celltype = Var1, group = Var2, Expected = Freq)

observed_table <- as.data.frame(as.table(tab)) %>%
  dplyr::rename(celltype = Var1, group = Var2, Observed = Freq)

roe_out <- observed_table %>%
  dplyr::left_join(expected_table, by = c("celltype", "group")) %>%
  dplyr::left_join(roe_table, by = c("celltype", "group")) %>%
  dplyr::mutate(RoE_symbol = map_roe_to_symbol(RoE))

write.csv(observed_table, "6.RoE_observed_count.csv", row.names = FALSE)
write.csv(expected_table, "7.RoE_expected_count.csv", row.names = FALSE)
write.csv(roe_out, "8.RoE_result_long.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(roe_values), "9.RoE_matrix.csv")

roe_log2 <- log2(roe_values + 1e-6)
roe_cap <- max(1, as.numeric(stats::quantile(abs(roe_log2), probs = 0.98, na.rm = TRUE)))
roe_plot <- pmax(pmin(roe_log2, roe_cap), -roe_cap)

roe_breaks <- seq(-roe_cap, roe_cap, length.out = 101)
roe_colors <- colorRampPalette(c("#053061", "white", "#B2182B"))(100)

pheatmap(
  roe_plot,
  display_numbers = symbol_matrix,
  color = roe_colors,
  breaks = roe_breaks,
  main = "Ro/E heatmap",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 11,
  fontsize_number = 12,
  border_color = "grey70",
  angle_col = 0,
  filename = "10.Fig2_RoE_heatmap_symbol.pdf",
  width = 7,
  height = max(4.5, 0.28 * length(celltype_levels) + 2)
)

pheatmap(
  roe_plot,
  display_numbers = symbol_matrix,
  color = roe_colors,
  breaks = roe_breaks,
  main = "Ro/E heatmap",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 11,
  fontsize_number = 12,
  border_color = "grey70",
  angle_col = 0,
  filename = "10.Fig2_RoE_heatmap_symbol.png",
  width = 7,
  height = max(4.5, 0.28 * length(celltype_levels) + 2)
)

# ================================
# 9. Fig3：细胞比例堆叠趋势图
#    每个分组是一根 100% 堆叠柱；相邻分组之间用半透明多边形连接同一种 celltype。
#    默认用于 Healthy -> AML；如果 GROUP_LABELS 为 4 组，则会自动连接第1-2、2-3、3-4组。
# ================================
stack_order <- rev(celltype_levels)
bar_half_width <- 0.28
trend_df <- group_celltype_ratio %>%
  dplyr::mutate(
    group = factor(as.character(group), levels = GROUP_LABELS),
    celltype = factor(as.character(celltype), levels = stack_order)
  ) %>%
  dplyr::arrange(group, celltype) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    ymin = cumsum(dplyr::lag(Percent, default = 0)),
    ymax = cumsum(Percent),
    label_y = (ymin + ymax) / 2,
    x = as.numeric(group)
  ) %>%
  dplyr::ungroup()

make_polygon_df <- function(df, group_levels, bar_half_width = 0.28) {
  if (length(group_levels) < 2) return(NULL)
  polygon_list <- vector("list", length(group_levels) - 1)
  
  for (i in seq_len(length(group_levels) - 1)) {
    g1 <- group_levels[i]
    g2 <- group_levels[i + 1]
    
    polygon_list[[i]] <- df %>%
      dplyr::filter(group %in% c(g1, g2)) %>%
      dplyr::group_by(celltype) %>%
      dplyr::reframe(
        x = c(i + bar_half_width, i + 1 - bar_half_width,
              i + 1 - bar_half_width, i + bar_half_width),
        y = c(
          ymax[group == g1][1],
          ymax[group == g2][1],
          ymin[group == g2][1],
          ymin[group == g1][1]
        ),
        point_id = 1:4,
        pair = paste0(g1, "_vs_", g2)
      )
  }
  
  dplyr::bind_rows(polygon_list)
}

polygon_df <- make_polygon_df(
  df = trend_df,
  group_levels = GROUP_LABELS,
  bar_half_width = bar_half_width
)

p_trend <- ggplot() +
  geom_polygon(
    data = polygon_df,
    aes(
      x = x,
      y = y,
      group = interaction(celltype, pair),
      fill = celltype
    ),
    alpha = 0.25,
    color = "grey45",
    linewidth = 0.25
  ) +
  geom_rect(
    data = trend_df,
    aes(
      xmin = x - bar_half_width,
      xmax = x + bar_half_width,
      ymin = ymin,
      ymax = ymax,
      fill = celltype
    ),
    color = "grey25",
    linewidth = 0.35,
    alpha = 0.90
  ) +
  geom_text(
    data = trend_df %>% dplyr::filter(Percent >= TREND_LABEL_MIN_PERCENT),
    aes(
      x = x,
      y = label_y,
      label = sprintf("%.1f%%", Percent)
    ),
    size = 3.4,
    color = "black",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = celltype_colors,
    breaks = celltype_levels
  ) +
  scale_x_continuous(
    breaks = seq_along(GROUP_LABELS),
    labels = GROUP_LABELS,
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Cell type",
    title = "Cell-type proportion trend"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      size = 11,
      face = "bold"
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.title.y = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.35, "cm"),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "11.Fig3_celltype_proportion_stacked_trend.pdf",
  plot = p_trend,
  width = max(6.5, 1.4 * length(GROUP_LABELS) + 3.5),
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = "11.Fig3_celltype_proportion_stacked_trend.png",
  plot = p_trend,
  width = max(6.5, 1.4 * length(GROUP_LABELS) + 3.5),
  height = 5.5,
  dpi = 300
)

print(p_stacked)
print(p_trend)

cat("\n细胞比例分析完成。\n")
cat("\n输出表格：\n")
cat("0.cell_metadata_for_ratio.csv\n")
cat("1.group_celltype_count.csv\n")
cat("2.group_celltype_ratio.csv\n")
cat("3.group_celltype_count_wide.csv\n")
cat("4.group_celltype_ratio_wide.csv\n")
cat("6.RoE_observed_count.csv\n")
cat("7.RoE_expected_count.csv\n")
cat("8.RoE_result_long.csv\n")
cat("9.RoE_matrix.csv\n")

cat("\n输出图片：\n")
cat("5.Fig1_celltype_stacked_proportion.pdf/png\n")
cat("10.Fig2_RoE_heatmap_symbol.pdf/png\n")
cat("11.Fig3_celltype_proportion_stacked_trend.pdf/png\n")













# END ORIGINAL SCRIPT P12

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P13: 08_SC_cell_ratio/plot_cell_ratio_pub.R
# ==============================================================================
############################################################
## 02_plot_cell_proportion_pub_core.R
## Core publication-style cell proportion figures
## Matched to the final annotation UMAP style/palette
############################################################

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(pheatmap)
  library(grid)
  library(forcats)
})

source("00_pub_theme_colors.R")

############################################################
## Compatibility function: standardize disease/group labels
############################################################
if (!exists("standardize_group")) {
  standardize_group <- function(x) {
    x <- as.character(x)
    
    dplyr::case_when(
      grepl("^AML", x, ignore.case = TRUE) ~ "AML",
      grepl("healthy|normal|control|^BM|^NBM", x, ignore.case = TRUE) ~ "Healthy",
      x %in% c("AML", "Healthy") ~ x,
      x %in% c("Normal", "NORMAL") ~ "Healthy",
      TRUE ~ x
    )
  }
}

############################################################
## 0. Basic settings
############################################################
OUT_DIR <- RATIO_DIR
safe_mkdir(OUT_DIR)

SEURAT_FILE  <- SEURAT_ANNO_FILE
SAMPLE_COL   <- "orig.ident"
CELLTYPE_COL <- "celltype"
GROUP_CANDIDATES <- c("type", "sample_group", "group")
TREND_LABEL_MIN_PERCENT <- 3

############################################################
## 1. Palette matched to your final UMAP
############################################################
celltype_palette_short <- c(
  "HSC-like"      = "#1F4E79",
  "GranuleProg"   = "#2F75B5",
  "Promyelo"      = "#5B9BD5",
  "MyeloidProg"   = "#7EC0EE",
  "Progenitor"    = "#4A86CF",
  "CyclingProg"   = "#7B68B0",
  "GMP-like"      = "#18B7C9",
  "LymphoidProg"  = "#68C7A7",
  "cMono"         = "#20B07A",
  "intMono"       = "#4CB648",
  "ncMono"        = "#9CCD4A",
  "Macrophage"    = "#9B772B",
  "cDC"           = "#1EA5C0",
  "moDC"          = "#47B4CF",
  "pDC"           = "#00A38B",
  "NaiveT"        = "#6F4AA2",
  "NK_T"          = "#9564B7",
  "Bcell"         = "#B247A3",
  "Plasma"        = "#E62B8A",
  "Mast"          = "#E6AA00",
  "Erythrocyte"   = "#D25759"
)

preferred_celltype_order <- c(
  "HSC-like", "GranuleProg", "Promyelo", "MyeloidProg", "Progenitor",
  "CyclingProg", "GMP-like", "LymphoidProg", "cMono", "intMono",
  "ncMono", "Macrophage", "cDC", "moDC", "pDC", "NaiveT",
  "NK_T", "Bcell", "Plasma", "Mast", "Erythrocyte"
)

disease_cols <- c(
  "Healthy" = "#8DB6E2",
  "AML"     = "#E58C8A"
)

############################################################
## 2. Read object and prepare metadata
############################################################
if (!file.exists(SEURAT_FILE)) {
  stop("Cannot find annotated Seurat object: ", SEURAT_FILE)
}

seurat_obj <- qread(SEURAT_FILE)

if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[[CELLTYPE_COL]] <- as.character(Idents(seurat_obj))
}
if (!"celltype_short" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj <- add_celltype_short(seurat_obj, celltype_col = CELLTYPE_COL, short_col = "celltype_short")
}

group_col <- GROUP_CANDIDATES[GROUP_CANDIDATES %in% colnames(seurat_obj@meta.data)][1]
if (is.na(group_col)) {
  seurat_obj$type <- ifelse(grepl("^AML", seurat_obj@meta.data[[SAMPLE_COL]], ignore.case = TRUE), "AML", "Healthy")
  group_col <- "type"
}

meta <- seurat_obj@meta.data %>%
  mutate(
    sample_id = as.character(.data[[SAMPLE_COL]]),
    celltype = as.character(.data[["celltype_short"]]),
    disease  = standardize_group(.data[[group_col]])
  ) %>%
  filter(!is.na(sample_id), sample_id != "",
         !is.na(celltype), celltype != "",
         !is.na(disease), disease != "")

celltype_levels <- c(
  intersect(preferred_celltype_order, unique(meta$celltype)),
  setdiff(unique(meta$celltype), preferred_celltype_order)
)
sample_levels <- unique(meta$sample_id)
disease_levels <- unique(meta$disease)
if (all(c("Healthy", "AML") %in% disease_levels)) {
  disease_levels <- c("Healthy", "AML", setdiff(disease_levels, c("Healthy", "AML")))
}

meta <- meta %>%
  mutate(
    celltype = factor(celltype, levels = celltype_levels),
    sample_id = factor(sample_id, levels = sample_levels),
    disease = factor(disease, levels = disease_levels)
  )

cell_cols <- get_named_colors(celltype_levels, palette_named = celltype_palette_short)
dis_cols  <- get_named_colors(disease_levels, palette_named = c(disease_cols, group_colors))

write.csv(meta, file.path(OUT_DIR, "00_core_cell_metadata.csv"), row.names = FALSE)

############################################################
## 3. Summary tables needed for plotting
############################################################
# sample x celltype
sample_celltype_df <- meta %>%
  count(sample_id, disease, celltype, name = "CellNumber") %>%
  group_by(sample_id) %>%
  mutate(
    SampleTotal = sum(CellNumber),
    SamplePercent = ifelse(SampleTotal > 0, CellNumber / SampleTotal * 100, 0)
  ) %>%
  ungroup()

# disease x celltype
Disease_celltype_df <- meta %>%
  count(disease, celltype, name = "CellNumber") %>%
  group_by(disease) %>%
  mutate(
    DiseaseTotal = sum(CellNumber),
    DiseasePercent = ifelse(DiseaseTotal > 0, CellNumber / DiseaseTotal * 100, 0)
  ) %>%
  ungroup()

# within each celltype, contribution from different sample sources
within_celltype_sample_df <- meta %>%
  count(celltype, sample_id, name = "CellNumber") %>%
  group_by(celltype) %>%
  mutate(
    CelltypeTotal = sum(CellNumber),
    WithinCelltypePercent = ifelse(CelltypeTotal > 0, CellNumber / CelltypeTotal * 100, 0)
  ) %>%
  ungroup()

# total cells per sample
sample_total_df <- meta %>%
  count(sample_id, disease, name = "TotalCells")

write.csv(sample_celltype_df, file.path(OUT_DIR, "01_sample_celltype_long.csv"), row.names = FALSE)
write.csv(Disease_celltype_df, file.path(OUT_DIR, "02_disease_celltype_long.csv"), row.names = FALSE)
write.csv(within_celltype_sample_df, file.path(OUT_DIR, "03_within_celltype_sample_long.csv"), row.names = FALSE)
write.csv(sample_total_df, file.path(OUT_DIR, "04_sample_total_cells.csv"), row.names = FALSE)

############################################################
## 4. Shared plot theme
############################################################
prop_theme <- function(base_size = 12) {
  theme_sc_pub(base_size = base_size) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_blank()
    )
}

############################################################
## 5. 1. 横轴为样本纵轴为细胞绝对数量
############################################################
p_sample_total <- ggplot(sample_total_df, aes(x = sample_id, y = TotalCells, fill = disease)) +
  geom_col(width = 0.76, color = "black", linewidth = 0.16) +
  scale_fill_manual(values = dis_cols, drop = FALSE) +
  labs(x = "Sample", y = "Absolute cell count", fill = "Disease status") +
  prop_theme(12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 9)
  )

save_plot_pub(
  p_sample_total,
  file.path(OUT_DIR, "01_横轴为样本纵轴为细胞绝对数量"),
  width = max(7.2, 0.42 * length(sample_levels) + 4.0),
  height = 5.3
)

############################################################
## 6. 2. 横轴为细胞类型纵轴为不同样本来源细胞比例
############################################################
p_celltype_sample_percent <- ggplot(sample_celltype_df, aes(x = celltype, y = SamplePercent, fill = sample_id)) +
  geom_col(position = position_dodge(width = 0.80), width = 0.72,
           color = "black", linewidth = 0.15) +
  scale_y_continuous(
    labels = label_percent(scale = 1, accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = "Cell type", y = "Cell proportion (%)", fill = "Sample") +
  prop_theme(12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 8.2)
  )

save_plot_pub(
  p_celltype_sample_percent,
  file.path(OUT_DIR, "02_横轴为细胞类型纵轴为不同样本来源细胞比例"),
  width = max(9.0, 0.36 * length(celltype_levels) + 5.0),
  height = 5.8
)

############################################################
## 7. 3. 横轴为细胞类型纵轴为不同样本来源细胞绝对数量
############################################################
p_celltype_sample_count <- ggplot(sample_celltype_df, aes(x = celltype, y = CellNumber, fill = sample_id)) +
  geom_col(position = position_dodge(width = 0.80), width = 0.72,
           color = "black", linewidth = 0.15) +
  labs(x = "Cell type", y = "Absolute cell count", fill = "Sample") +
  prop_theme(12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 8.2)
  )

save_plot_pub(
  p_celltype_sample_count,
  file.path(OUT_DIR, "03_横轴为细胞类型纵轴为不同样本来源细胞绝对数量"),
  width = max(9.0, 0.36 * length(celltype_levels) + 5.0),
  height = 5.8
)

############################################################
## 8. 4. 横轴为不同样本来源细胞比例纵轴为细胞类型
############################################################
p_horizontal_prop <- ggplot(within_celltype_sample_df,
                            aes(x = WithinCelltypePercent / 100, y = celltype, fill = sample_id)) +
  geom_col(width = 0.76, color = "black", linewidth = 0.14) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Proportion of different sample sources", y = "Cell type", fill = "Sample") +
  prop_theme(12) +
  theme(
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 8.2)
  )

save_plot_pub(
  p_horizontal_prop,
  file.path(OUT_DIR, "04_横轴为不同样本来源细胞比例纵轴为细胞类型"),
  width = max(8.6, 0.36 * length(sample_levels) + 3.8),
  height = max(6.4, 0.24 * length(celltype_levels) + 2.2)
)

############################################################
## 9. 5. Fig1_celltype_stacked_proportion（重点）
############################################################
p_fig1 <- ggplot(Disease_celltype_df, aes(x = disease, y = DiseasePercent / 100, fill = celltype)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.18) +
  scale_fill_manual(values = cell_cols, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Cell proportion", fill = "Cell type") +
  prop_theme(12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.32, "cm"),
    legend.text = element_text(size = 8.4)
  )

save_plot_pub(
  p_fig1,
  file.path(OUT_DIR, "Fig1_celltype_stacked_proportion"),
  width = max(6.3, 1.0 * length(disease_levels) + 4.2),
  height = 5.5
)

############################################################
## 10. 6. Fig2_Ro/e_heatmap_symbol（重点）
############################################################
observed <- as.matrix(table(meta$celltype, meta$disease))
chisq_res <- suppressWarnings(chisq.test(observed))
expected <- chisq_res$expected
roe <- observed / expected

map_roe_to_symbol <- function(x) {
  ifelse(x > 1, "+++",
         ifelse(x > 0.8 & x <= 1, "++",
                ifelse(x > 0.2 & x <= 0.8, "+",
                       ifelse(x > 0 & x <= 0.2, "+/-", ""))))
}

symbol_mat <- apply(roe, c(1, 2), map_roe_to_symbol)
roe_long <- as.data.frame(as.table(roe))
colnames(roe_long) <- c("celltype", "disease", "RoE")
roe_long$RoE_symbol <- map_roe_to_symbol(roe_long$RoE)
write.csv(roe_long, file.path(OUT_DIR, "Fig2_RoE_long_table.csv"), row.names = FALSE)

roe_log2 <- log2(roe + 1e-6)
roe_cap <- max(1, as.numeric(quantile(abs(roe_log2), 0.98, na.rm = TRUE)))
roe_plot <- pmax(pmin(roe_log2, roe_cap), -roe_cap)
roe_breaks <- seq(-roe_cap, roe_cap, length.out = 101)
heat_h <- max(4.8, 0.28 * nrow(roe_plot) + 2.0)

pheatmap(
  roe_plot,
  display_numbers = symbol_mat,
  color = roe_colors,
  breaks = roe_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80",
  fontsize = 10.5,
  fontsize_number = 12,
  angle_col = 0,
  main = "Ro/E heatmap",
  cellwidth = 60,
  cellheight = 22,
  filename = file.path(OUT_DIR, "Fig2_RoE_heatmap_symbol.pdf"),
  width = 4,
  height = heat_h
)

pheatmap(
  roe_plot,
  display_numbers = symbol_mat,
  color = roe_colors,
  breaks = roe_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80",
  fontsize = 10.5,
  fontsize_number = 12,
  angle_col = 0,
  main = "Ro/E heatmap",
  cellwidth = 60,
  cellheight = 22,
  filename = file.path(OUT_DIR, "Fig2_RoE_heatmap_symbol.png"),
  width = 4,
  height = heat_h
)

############################################################
## 11. 7. Fig3_celltype_proportion_stacked_trend（重点，窄长）
############################################################
trend_levels <- disease_levels
stack_order <- rev(celltype_levels)
bar_half_width <- 0.30

trend_df <- Disease_celltype_df %>%
  mutate(
    disease = factor(as.character(disease), levels = trend_levels),
    celltype = factor(as.character(celltype), levels = stack_order)
  ) %>%
  arrange(disease, celltype) %>%
  group_by(disease) %>%
  mutate(
    ymin = cumsum(lag(DiseasePercent, default = 0)),
    ymax = cumsum(DiseasePercent),
    label_y = (ymin + ymax) / 2,
    x = as.numeric(disease)
  ) %>%
  ungroup()

polygon_df <- if (length(trend_levels) >= 2) {
  lapply(seq_len(length(trend_levels) - 1), function(i) {
    g1 <- trend_levels[i]
    g2 <- trend_levels[i + 1]
    trend_df %>%
      filter(disease %in% c(g1, g2)) %>%
      group_by(celltype) %>%
      reframe(
        x = c(i + bar_half_width, i + 1 - bar_half_width, i + 1 - bar_half_width, i + bar_half_width),
        y = c(ymax[disease == g1][1], ymax[disease == g2][1],
              ymin[disease == g2][1], ymin[disease == g1][1]),
        pair = paste0(g1, "_vs_", g2)
      )
  }) %>% bind_rows()
} else {
  data.frame()
}

p_trend <- ggplot() +
  {if (nrow(polygon_df) > 0) geom_polygon(
    data = polygon_df,
    aes(x = x, y = y, group = interaction(celltype, pair), fill = celltype),
    alpha = 0.25,
    color = "grey50",
    linewidth = 0.18
  )} +
  geom_rect(
    data = trend_df,
    aes(xmin = x - bar_half_width, xmax = x + bar_half_width, ymin = ymin, ymax = ymax, fill = celltype),
    color = "grey25",
    linewidth = 0.28,
    alpha = 0.94
  ) +
  geom_text(
    data = trend_df %>% filter(DiseasePercent >= TREND_LABEL_MIN_PERCENT),
    aes(x = x, y = label_y, label = sprintf("%.1f%%", DiseasePercent)),
    size = 3.0,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = cell_cols, drop = FALSE) +
  guides(fill = "none") +
  scale_x_continuous(
    breaks = seq_along(trend_levels),
    labels = trend_levels,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell proportion (%)", fill = "Cell type") +
  theme_sc_pub(base_size = 12) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.key.size = unit(0.30, "cm"),
    legend.text = element_text(size = 8.2),
    legend.title = element_text(size = 10.2, face = "bold")
  )

save_plot_pub(
  p_trend,
  file.path(OUT_DIR, "Fig3_celltype_proportion_stacked_trend"),
  width = max(3.6, 0.75 * length(trend_levels) + 1.6),
  height = 7.0
)

message("Core cell proportion figures were exported to: ", OUT_DIR)

# END ORIGINAL SCRIPT P13

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P14: 09_SC_CYTOTRACE2.R
# ==============================================================================
####CytoTRACE2单细胞干性分析####

####旧版本，新版在下面####

#初始化 ----
rm(list = ls())

library(Seurat)
library(CytoTRACE2)
library(qs)
library(ggplot2)
# 根目录(按需修改)
ANNO_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/2.Sc_annotation/new"
TRACE_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/12.CytoTRACE2/0.new"


if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(TRACE_DIR)) dir.create(TRACE_DIR, recursive = TRUE)
# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

setwd(ANNO_DIR)
# 读取Seurat对象
seurat_obj <- qread("2.resolution_0.5_annotation.qs")


# 设置工作目录
setwd(TRACE_DIR)
cytotrace2_result <- cytotrace2(seurat_obj, is_seurat = TRUE,
                                slot_type = "counts",
                                species = "human",
                                ncores = 5)


phe.anno = as.data.frame(cytotrace2_result$celltype)
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  expression_data = NULL, 
                  is_seurat = TRUE,
                  annotation = phe.anno, 
                  slot_type = "counts",
                  pc_dims = 50)


p2 <- plots$CytoTRACE2_Boxplot_byPheno + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
  )
pdf("CytoTRACE2_Boxplot.pdf", width = 8, height = 5)
p2
dev.off()

library(ggplot2)
library(ggsignif)
library(patchwork)
library(dplyr)
# 发表级统一主题
pub_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    legend.position = "none"
  )

# 分化潜能 UMAP 图
p_umap <- plots$CytoTRACE2_Potency_UMAP +
  ggtitle("CytoTRACE2 Differentiation Potency") +
  pub_theme +
  theme(legend.position = "right")

#细胞类型分组箱线图-不带显著性
# p_box <- plots$CytoTRACE2_Boxplot_byPheno +
#   ggtitle("CytoTRACE2 Score by Cell Type") +
#   pub_theme +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   ) +
#   ylab("CytoTRACE2 Potency Score") + xlab("Cell Type")

# ========== 绘制带显著性星号的箱线图 ==========
p_box <- plots$CytoTRACE2_Boxplot_byPheno +
  ggtitle("CytoTRACE2 Potency Score by Cell Type") +
  xlab("Cell Type") +
  ylab("CytoTRACE2 Potency Score") +
  pub_theme +
  
  # ✨ 自动添加显著性星号（核心代码）
  geom_signif(
    test = "wilcox.test",        # 常用非参数检验
    map_signif_level = TRUE,     # 自动显示 * / ** / ***
    textsize = 3.8,              # 星号文字大小
    vjust = 0.4                  # 星号位置微调
  )


#基因特征 UMAP 图 如果有
# if(!is.null(plots$CytoTRACE2_Gene_UMAP)){
#   p_gene <- plots$CytoTRACE2_Gene_UMAP +
#     ggtitle("CytoTRACE2 Gene Signature") +
#     pub_theme
#   print(p_gene)
# }
#密度分布图
if(!is.null(plots$CytoTRACE2_Density)){
  p_density <- plots$CytoTRACE2_Density +
    ggtitle("Potency Score Distribution") +
    pub_theme
  print(p_density)
}
ggsave("CytoTRACE2_UMAP.png", p_umap, width = 8, height = 6, dpi = 300)
ggsave("CytoTRACE2_Boxplot.pdf", p_box, width = 10, height = 6, dpi = 300)









############################################
## CytoTRACE2 单细胞干性分析全流程（最终总修正版）
## 已修复：
## 1) CytoTRACE2 score 提取与NA处理
## 2) cytotrace2_result 为 Seurat object 时的 Potency category 提取
## 3) UMAP / tSNE 离散分类绘图
## 4) 单细胞层面与样本层面统计 left_join 报错
############################################

rm(list = ls())
options(stringsAsFactors = FALSE)

#############################
## 0. 加载依赖包
#############################
required_pkgs <- c(
  "Seurat",
  "CytoTRACE2",
  "qs",
  "ggplot2",
  "dplyr",
  "tibble",
  "patchwork",
  "rstatix",
  "ggpubr",
  "viridis"
)

need_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(need_install) > 0) {
  install.packages(need_install)
}

library(Seurat)
library(CytoTRACE2)
library(qs)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(rstatix)
library(ggpubr)
library(viridis)

#############################
## 1. 参数设置
#############################
set.seed(1234)

ANNO_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/2.Sc_annotation/new"
TRACE_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/12.CytoTRACE2/0.new"

SEURAT_FILE <- file.path(ANNO_DIR, "2.resolution_0.5_annotation.qs")

CYTO_QS_FILE      <- file.path(TRACE_DIR, "cytotrace2_result.qs")
CYTO_RDS_FILE     <- file.path(TRACE_DIR, "cytotrace2_result.rds")
PLOTS_RDS_FILE    <- file.path(TRACE_DIR, "cytotrace2_plots.rds")
CELL_META_FILE    <- file.path(TRACE_DIR, "CytoTRACE2_cell_metadata.csv")
SAMPLE_SUM_FILE   <- file.path(TRACE_DIR, "CytoTRACE2_sample_summary.csv")
CELL_STATS_FILE   <- file.path(TRACE_DIR, "CytoTRACE2_cell_level_stats.csv")
SAMPLE_STATS_FILE <- file.path(TRACE_DIR, "CytoTRACE2_sample_level_stats.csv")
SESSIONINFO_FILE  <- file.path(TRACE_DIR, "sessionInfo.txt")

SPECIES   <- "human"
NCORES    <- 5
SLOT_TYPE <- "counts"
PC_DIMS   <- 50
TSNE_DIMS <- 1:min(30, PC_DIMS)

# 按你的对象实际列名设置
CELLTYPE_COL <- "celltype"
SAMPLE_COL   <- "orig.ident"
GROUP_COL    <- "type"

MIN_CELLS_PER_GROUP   <- 10
MIN_SAMPLES_PER_GROUP <- 2

#############################
## 2. 工具函数
#############################
safe_mkdir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
}

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

check_file_exists <- function(x) {
  if (!file.exists(x)) stop("文件不存在: ", x)
}

save_pdf <- function(file, plot, w = 8, h = 6) {
  ggsave(filename = file, plot = plot, width = w, height = h, device = cairo_pdf)
}

pub_theme <- theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8, color = "black"),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

check_seurat_object <- function(seurat_obj, celltype_col, sample_col = NULL, group_col = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("输入对象不是Seurat对象。")
  }
  
  if (!"RNA" %in% Assays(seurat_obj)) {
    stop("Seurat对象中未找到RNA assay。")
  }
  
  if (ncol(seurat_obj) == 0) {
    stop("Seurat对象中没有细胞。")
  }
  
  meta_cols <- colnames(seurat_obj@meta.data)
  
  if (!celltype_col %in% meta_cols) {
    stop("meta.data中缺少细胞类型列: ", celltype_col)
  }
  
  if (!is.null(sample_col) && !sample_col %in% meta_cols) {
    warning("meta.data中未找到样本列: ", sample_col, "，将跳过样本层面统计。")
  }
  
  if (!is.null(group_col) && !group_col %in% meta_cols) {
    warning("meta.data中未找到分组列: ", group_col, "，将跳过group相关统计。")
  }
  
  invisible(TRUE)
}

extract_cytotrace2_score <- function(cytotrace2_result, seurat_obj) {
  possible_slots <- c(
    "CytoTRACE2_Score",
    "cytotrace2_score",
    "score",
    "CytoTRACE2",
    "preKNN_CytoTRACE2"
  )
  
  out <- NULL
  
  # 如果是 Seurat 对象，优先从 meta.data / misc 中找
  if (inherits(cytotrace2_result, "Seurat")) {
    meta_cols <- colnames(cytotrace2_result@meta.data)
    hit_meta <- intersect(possible_slots, meta_cols)
    
    if (length(hit_meta) > 0) {
      nm <- hit_meta[1]
      msg("检测到CytoTRACE2 score字段（meta.data）: ", nm)
      out <- cytotrace2_result@meta.data[[nm]]
    }
    
    if (is.null(out) && !is.null(cytotrace2_result@misc) && length(cytotrace2_result@misc) > 0) {
      misc_names <- names(cytotrace2_result@misc)
      hit_misc <- intersect(possible_slots, misc_names)
      
      if (length(hit_misc) > 0) {
        nm <- hit_misc[1]
        obj <- cytotrace2_result@misc[[nm]]
        msg("检测到CytoTRACE2 score字段（misc）: ", nm)
        
        if (is.numeric(obj) && length(obj) == ncol(seurat_obj)) {
          out <- obj
        } else if (is.data.frame(obj) && nrow(obj) == ncol(seurat_obj)) {
          num_cols <- sapply(obj, is.numeric)
          if (any(num_cols)) out <- obj[[which(num_cols)[1]]]
        } else if (is.matrix(obj) && nrow(obj) == ncol(seurat_obj) && is.numeric(obj)) {
          out <- obj[, 1]
        }
      }
    }
  } else {
    # 普通 list / 对象
    for (nm in possible_slots) {
      obj <- tryCatch(cytotrace2_result[[nm]], error = function(e) NULL)
      if (!is.null(obj)) {
        msg("检测到CytoTRACE2 score字段: ", nm)
        
        if (is.numeric(obj) && length(obj) == ncol(seurat_obj)) {
          out <- obj
          break
        }
        
        if (is.data.frame(obj) && nrow(obj) == ncol(seurat_obj)) {
          num_cols <- sapply(obj, is.numeric)
          if (any(num_cols)) {
            out <- obj[[which(num_cols)[1]]]
            break
          }
        }
        
        if (is.matrix(obj) && nrow(obj) == ncol(seurat_obj) && is.numeric(obj)) {
          out <- obj[, 1]
          break
        }
      }
    }
  }
  
  if (is.null(out)) {
    stop("无法自动提取CytoTRACE2 score，请检查 cytotrace2_result 的结构。")
  }
  
  out <- as.numeric(out)
  
  if (length(out) != ncol(seurat_obj)) {
    stop("提取出的CytoTRACE2 score长度与Seurat对象细胞数不一致。")
  }
  
  return(out)
}

extract_potency_category <- function(cytotrace2_result, seurat_obj, cyto_score = NULL) {
  possible_slots <- c(
    "potency_category",
    "Potency_Category",
    "CytoTRACE2_Category",
    "CytoTRACE_Category",
    "label",
    "labels",
    "category"
  )
  
  out <- NULL
  ncell <- ncol(seurat_obj)
  
  # 1. 如果 cytotrace2_result 是 Seurat 对象
  if (inherits(cytotrace2_result, "Seurat")) {
    msg("cytotrace2_result 为 Seurat 对象，优先从 meta.data 中查找 Potency category ...")
    
    meta_cols <- colnames(cytotrace2_result@meta.data)
    hit_cols <- intersect(possible_slots, meta_cols)
    
    if (length(hit_cols) > 0) {
      nm <- hit_cols[1]
      msg("在 meta.data 中检测到 Potency category 字段: ", nm)
      out <- as.character(cytotrace2_result@meta.data[[nm]])
    }
    
    if (is.null(out) && !is.null(cytotrace2_result@misc) && length(cytotrace2_result@misc) > 0) {
      misc_names <- names(cytotrace2_result@misc)
      hit_misc <- intersect(possible_slots, misc_names)
      
      if (length(hit_misc) > 0) {
        nm <- hit_misc[1]
        msg("在 misc 中检测到 Potency category 字段: ", nm)
        obj <- cytotrace2_result@misc[[nm]]
        
        if (is.vector(obj) && length(obj) == ncell) {
          out <- as.character(obj)
        } else if (is.data.frame(obj) && nrow(obj) == ncell) {
          char_cols <- sapply(obj, function(x) is.character(x) || is.factor(x))
          if (any(char_cols)) {
            out <- as.character(obj[[which(char_cols)[1]]])
          }
        }
      }
    }
  } else {
    # 2. 普通 list / 对象
    msg("cytotrace2_result 不是 Seurat 对象，按普通对象查找 Potency category ...")
    obj_names <- names(cytotrace2_result)
    hit_cols <- intersect(possible_slots, obj_names)
    
    if (length(hit_cols) > 0) {
      nm <- hit_cols[1]
      msg("检测到 Potency category 字段: ", nm)
      obj <- tryCatch(cytotrace2_result[[nm]], error = function(e) NULL)
      
      if (!is.null(obj)) {
        if (is.vector(obj) && length(obj) == ncell) {
          out <- as.character(obj)
        } else if (is.data.frame(obj) && nrow(obj) == ncell) {
          char_cols <- sapply(obj, function(x) is.character(x) || is.factor(x))
          if (any(char_cols)) {
            out <- as.character(obj[[which(char_cols)[1]]])
          }
        }
      }
    }
  }
  
  # 3. 如果没找到，按 score 自动分组
  if (is.null(out)) {
    msg("未直接检测到 Potency category 字段，按 CytoTRACE2 score 自动分组...")
    
    if (is.null(cyto_score)) {
      stop("未检测到 Potency category 字段，且未提供 cyto_score 用于自动分组。")
    }
    
    score_vec <- as.numeric(cyto_score)
    score_vec[is.infinite(score_vec)] <- NA
    score_vec[is.nan(score_vec)] <- NA
    
    if (all(is.na(score_vec))) {
      stop("cyto_score 全为 NA，无法自动生成 Potency category。")
    }
    
    qs <- quantile(score_vec, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    
    out <- dplyr::case_when(
      score_vec <= qs[1] ~ "Differentiated",
      score_vec <= qs[2] ~ "Unipotent",
      score_vec <= qs[3] ~ "Oligopotent",
      score_vec >  qs[3] ~ "Multipotent",
      TRUE ~ NA_character_
    )
  }
  
  out <- as.character(out)
  
  if (length(out) != ncell) {
    stop("提取/构造出的 Potency category 长度与细胞数不一致。")
  }
  
  out[is.na(out)] <- "Unknown"
  
  out <- factor(
    out,
    levels = c("Multipotent", "Oligopotent", "Unipotent", "Differentiated", "Unknown")
  )
  
  return(out)
}

format_p <- function(p) {
  ifelse(is.na(p), NA,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

#############################
## 2.1 Potency图配色和主题
#############################
potency_cols <- c(
  "Multipotent"    = "#F0CF6A",
  "Oligopotent"    = "#D9E97A",
  "Unipotent"      = "#58B89A",
  "Differentiated" = "#6153A6",
  "Unknown"        = "grey70"
)

potency_theme <- theme_classic(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = "#EBEBEB", color = "black", linewidth = 0.8),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_rect(fill = NA, color = "black", linewidth = 0.8),
    axis.line         = element_blank(),
    axis.ticks        = element_blank(),
    axis.text.x       = element_blank(),
    axis.text.y       = element_blank(),
    axis.title.x      = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y      = element_text(size = 13, face = "bold", color = "black"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    legend.title      = element_text(size = 13, face = "bold"),
    legend.text       = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA)
  )

make_potency_dimplot <- function(seurat_obj, reduction_use = "umap", title_txt = NULL, pt_size = 0.35) {
  if (is.null(title_txt)) {
    title_txt <- "CytoTRACE2 Differentiation Potency"
  }
  
  axis_prefix <- toupper(reduction_use)
  
  p <- DimPlot(
    object = seurat_obj,
    reduction = reduction_use,
    group.by = "Potency_Category",
    cols = potency_cols,
    pt.size = pt_size,
    shuffle = TRUE
  ) +
    labs(
      title = title_txt,
      x = paste0(axis_prefix, "1"),
      y = paste0(axis_prefix, "2"),
      color = "Potency category"
    ) +
    potency_theme
  
  return(p)
}

#############################
## 3. 创建目录并读取对象
#############################
safe_mkdir(ANNO_DIR)
safe_mkdir(TRACE_DIR)

check_file_exists(SEURAT_FILE)

msg("读取Seurat对象...")
seurat_obj <- qread(SEURAT_FILE)

check_seurat_object(
  seurat_obj = seurat_obj,
  celltype_col = CELLTYPE_COL,
  sample_col   = SAMPLE_COL,
  group_col    = GROUP_COL
)

msg("Seurat对象读取完成")
msg("细胞数: ", ncol(seurat_obj))
msg("基因数: ", nrow(seurat_obj))

#############################
## 4. 运行 CytoTRACE2
#############################
if (!file.exists(CYTO_QS_FILE)) {
  msg("开始运行CytoTRACE2 ...")
  setwd(TRACE_DIR)
  
  cytotrace2_result <- cytotrace2(
    input     = seurat_obj,
    species   = SPECIES,
    is_seurat = TRUE,
    slot_type = SLOT_TYPE,
    ncores    = NCORES,
    seed      = 1234
  )
  
  msg("CytoTRACE2运行完成，保存结果...")
  qsave(cytotrace2_result, CYTO_QS_FILE)
  saveRDS(cytotrace2_result, CYTO_RDS_FILE)
} else {
  msg("检测到已有CytoTRACE2结果，直接读取...")
  cytotrace2_result <- qread(CYTO_QS_FILE)
}

#############################
## 5. 生成plotData图对象
#############################
msg("生成plotData图对象...")

phe.anno <- data.frame(
  celltype = seurat_obj@meta.data[[CELLTYPE_COL]],
  row.names = colnames(seurat_obj)
)

plots <- plotData(
  cytotrace2_result = cytotrace2_result,
  expression_data   = NULL,
  is_seurat         = TRUE,
  annotation        = phe.anno,
  slot_type         = SLOT_TYPE,
  pc_dims           = PC_DIMS
)

qsave(plots, PLOTS_RDS_FILE)

#############################
## 6. 整理CytoTRACE2元数据
#############################
msg("提取CytoTRACE2 score并整理元数据...")

cyto_score <- extract_cytotrace2_score(cytotrace2_result, seurat_obj)
cyto_score <- as.numeric(cyto_score)
cyto_score[is.infinite(cyto_score)] <- NA
cyto_score[is.nan(cyto_score)] <- NA

if (all(is.na(cyto_score))) {
  stop("CytoTRACE2_Score 全部为 NA，请检查 cytotrace2_result 的实际结构。")
}

if (any(is.na(cyto_score))) {
  msg("检测到部分 NA，进行填补...")
  fill_value <- min(cyto_score[!is.na(cyto_score)])
  cyto_score[is.na(cyto_score)] <- fill_value
}

seurat_obj$CytoTRACE2_Score <- cyto_score

potency_cat <- extract_potency_category(
  cytotrace2_result = cytotrace2_result,
  seurat_obj = seurat_obj,
  cyto_score = seurat_obj$CytoTRACE2_Score
)

seurat_obj$Potency_Category <- potency_cat

msg("CytoTRACE2_Score summary:")
print(summary(seurat_obj$CytoTRACE2_Score))

msg("Potency category分布:")
print(table(seurat_obj$Potency_Category, useNA = "ifany"))

cell_df <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    CytoTRACE2_Score = .data$CytoTRACE2_Score,
    Potency_Category = .data$Potency_Category,
    celltype = .data[[CELLTYPE_COL]],
    sample_id = .data[[SAMPLE_COL]],
    group = .data[[GROUP_COL]]
  )

write.csv(cell_df, CELL_META_FILE, row.names = FALSE)

#############################
## 6.1 准备UMAP/tSNE
#############################
if (!"pca" %in% names(seurat_obj@reductions)) {
  stop("Seurat对象中未找到pca降维结果。请先确保对象已完成RunPCA。")
}

if (!"umap" %in% names(seurat_obj@reductions)) {
  msg("未检测到umap降维结果，开始计算RunUMAP ...")
  seurat_obj <- RunUMAP(
    object = seurat_obj,
    reduction = "pca",
    dims = 1:min(30, PC_DIMS),
    seed.use = 1234
  )
  msg("RunUMAP完成")
} else {
  msg("检测到已有umap降维结果，直接使用")
}

if (!"tsne" %in% names(seurat_obj@reductions)) {
  msg("未检测到tsne降维结果，开始计算RunTSNE ...")
  seurat_obj <- RunTSNE(
    object = seurat_obj,
    reduction = "pca",
    dims = TSNE_DIMS,
    check_duplicates = FALSE,
    seed.use = 1234
  )
  msg("RunTSNE完成")
} else {
  msg("检测到已有tsne降维结果，直接使用")
}

#############################
## 7. 基础图件
#############################
msg("输出基础图件...")

# 7.1 UMAP
p_umap <- make_potency_dimplot(
  seurat_obj = seurat_obj,
  reduction_use = "umap",
  title_txt = "CytoTRACE2 Differentiation Potency",
  pt_size = 0.35
)

ggsave(
  filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP.png"),
  plot = p_umap, width = 8, height = 6, dpi = 600, bg = "white"
)
save_pdf(file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP.pdf"), p_umap, 8, 6)

# 7.2 tSNE
p_tsne <- make_potency_dimplot(
  seurat_obj = seurat_obj,
  reduction_use = "tsne",
  title_txt = "CytoTRACE2 Differentiation Potency",
  pt_size = 0.35
)

ggsave(
  filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_tSNE.png"),
  plot = p_tsne, width = 8, height = 6, dpi = 600, bg = "white"
)
save_pdf(file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_tSNE.pdf"), p_tsne, 8, 6)

# 7.3 UMAP + tSNE 拼图
p_umap_tsne <- p_umap + p_tsne + plot_layout(ncol = 2)

ggsave(
  filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP_tSNE.png"),
  plot = p_umap_tsne, width = 14, height = 6, dpi = 600, bg = "white"
)
save_pdf(file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP_tSNE.pdf"), p_umap_tsne, 14, 6)

# 7.4 原始 boxplot
p_box_raw <- plots$CytoTRACE2_Boxplot_byPheno +
  ggtitle("CytoTRACE2 Potency Score by Cell Type") +
  xlab("Cell Type") +
  ylab("CytoTRACE2 Potency Score") +
  pub_theme +
  theme(legend.position = "none")

ggsave(
  filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_raw.png"),
  plot = p_box_raw, width = 10, height = 6, dpi = 600
)
save_pdf(file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_raw.pdf"), p_box_raw, 10, 6)

# 7.5 密度图
if (!is.null(plots$CytoTRACE2_Density)) {
  p_density <- plots$CytoTRACE2_Density +
    ggtitle("Distribution of CytoTRACE2 Potency Score") +
    xlab("CytoTRACE2 Score") +
    ylab("Density") +
    pub_theme +
    theme(legend.position = "right")
  
  ggsave(
    filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Density.png"),
    plot = p_density, width = 8, height = 6, dpi = 600
  )
  save_pdf(file.path(TRACE_DIR, "Fig_CytoTRACE2_Density.pdf"), p_density, 8, 6)
}

#############################
## 8. 单细胞层面统计
#############################
msg("进行单细胞层面统计...")

tmp_cell <- cell_df %>%
  filter(!is.na(celltype), !is.na(group), !is.na(CytoTRACE2_Score)) %>%
  group_by(celltype, group) %>%
  mutate(n_in_group = n()) %>%
  ungroup() %>%
  filter(n_in_group >= MIN_CELLS_PER_GROUP)

valid_ct <- tmp_cell %>%
  dplyr::group_by(.data$celltype) %>%
  dplyr::summarise(n_group = dplyr::n_distinct(.data$group), .groups = "drop") %>%
  dplyr::filter(.data$n_group >= 2) %>%
  dplyr::pull("celltype")

tmp_cell <- tmp_cell %>% dplyr::filter(.data$celltype %in% valid_ct)

if (nrow(tmp_cell) > 0) {
  celltype_list <- split(tmp_cell, tmp_cell$celltype)
  
  cell_stats_list <- lapply(names(celltype_list), function(ct) {
    dat_ct <- celltype_list[[ct]]
    if (length(unique(dat_ct$group)) < 2) return(NULL)
    
    wt <- rstatix::wilcox_test(dat_ct, CytoTRACE2_Score ~ group)
    wt$celltype <- ct
    wt
  })
  
  cell_stats <- dplyr::bind_rows(cell_stats_list)
  
  if (nrow(cell_stats) > 0) {
    cell_stats <- cell_stats %>%
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance("p.adj") %>%
      dplyr::mutate(
        p.raw = p,
        p.adj.raw = p.adj,
        p = format_p(p),
        p.adj = format_p(p.adj)
      )
    
    write.csv(cell_stats, CELL_STATS_FILE, row.names = FALSE)
    
    stat_y <- tmp_cell %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(
        y.position = max(CytoTRACE2_Score, na.rm = TRUE) * 1.08,
        .groups = "drop"
      )
    
    cell_stats_plot <- cell_stats %>%
      dplyr::left_join(stat_y, by = "celltype") %>%
      dplyr::mutate(
        xmin = celltype,
        xmax = celltype
      )
    
    p_box_cell_level <- ggplot(
      tmp_cell,
      aes(x = celltype, y = CytoTRACE2_Score, fill = group)
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.75,
        alpha = 0.85,
        position = position_dodge(width = 0.8)
      ) +
      geom_point(
        aes(color = group),
        size = 0.6,
        alpha = 0.35,
        position = position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = 0.8
        ),
        show.legend = FALSE
      ) +
      ggpubr::stat_pvalue_manual(
        cell_stats_plot,
        label = "p.adj.signif",
        xmin = "xmin",
        xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        inherit.aes = FALSE
      ) +
      scale_fill_viridis_d(option = "D", end = 0.85) +
      scale_color_viridis_d(option = "D", end = 0.85) +
      labs(
        title = "CytoTRACE2 Score by Cell Type (Cell-level comparison)",
        subtitle = "Cell-level statistics; interpret cautiously",
        x = "Cell Type",
        y = "CytoTRACE2 Potency Score"
      ) +
      pub_theme +
      theme(legend.position = "right")
    
    ggsave(
      filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_cellLevelStats.png"),
      plot = p_box_cell_level,
      width = 11,
      height = 6.5,
      dpi = 600
    )
    
    save_pdf(
      file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_cellLevelStats.pdf"),
      p_box_cell_level,
      11,
      6.5
    )
  } else {
    msg("单细胞层面统计跳过：没有可用于检验的 cell type。")
  }
} else {
  msg("单细胞层面统计跳过：没有满足条件的细胞。")
}
#############################
## 9. 样本层面统计（推荐）
#############################
msg("进行样本层面统计...")

sample_input <- cell_df %>%
  dplyr::transmute(
    sample_id = .data$sample_id,
    group = .data$group,
    celltype = .data$celltype,
    CytoTRACE2_Score = .data$CytoTRACE2_Score
  ) %>%
  dplyr::filter(
    !is.na(.data$sample_id),
    !is.na(.data$group),
    !is.na(.data$celltype),
    !is.na(.data$CytoTRACE2_Score)
  )

sample_summary <- sample_input %>%
  dplyr::group_by(.data$sample_id, .data$group, .data$celltype) %>%
  dplyr::summarise(
    mean_score = mean(.data$CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(.data$CytoTRACE2_Score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$n_cells >= MIN_CELLS_PER_GROUP)

write.csv(sample_summary, SAMPLE_SUM_FILE, row.names = FALSE)

sample_group_count <- sample_summary %>%
  dplyr::group_by(.data$celltype, .data$group) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(.data$sample_id),
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$n_samples >= MIN_SAMPLES_PER_GROUP)

valid_ct_sample <- sample_group_count %>%
  dplyr::group_by(.data$celltype) %>%
  dplyr::summarise(
    n_group = dplyr::n_distinct(.data$group),
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$n_group >= 2) %>%
  dplyr::pull("celltype")

sample_summary2 <- sample_summary %>%
  dplyr::filter(.data$celltype %in% valid_ct_sample)

if (nrow(sample_summary2) > 0) {
  sampletype_list <- split(sample_summary2, sample_summary2$celltype)
  
  sample_stats_list <- lapply(names(sampletype_list), function(ct) {
    dat_ct <- sampletype_list[[ct]]
    if (length(unique(dat_ct$group)) < 2) return(NULL)
    wt <- rstatix::wilcox_test(dat_ct, mean_score ~ group)
    wt$celltype <- ct
    wt
  })
  
  sample_stats <- dplyr::bind_rows(sample_stats_list)
  
  if (nrow(sample_stats) > 0) {
    sample_stats <- sample_stats %>%
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance("p.adj") %>%
      dplyr::mutate(
        p.raw = p,
        p.adj.raw = p.adj,
        p = format_p(p),
        p.adj = format_p(p.adj)
      )
    
    write.csv(sample_stats, SAMPLE_STATS_FILE, row.names = FALSE)
    
    stat_y2 <- sample_summary2 %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(
        y.position = max(mean_score, na.rm = TRUE) * 1.12,
        .groups = "drop"
      )
    
    sample_stats_plot <- sample_stats %>%
      dplyr::left_join(stat_y2, by = "celltype") %>%
      dplyr::mutate(
        xmin = celltype,
        xmax = celltype
      )
    
    p_box_sample_level <- ggplot(
      sample_summary2,
      aes(x = celltype, y = mean_score, fill = group)
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.72,
        alpha = 0.85,
        position = position_dodge(width = 0.8)
      ) +
      geom_point(
        aes(color = group),
        size = 2.2,
        alpha = 0.85,
        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8)
      ) +
      ggpubr::stat_pvalue_manual(
        sample_stats_plot,
        label = "p.adj.signif",
        xmin = "xmin",
        xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        inherit.aes = FALSE
      ) +
      scale_fill_viridis_d(option = "D", end = 0.85) +
      scale_color_viridis_d(option = "D", end = 0.85) +
      labs(
        title = "CytoTRACE2 Score by Cell Type (Sample-level comparison)",
        subtitle = "Recommended for manuscript-level inference",
        x = "Cell Type",
        y = "Mean CytoTRACE2 Score per Sample"
      ) +
      pub_theme +
      theme(legend.position = "right")
    
    ggsave(
      filename = file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_sampleLevelStats.png"),
      plot = p_box_sample_level, width = 11, height = 6.5, dpi = 600
    )
    save_pdf(
      file.path(TRACE_DIR, "Fig_CytoTRACE2_Boxplot_byCellType_sampleLevelStats.pdf"),
      p_box_sample_level, 11, 6.5
    )
  } else {
    msg("样本层面统计跳过：没有可用于检验的 cell type。")
  }
} else {
  msg("样本层面统计跳过：没有满足“每组至少2个样本且每个样本至少10个细胞”的cell type。")
}

#############################
## 10. sessionInfo
#############################
writeLines(capture.output(sessionInfo()), SESSIONINFO_FILE)

msg("分析完成")








rm(list = ls())

############################################################
## Redraw CytoTRACE2 Potency UMAP / tSNE
## Only change grey background to white
## Do NOT rerun CytoTRACE2
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

############################################################
## 1. Parameters
############################################################

ANNO_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/2.Sc_annotation/new"
TRACE_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL/12.CytoTRACE2/0.new"

SEURAT_FILE <- file.path(ANNO_DIR, "2.resolution_0.5_annotation.qs")
PLOTS_RDS_FILE <- file.path(TRACE_DIR, "cytotrace2_plots.rds")
CELL_META_FILE <- file.path(TRACE_DIR, "CytoTRACE2_cell_metadata.csv")
CYTO_QS_FILE <- file.path(TRACE_DIR, "cytotrace2_result.qs")

## TRUE  = 覆盖原来的 Fig_CytoTRACE2_Potency_UMAP / tSNE / UMAP_tSNE
## FALSE = 输出带 _white_background 后缀的新文件
OVERWRITE_ORIGINAL <- FALSE

############################################################
## 2. Helper functions
############################################################

check_file_exists <- function(x) {
  if (!file.exists(x)) stop("File not found: ", x)
}

read_saved_object <- function(file) {
  obj <- tryCatch(qread(file), error = function(e) NULL)
  if (is.null(obj)) {
    obj <- tryCatch(readRDS(file), error = function(e) NULL)
  }
  if (is.null(obj)) {
    stop("Cannot read object file by qread() or readRDS(): ", file)
  }
  obj
}

save_plot_both <- function(plot_obj,
                           file_prefix,
                           width = 8,
                           height = 6,
                           dpi = 600) {
  
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

## 只改背景，不改颜色、图例、点大小等
make_panel_white_only <- function(p) {
  p +
    theme(
      panel.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.8
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      legend.background = element_rect(
        fill = "white",
        color = NA
      ),
      legend.key = element_rect(
        fill = "white",
        color = NA
      )
    )
}

############################################################
## 3. Original colors and theme from your previous script
############################################################

## 这一组颜色来自你原始 CytoTRACE2 代码
## 保证颜色与之前保持一致
potency_cols <- c(
  "Multipotent"    = "#F0CF6A",
  "Oligopotent"    = "#D9E97A",
  "Unipotent"      = "#58B89A",
  "Differentiated" = "#6153A6",
  "Unknown"        = "grey70"
)

potency_theme_white <- theme_classic(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_rect(fill = NA, color = "black", linewidth = 0.8),
    axis.line         = element_blank(),
    axis.ticks        = element_blank(),
    axis.text.x       = element_blank(),
    axis.text.y       = element_blank(),
    axis.title.x      = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y      = element_text(size = 13, face = "bold", color = "black"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    legend.title      = element_text(size = 13, face = "bold"),
    legend.text       = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA)
  )

make_potency_dimplot_white <- function(seurat_obj,
                                       reduction_use = "umap",
                                       title_txt = "CytoTRACE2 Differentiation Potency",
                                       pt_size = 0.35) {
  
  if (!reduction_use %in% names(seurat_obj@reductions)) {
    stop(
      "Reduction not found: ",
      reduction_use,
      "\nAvailable reductions: ",
      paste(names(seurat_obj@reductions), collapse = ", ")
    )
  }
  
  axis_prefix <- toupper(reduction_use)
  
  present_levels <- levels(factor(seurat_obj$Potency_Category))
  cols_use <- potency_cols[names(potency_cols) %in% present_levels]
  
  DimPlot(
    object = seurat_obj,
    reduction = reduction_use,
    group.by = "Potency_Category",
    cols = cols_use,
    pt.size = pt_size,
    shuffle = TRUE
  ) +
    labs(
      title = title_txt,
      x = paste0(axis_prefix, "1"),
      y = paste0(axis_prefix, "2"),
      color = "Potency category"
    ) +
    potency_theme_white
}

############################################################
## 4. Read existing plots and metadata
############################################################

check_file_exists(PLOTS_RDS_FILE)
check_file_exists(SEURAT_FILE)
check_file_exists(CELL_META_FILE)

plots <- read_saved_object(PLOTS_RDS_FILE)

cat("plots 中包含的对象名称：\n")
print(names(plots))

seurat_obj <- qread(SEURAT_FILE)

cell_meta <- read.csv(
  CELL_META_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"cell_id" %in% colnames(cell_meta)) {
  stop("CytoTRACE2_cell_metadata.csv 中未找到 cell_id 列。")
}

if (!"Potency_Category" %in% colnames(cell_meta)) {
  stop("CytoTRACE2_cell_metadata.csv 中未找到 Potency_Category 列。")
}

common_cells <- intersect(colnames(seurat_obj), cell_meta$cell_id)

if (length(common_cells) == 0) {
  stop("Seurat 对象细胞名与 CytoTRACE2_cell_metadata.csv 中 cell_id 无交集。")
}

cell_meta_use <- cell_meta %>%
  filter(cell_id %in% common_cells) %>%
  distinct(cell_id, .keep_all = TRUE)

rownames(cell_meta_use) <- cell_meta_use$cell_id

seurat_obj <- subset(seurat_obj, cells = common_cells)
cell_meta_use <- cell_meta_use[colnames(seurat_obj), , drop = FALSE]

seurat_obj$Potency_Category <- cell_meta_use$Potency_Category

seurat_obj$Potency_Category <- factor(
  seurat_obj$Potency_Category,
  levels = c(
    "Multipotent",
    "Oligopotent",
    "Unipotent",
    "Differentiated",
    "Unknown"
  )
)

############################################################
## 5. Recover tSNE / UMAP reductions if needed
############################################################

## 如果 annotation 对象里没有 tSNE，尝试从 cytotrace2_result.qs 里找
if (!"tsne" %in% names(seurat_obj@reductions) ||
    !"umap" %in% names(seurat_obj@reductions)) {
  
  if (file.exists(CYTO_QS_FILE)) {
    cytotrace2_result <- qread(CYTO_QS_FILE)
    
    if (inherits(cytotrace2_result, "Seurat")) {
      
      if (!"umap" %in% names(seurat_obj@reductions) &&
          "umap" %in% names(cytotrace2_result@reductions)) {
        seurat_obj[["umap"]] <- cytotrace2_result[["umap"]]
      }
      
      if (!"tsne" %in% names(seurat_obj@reductions) &&
          "tsne" %in% names(cytotrace2_result@reductions)) {
        seurat_obj[["tsne"]] <- cytotrace2_result[["tsne"]]
      }
    }
  }
}

############################################################
## 6. Draw UMAP and tSNE
############################################################

## UMAP：如果 plots 中有原始 UMAP 图，就只改白底；
## 否则使用 Seurat reduction 重新绘制，但颜色仍保持原代码颜色。
if ("CytoTRACE2_Potency_UMAP" %in% names(plots)) {
  p_umap <- make_panel_white_only(
    plots$CytoTRACE2_Potency_UMAP
  )
} else {
  p_umap <- make_potency_dimplot_white(
    seurat_obj = seurat_obj,
    reduction_use = "umap",
    title_txt = "CytoTRACE2 Differentiation Potency",
    pt_size = 0.35
  )
}

## tSNE：plots 中通常没有，所以这里用 Seurat reduction 重画
if ("CytoTRACE2_Potency_tSNE" %in% names(plots)) {
  p_tsne <- make_panel_white_only(
    plots$CytoTRACE2_Potency_tSNE
  )
} else {
  p_tsne <- make_potency_dimplot_white(
    seurat_obj = seurat_obj,
    reduction_use = "tsne",
    title_txt = "CytoTRACE2 Differentiation Potency",
    pt_size = 0.35
  )
}

p_umap_tsne <- p_umap + p_tsne + patchwork::plot_layout(ncol = 2)

############################################################
## 7. Save figures
############################################################

if (OVERWRITE_ORIGINAL) {
  prefix_umap <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP")
  prefix_tsne <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_tSNE")
  prefix_both <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP_tSNE")
} else {
  prefix_umap <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP_white_background")
  prefix_tsne <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_tSNE_white_background")
  prefix_both <- file.path(TRACE_DIR, "Fig_CytoTRACE2_Potency_UMAP_tSNE_white_background")
}

save_plot_both(
  p_umap,
  prefix_umap,
  width = 8,
  height = 6,
  dpi = 600
)

save_plot_both(
  p_tsne,
  prefix_tsne,
  width = 8,
  height = 6,
  dpi = 600
)

save_plot_both(
  p_umap_tsne,
  prefix_both,
  width = 14,
  height = 6,
  dpi = 600
)

cat("CytoTRACE2 white-background plots have been exported:\n")
cat(" - ", paste0(prefix_umap, ".pdf / .png"), "\n", sep = "")
cat(" - ", paste0(prefix_tsne, ".pdf / .png"), "\n", sep = "")
cat(" - ", paste0(prefix_both, ".pdf / .png"), "\n", sep = "")

# END ORIGINAL SCRIPT P14

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P15: 09_SC_CYTOTRACE2/plot_cytotrace2_pub.R
# ==============================================================================
############################################################
## 04_plot_cytotrace2_pub_original_order.R
## Re-plot CytoTRACE2 stemness analysis results
## Matched with final annotation / proportion / CellChat styles
############################################################

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(grid)
})

source("00_pub_theme_colors.R")

############################################################
## 0. Compatibility helpers
############################################################
if (!exists("standardize_group")) {
  standardize_group <- function(x) {
    x <- as.character(x)
    dplyr::case_when(
      grepl("^AML", x, ignore.case = TRUE) ~ "AML",
      grepl("healthy|normal|control|^BM|^NBM", x, ignore.case = TRUE) ~ "Healthy",
      x %in% c("AML", "Healthy") ~ x,
      x %in% c("Normal", "NORMAL") ~ "Healthy",
      TRUE ~ x
    )
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


############################################################
## 0.1 Helpers for preserving original CytoTRACE2 plot order
############################################################
read_cytotrace2_plots <- function(file) {
  if (!file.exists(file)) return(NULL)
  obj <- tryCatch(qs::qread(file), error = function(e) NULL)
  if (is.null(obj)) obj <- tryCatch(readRDS(file), error = function(e) NULL)
  obj
}

map_celltype_to_short <- function(x) {
  x <- as.character(x)
  mapped <- unname(celltype_map[x])
  mapped[is.na(mapped) | mapped == ""] <- x[is.na(mapped) | mapped == ""]
  mapped
}

extract_cytotrace2_original_order <- function(plots_obj, available_celltypes) {
  available_celltypes <- as.character(available_celltypes)
  if (is.null(plots_obj) || is.null(plots_obj$CytoTRACE2_Boxplot_byPheno)) {
    return(available_celltypes)
  }
  p <- plots_obj$CytoTRACE2_Boxplot_byPheno

  # 1) Try the data embedded in the original CytoTRACE2 ggplot object.
  if (!is.null(p$data) && is.data.frame(p$data) && nrow(p$data) > 0) {
    for (cc in colnames(p$data)) {
      v <- p$data[[cc]]
      if (is.factor(v)) {
        cand <- map_celltype_to_short(levels(v))
        cand <- cand[cand %in% available_celltypes]
        if (length(cand) >= 2) {
          return(c(cand, setdiff(available_celltypes, cand)))
        }
      }
    }
    for (cc in colnames(p$data)) {
      v <- p$data[[cc]]
      if (is.character(v) || is.factor(v)) {
        cand <- map_celltype_to_short(unique(as.character(v)))
        cand <- cand[cand %in% available_celltypes]
        if (length(cand) >= 2) {
          return(c(cand, setdiff(available_celltypes, cand)))
        }
      }
    }
  }

  # 2) Try the built ggplot scale labels.
  cand <- tryCatch({
    gb <- ggplot2::ggplot_build(p)
    pp <- gb$layout$panel_params[[1]]
    labs <- NULL
    if (!is.null(pp$x$get_labels)) labs <- pp$x$get_labels()
    if (is.null(labs) && !is.null(pp$x$labels)) labs <- pp$x$labels
    map_celltype_to_short(labs)
  }, error = function(e) character(0))
  cand <- cand[cand %in% available_celltypes]
  if (length(cand) >= 2) return(c(cand, setdiff(available_celltypes, cand)))

  available_celltypes
}

############################################################
## 1. Paths and core settings
############################################################
OUT_DIR <- file.path(TRACE_DIR, "pub_figures_cytotrace2")
safe_mkdir(OUT_DIR)

SEURAT_FILE <- SEURAT_ANNO_FILE
CELL_META_FILE <- file.path(TRACE_DIR, "CytoTRACE2_cell_metadata.csv")
CYTO_QS_FILE <- file.path(TRACE_DIR, "cytotrace2_result.qs")
CYTO_PLOTS_FILE <- file.path(TRACE_DIR, "cytotrace2_plots.rds")

CELLTYPE_COL <- "celltype"
SAMPLE_COL   <- "orig.ident"
GROUP_CANDIDATES <- c("type", "sample_group", "group")

POINT_SIZE_DIM <- 0.42
POINT_ALPHA_DIM <- 0.95
SHORT_AXIS_TEXT_SIZE <- 4.0
SHORT_AXIS_LINE_WIDTH <- 0.75

############################################################
## 2. Read annotated Seurat object and CytoTRACE2 metadata
############################################################
if (!file.exists(SEURAT_FILE)) {
  stop("Cannot find annotated Seurat object: ", SEURAT_FILE)
}
if (!file.exists(CELL_META_FILE)) {
  stop("Cannot find CytoTRACE2 cell metadata: ", CELL_META_FILE,
       "\n请先确认 16.单细胞cytoTRACE干性分析.R 已输出 CytoTRACE2_cell_metadata.csv。")
}

seurat_obj <- qread(SEURAT_FILE)

if (!CELLTYPE_COL %in% colnames(seurat_obj@meta.data)) {
  seurat_obj@meta.data[[CELLTYPE_COL]] <- as.character(Idents(seurat_obj))
}
if (!"celltype_short" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj <- add_celltype_short(seurat_obj, celltype_col = CELLTYPE_COL, short_col = "celltype_short")
}

group_col <- GROUP_CANDIDATES[GROUP_CANDIDATES %in% colnames(seurat_obj@meta.data)][1]
if (is.na(group_col)) {
  if (SAMPLE_COL %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$type <- ifelse(grepl("^AML", seurat_obj@meta.data[[SAMPLE_COL]], ignore.case = TRUE), "AML", "Healthy")
    group_col <- "type"
  } else {
    seurat_obj@meta.data$type <- "All"
    group_col <- "type"
  }
}

cytotrace2_plots <- read_cytotrace2_plots(CYTO_PLOTS_FILE)

cell_meta <- read.csv(CELL_META_FILE, stringsAsFactors = FALSE, check.names = FALSE)

## infer cell_id column robustly
if (!"cell_id" %in% colnames(cell_meta)) {
  candidate_cell_cols <- c("Cell", "cell", "barcode", "Barcode", "barcodes", "cells", "X", "...1")
  hit <- candidate_cell_cols[candidate_cell_cols %in% colnames(cell_meta)][1]
  if (is.na(hit)) {
    overlap_n <- sapply(colnames(cell_meta), function(cc) {
      sum(as.character(cell_meta[[cc]]) %in% colnames(seurat_obj), na.rm = TRUE)
    })
    if (max(overlap_n, na.rm = TRUE) > 0) {
      hit <- names(which.max(overlap_n))
    }
  }
  if (!is.na(hit)) {
    cell_meta$cell_id <- as.character(cell_meta[[hit]])
  } else {
    stop("Cannot identify cell barcode column in CytoTRACE2_cell_metadata.csv. Available columns: ",
         paste(colnames(cell_meta), collapse = ", "))
  }
}
cell_meta$cell_id <- as.character(cell_meta$cell_id)

score_candidates <- c(
  "CytoTRACE2_Score", "CytoTRACE2_score", "cytotrace2_score",
  "CytoTRACE_score", "CytoTRACE", "preKNN_CytoTRACE2",
  "Score", "score"
)
score_col <- score_candidates[score_candidates %in% colnames(cell_meta)][1]
if (is.na(score_col)) {
  numeric_cols <- colnames(cell_meta)[sapply(cell_meta, is.numeric)]
  stop("Cannot find CytoTRACE2 score column. Numeric columns detected: ",
       paste(numeric_cols, collapse = ", "),
       "\nAvailable columns: ", paste(colnames(cell_meta), collapse = ", "))
}

potency_candidates <- c(
  "Potency_Category", "potency_category", "CytoTRACE2_Potency",
  "CytoTRACE2_Category", "CytoTRACE_Category", "Potency", "potency",
  "category", "label", "labels"
)
potency_col <- potency_candidates[potency_candidates %in% colnames(cell_meta)][1]

common_cells <- intersect(colnames(seurat_obj), cell_meta$cell_id)
if (length(common_cells) == 0) {
  stop("No shared cell barcodes between Seurat object and CytoTRACE2 metadata.")
}

seurat_obj <- subset(seurat_obj, cells = common_cells)
cell_meta_use <- cell_meta %>%
  filter(cell_id %in% common_cells) %>%
  distinct(cell_id, .keep_all = TRUE)
rownames(cell_meta_use) <- cell_meta_use$cell_id
cell_meta_use <- cell_meta_use[colnames(seurat_obj), , drop = FALSE]

seurat_obj@meta.data$CytoTRACE2_Score <- as.numeric(cell_meta_use[[score_col]])

if (!is.na(potency_col)) {
  potency_raw <- as.character(cell_meta_use[[potency_col]])
} else {
  score_vec <- seurat_obj@meta.data$CytoTRACE2_Score
  qs <- quantile(score_vec, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  potency_raw <- dplyr::case_when(
    score_vec <= qs[1] ~ "Differentiated",
    score_vec <= qs[2] ~ "Unipotent",
    score_vec <= qs[3] ~ "Oligopotent",
    score_vec >  qs[3] ~ "Multipotent",
    TRUE ~ "Unknown"
  )
}
potency_raw[is.na(potency_raw) | potency_raw == ""] <- "Unknown"
seurat_obj@meta.data$Potency_Category <- factor(
  potency_raw,
  levels = c("Multipotent", "Oligopotent", "Unipotent", "Differentiated", "Unknown")
)

seurat_obj@meta.data$disease <- standardize_group(seurat_obj@meta.data[[group_col]])
seurat_obj@meta.data$sample_id <- if (SAMPLE_COL %in% colnames(seurat_obj@meta.data)) {
  as.character(seurat_obj@meta.data[[SAMPLE_COL]])
} else {
  "Sample"
}
seurat_obj <- add_celltype_short(seurat_obj, celltype_col = CELLTYPE_COL, short_col = "celltype_short")

celltype_levels <- c(
  intersect(celltype_levels_short, unique(as.character(seurat_obj@meta.data$celltype_short))),
  setdiff(unique(as.character(seurat_obj@meta.data$celltype_short)), celltype_levels_short)
)
seurat_obj@meta.data$celltype_short <- factor(as.character(seurat_obj@meta.data$celltype_short), levels = celltype_levels)

disease_levels <- unique(as.character(seurat_obj@meta.data$disease))
if (all(c("Healthy", "AML") %in% disease_levels)) {
  disease_levels <- c("Healthy", "AML", setdiff(disease_levels, c("Healthy", "AML")))
}
seurat_obj@meta.data$disease <- factor(as.character(seurat_obj@meta.data$disease), levels = disease_levels)

cell_cols <- get_celltype_colors(celltype_levels, use_short = TRUE)
group_cols <- get_named_colors(disease_levels, palette_named = group_colors)
pot_cols <- potency_category_colors[names(potency_category_colors) %in% levels(droplevels(seurat_obj@meta.data$Potency_Category))]

meta_out <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(any_of(c(
    "cell_id", "sample_id", "disease", CELLTYPE_COL, "celltype_short",
    "CytoTRACE2_Score", "Potency_Category"
  )))
write.csv(meta_out, file.path(OUT_DIR, "Tab01_CytoTRACE2_cell_metadata_pub.csv"), row.names = FALSE)

############################################################
## 3. Short-axis UMAP/tSNE helpers
############################################################
get_reduction_df <- function(obj, reduction_use = "umap") {
  if (!reduction_use %in% names(obj@reductions)) return(NULL)
  emb <- Embeddings(obj, reduction = reduction_use)
  if (ncol(emb) < 2) return(NULL)
  df <- as.data.frame(emb[, 1:2, drop = FALSE])
  colnames(df) <- c("dim1", "dim2")
  df$cell_id <- rownames(df)
  meta_use <- obj@meta.data %>% tibble::rownames_to_column("cell_id")
  df <- df %>% left_join(meta_use, by = "cell_id")
  df
}

add_short_axes <- function(p,
                           plot_df,
                           x_lab = "UMAP1",
                           y_lab = "UMAP2",
                           axis_frac = 0.115,
                           offset_frac = -0.055,
                           text_size = SHORT_AXIS_TEXT_SIZE,
                           line_width = SHORT_AXIS_LINE_WIDTH) {
  xr <- range(plot_df$dim1, na.rm = TRUE)
  yr <- range(plot_df$dim2, na.rm = TRUE)
  x_span <- diff(xr)
  y_span <- diff(yr)
  axis_len <- min(x_span, y_span) * axis_frac
  x0 <- xr[1] + x_span * offset_frac
  y0 <- yr[1] + y_span * offset_frac
  arrow_obj <- grid::arrow(length = grid::unit(0.16, "inches"), type = "closed")
  p +
    annotate("segment", x = x0, xend = x0 + axis_len, y = y0, yend = y0,
             linewidth = line_width, color = "black", arrow = arrow_obj) +
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + axis_len,
             linewidth = line_width, color = "black", arrow = arrow_obj) +
    annotate("text", x = x0 + axis_len * 0.52, y = y0 - y_span * 0.035,
             label = x_lab, size = text_size, fontface = "bold", hjust = 0.5, color = "black") +
    annotate("text", x = x0 - x_span * 0.035, y = y0 + axis_len * 0.55,
             label = y_lab, size = text_size, fontface = "bold", angle = 90, hjust = 0.5, color = "black") +
    coord_equal(clip = "off")
}

theme_dim_cyto <- function(base_size = 12) {
  theme_dim_shortaxis(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(size = 11, face = "bold", color = "black"),
      legend.text = element_text(size = 9.5, color = "black"),
      legend.key.size = unit(0.38, "cm"),
      plot.margin = margin(8, 8, 14, 14)
    )
}

plot_score_dim <- function(obj, reduction_use = "umap") {
  df <- get_reduction_df(obj, reduction_use)
  if (is.null(df)) return(NULL)
  axis_prefix <- toupper(reduction_use)
  p <- ggplot(df, aes(x = dim1, y = dim2, color = CytoTRACE2_Score)) +
    geom_point(size = POINT_SIZE_DIM, alpha = POINT_ALPHA_DIM, stroke = 0) +
    scale_color_gradientn(colors = cytotrace_score_cols, name = "CytoTRACE2\nscore") +
    labs(x = NULL, y = NULL) +
    theme_dim_cyto(12)
  add_short_axes(p, df, x_lab = paste0(axis_prefix, "1"), y_lab = paste0(axis_prefix, "2"))
}

plot_potency_dim <- function(obj, reduction_use = "umap") {
  df <- get_reduction_df(obj, reduction_use)
  if (is.null(df)) return(NULL)
  axis_prefix <- toupper(reduction_use)
  df$Potency_Category <- factor(as.character(df$Potency_Category),
                                levels = c("Multipotent", "Oligopotent", "Unipotent", "Differentiated", "Unknown"))
  cols_use <- potency_category_colors[names(potency_category_colors) %in% levels(droplevels(df$Potency_Category))]
  p <- ggplot(df, aes(x = dim1, y = dim2, color = Potency_Category)) +
    geom_point(size = POINT_SIZE_DIM, alpha = POINT_ALPHA_DIM, stroke = 0) +
    scale_color_manual(values = cols_use, drop = FALSE, name = "Potency") +
    guides(color = guide_legend(override.aes = list(size = 3.0, alpha = 1), ncol = 1)) +
    labs(x = NULL, y = NULL) +
    theme_dim_cyto(12)
  add_short_axes(p, df, x_lab = paste0(axis_prefix, "1"), y_lab = paste0(axis_prefix, "2"))
}

############################################################
## 4. Dimensional plots
############################################################
p_score_umap <- plot_score_dim(seurat_obj, "umap")
p_score_tsne <- plot_score_dim(seurat_obj, "tsne")
p_pot_umap <- plot_potency_dim(seurat_obj, "umap")
p_pot_tsne <- plot_potency_dim(seurat_obj, "tsne")

if (!is.null(p_score_umap)) save_plot_pub(p_score_umap, file.path(OUT_DIR, "Fig01_CytoTRACE2_score_UMAP_shortaxis"), width = 6.3, height = 5.8)
if (!is.null(p_score_tsne)) save_plot_pub(p_score_tsne, file.path(OUT_DIR, "Fig02_CytoTRACE2_score_tSNE_shortaxis"), width = 6.3, height = 5.8)
if (!is.null(p_pot_umap)) save_plot_pub(p_pot_umap, file.path(OUT_DIR, "Fig03_CytoTRACE2_potency_UMAP_shortaxis"), width = 6.6, height = 5.8)
if (!is.null(p_pot_tsne)) save_plot_pub(p_pot_tsne, file.path(OUT_DIR, "Fig04_CytoTRACE2_potency_tSNE_shortaxis"), width = 6.6, height = 5.8)

if (!is.null(p_score_umap) && !is.null(p_score_tsne)) {
  save_plot_pub(p_score_umap + p_score_tsne, file.path(OUT_DIR, "Fig05_CytoTRACE2_score_UMAP_tSNE"), width = 12.2, height = 5.8)
}
if (!is.null(p_pot_umap) && !is.null(p_pot_tsne)) {
  save_plot_pub(p_pot_umap + p_pot_tsne, file.path(OUT_DIR, "Fig06_CytoTRACE2_potency_UMAP_tSNE"), width = 12.6, height = 5.8)
}

############################################################
## 5. Score distributions by cell type and disease status
############################################################
plot_df <- seurat_obj@meta.data %>%
  mutate(
    celltype_short = factor(as.character(celltype_short), levels = celltype_levels),
    disease = factor(as.character(disease), levels = disease_levels)
  ) %>%
  filter(!is.na(CytoTRACE2_Score), !is.na(celltype_short), !is.na(disease))

celltype_summary <- plot_df %>%
  group_by(celltype_short) %>%
  summarise(
    mean_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(CytoTRACE2_Score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))
write.csv(celltype_summary, file.path(OUT_DIR, "Tab02_CytoTRACE2_score_by_celltype_summary.csv"), row.names = FALSE)

sample_summary <- plot_df %>%
  group_by(sample_id, disease, celltype_short) %>%
  summarise(
    mean_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(CytoTRACE2_Score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )
write.csv(sample_summary, file.path(OUT_DIR, "Tab03_CytoTRACE2_score_by_sample_celltype.csv"), row.names = FALSE)

## Fig07: stemness ranking by mean score
p_rank <- ggplot(celltype_summary,
                 aes(x = reorder(celltype_short, mean_score), y = mean_score, fill = celltype_short)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.18) +
  geom_text(aes(label = sprintf("%.3f", mean_score)), hjust = -0.08, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cell_cols, guide = "none") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Mean CytoTRACE2 score") +
  theme_sc_pub(base_size = 12) +
  theme(panel.border = element_blank())
save_plot_pub(p_rank, file.path(OUT_DIR, "Fig07_CytoTRACE2_celltype_stemness_ranking"), width = 6.6, height = max(5.8, 0.25 * length(celltype_levels) + 1.6))

## Fig08: boxplot by cell type
p_box_celltype <- ggplot(plot_df, aes(x = celltype_short, y = CytoTRACE2_Score, fill = celltype_short)) +
  geom_boxplot(width = 0.64, outlier.size = 0.15, linewidth = 0.32, alpha = 0.92) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  labs(x = NULL, y = "CytoTRACE2 score") +
  theme_sc_pub(base_size = 12) +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5)
  )
save_plot_pub(p_box_celltype, file.path(OUT_DIR, "Fig08_CytoTRACE2_score_boxplot_by_celltype"), width = 10.8, height = 5.5)

## Fig08b: original CytoTRACE2 package output order, with potency scale on right axis
## Important: this figure does NOT re-rank cell types by mean or median score.
## It follows the cell-type order from plots$CytoTRACE2_Boxplot_byPheno, i.e. the original CytoTRACE2 plotData output.
celltype_order_cytotrace2 <- extract_cytotrace2_original_order(
  plots_obj = cytotrace2_plots,
  available_celltypes = as.character(unique(plot_df$celltype_short))
)

write.csv(
  data.frame(order_index = seq_along(celltype_order_cytotrace2), celltype_short = celltype_order_cytotrace2),
  file.path(OUT_DIR, "Tab02b_CytoTRACE2_original_boxplot_celltype_order.csv"),
  row.names = FALSE
)

plot_df_cyto_order <- plot_df %>%
  mutate(celltype_short = factor(as.character(celltype_short), levels = celltype_order_cytotrace2))

cell_cols_cyto_order <- cell_cols[celltype_order_cytotrace2]
cell_cols_cyto_order <- cell_cols_cyto_order[!is.na(cell_cols_cyto_order)]

p_box_cyto_order <- ggplot(plot_df_cyto_order, aes(x = celltype_short, y = CytoTRACE2_Score, fill = celltype_short)) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.34,
    alpha = 0.55,
    color = "grey20"
  ) +
  geom_jitter(
    aes(color = celltype_short),
    width = 0.18,
    size = 0.28,
    alpha = 0.42,
    stroke = 0
  ) +
  scale_fill_manual(values = cell_cols_cyto_order, guide = "none") +
  scale_color_manual(values = cell_cols_cyto_order, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.03)),
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = c(0.08, 0.25, 0.43, 0.60, 0.75, 0.92),
      labels = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"),
      name = NULL
    )
  ) +
  labs(
    x = "Cell Type",
    y = "CytoTRACE2 Potency Score",
    title = "CytoTRACE2 Potency Score by Cell Type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55),
    axis.title = element_text(face = "bold", color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 9.2),
    axis.text.y = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black", size = 9.5),
    axis.ticks.y.right = element_line(color = "black"),
    axis.line.y.right = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 14),
    plot.margin = margin(8, 18, 8, 8)
  )

save_plot_pub(
  p_box_cyto_order,
  file.path(OUT_DIR, "Fig08b_CytoTRACE2_potency_score_original_order_boxplot_by_celltype"),
  width = 11.2,
  height = 6.2
)

## Optional: also save the raw CytoTRACE2 plotData boxplot with minimal styling.
if (!is.null(cytotrace2_plots) && !is.null(cytotrace2_plots$CytoTRACE2_Boxplot_byPheno)) {
  p_raw_cytotrace2_box <- cytotrace2_plots$CytoTRACE2_Boxplot_byPheno +
    ggtitle("CytoTRACE2 Potency Score by Cell Type") +
    xlab("Cell Type") +
    ylab("CytoTRACE2 Potency Score") +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 9.2),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 14),
      legend.position = "none"
    )
  save_plot_pub(
    p_raw_cytotrace2_box,
    file.path(OUT_DIR, "Fig08c_CytoTRACE2_raw_package_boxplot_by_celltype"),
    width = 11.2,
    height = 6.2
  )
}

############################################################
## Fig09: violin + boxplot by cell type
## Revised version:
## 1) Cell types ordered by mean CytoTRACE2 score, high to low
## 2) Violin colors mapped to dominant CytoTRACE2 potency category
##    rather than cell-type annotation colors
############################################################

## 计算每个细胞类型的平均/中位 CytoTRACE2 score，并按平均分从高到低排序
celltype_score_order_df <- plot_df %>%
  group_by(celltype_short) %>%
  summarise(
    mean_score = mean(CytoTRACE2_Score, na.rm = TRUE),
    median_score = median(CytoTRACE2_Score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

celltype_order_by_score <- as.character(celltype_score_order_df$celltype_short)

## 计算每个细胞类型中占比最高的 potency category
celltype_dominant_potency <- plot_df %>%
  mutate(
    celltype_short = as.character(celltype_short),
    Potency_Category = as.character(Potency_Category)
  ) %>%
  count(celltype_short, Potency_Category, name = "n") %>%
  group_by(celltype_short) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(celltype_short, desc(prop), desc(n)) %>%
  slice(1) %>%
  ungroup() %>%
  select(celltype_short, dominant_potency = Potency_Category, dominant_prop = prop)

## 整理绘图数据
plot_df_violin_potency <- plot_df %>%
  mutate(
    celltype_short = as.character(celltype_short),
    Potency_Category = as.character(Potency_Category)
  ) %>%
  left_join(celltype_dominant_potency, by = "celltype_short") %>%
  mutate(
    celltype_short = factor(celltype_short, levels = celltype_order_by_score),
    dominant_potency = factor(
      dominant_potency,
      levels = c("Multipotent", "Oligopotent", "Unipotent", "Differentiated", "Unknown")
    )
  )

## 为每个细胞类型分配其 dominant potency 对应的颜色
celltype_potency_fill_cols <- celltype_dominant_potency %>%
  mutate(
    celltype_short = factor(celltype_short, levels = celltype_order_by_score),
    fill_col = potency_category_colors[dominant_potency]
  ) %>%
  arrange(celltype_short)

celltype_potency_fill_cols <- setNames(
  celltype_potency_fill_cols$fill_col,
  as.character(celltype_potency_fill_cols$celltype_short)
)

## 保存排序和 dominant potency 信息，便于论文补充表或检查
write.csv(
  celltype_score_order_df,
  file.path(OUT_DIR, "Tab02c_CytoTRACE2_celltype_order_by_mean_score.csv"),
  row.names = FALSE
)

write.csv(
  celltype_dominant_potency,
  file.path(OUT_DIR, "Tab02d_CytoTRACE2_celltype_dominant_potency.csv"),
  row.names = FALSE
)

p_violin_celltype <- ggplot(
  plot_df_violin_potency,
  aes(x = celltype_short, y = CytoTRACE2_Score, fill = celltype_short)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "grey25",
    linewidth = 0.25,
    alpha = 0.88
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    linewidth = 0.25,
    color = "black",
    fill = "white",
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = celltype_potency_fill_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.03)),
    sec.axis = sec_axis(
      trans = ~ .,
      breaks = c(0.08, 0.25, 0.43, 0.60, 0.75, 0.92),
      labels = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"),
      name = NULL
    )
  ) +
  labs(
    x = NULL,
    y = "CytoTRACE2 score"
  ) +
  theme_sc_pub(base_size = 12) +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black", size = 9.5),
    axis.ticks.y.right = element_line(color = "black"),
    axis.line.y.right = element_line(color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.margin = margin(8, 18, 8, 8)
  )

save_plot_pub(
  p_violin_celltype,
  file.path(OUT_DIR, "Fig09_CytoTRACE2_score_violin_by_celltype"),
  width = 10.8,
  height = 5.5
)
## Fig10: group-level score distribution
if (length(disease_levels) > 1) {
  p_box_group <- ggplot(plot_df, aes(x = disease, y = CytoTRACE2_Score, fill = disease)) +
    geom_violin(width = 0.75, trim = TRUE, color = "grey25", linewidth = 0.25, alpha = 0.70) +
    geom_boxplot(width = 0.22, outlier.size = 0.20, linewidth = 0.32, color = "black", fill = "white") +
    scale_fill_manual(values = group_cols, guide = "none") +
    labs(x = NULL, y = "CytoTRACE2 score") +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_blank(), axis.text.x = element_text(face = "bold"))
  save_plot_pub(p_box_group, file.path(OUT_DIR, "Fig10_CytoTRACE2_score_by_disease_status"), width = 4.8, height = 5.2)
  
  p_box_cell_group <- ggplot(plot_df, aes(x = celltype_short, y = CytoTRACE2_Score, fill = disease)) +
    geom_boxplot(width = 0.60, outlier.size = 0.10, linewidth = 0.28, position = position_dodge(width = 0.72)) +
    scale_fill_manual(values = group_cols, name = "Disease status") +
    labs(x = NULL, y = "CytoTRACE2 score") +
    theme_sc_pub(base_size = 12) +
    theme(
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.0),
      legend.key.size = unit(0.34, "cm"),
      legend.text = element_text(size = 9.5)
    )
  save_plot_pub(p_box_cell_group, file.path(OUT_DIR, "Fig11_CytoTRACE2_score_by_celltype_and_disease"), width = 11.8, height = 5.8)
}

############################################################
## 6. Sample-level heatmap and potency composition
############################################################
sample_heat_df <- sample_summary %>%
  mutate(celltype_short = factor(as.character(celltype_short), levels = rev(celltype_levels)))

p_sample_heat <- ggplot(sample_heat_df, aes(x = sample_id, y = celltype_short, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_gradientn(colors = cytotrace_score_cols, name = "Mean\nscore") +
  labs(x = "Sample", y = "Cell type") +
  theme_sc_pub(base_size = 12) +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height = unit(1.2, "cm")
  )
save_plot_pub(p_sample_heat, file.path(OUT_DIR, "Fig12_CytoTRACE2_sample_celltype_mean_score_heatmap"),
              width = max(7.0, 0.42 * length(unique(plot_df$sample_id)) + 4.0),
              height = max(5.8, 0.24 * length(celltype_levels) + 2.0))

potency_celltype_df <- plot_df %>%
  count(celltype_short, Potency_Category, name = "CellNumber") %>%
  group_by(celltype_short) %>%
  mutate(Percent = CellNumber / sum(CellNumber) * 100) %>%
  ungroup()
write.csv(potency_celltype_df, file.path(OUT_DIR, "Tab04_Potency_category_by_celltype.csv"), row.names = FALSE)

p_potency_prop <- ggplot(potency_celltype_df,
                         aes(x = celltype_short, y = Percent / 100, fill = Potency_Category)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.16) +
  scale_fill_manual(values = potency_category_colors, drop = FALSE, name = "Potency") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cell proportion") +
  theme_sc_pub(base_size = 12) +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 9.5)
  )
save_plot_pub(p_potency_prop, file.path(OUT_DIR, "Fig13_Potency_category_composition_by_celltype"), width = 10.8, height = 5.5)

############################################################
## 7. HSC-like focused plots
############################################################
if ("HSC-like" %in% as.character(plot_df$celltype_short)) {
  hsc_df <- plot_df %>% filter(as.character(celltype_short) == "HSC-like")
  write.csv(hsc_df, file.path(OUT_DIR, "Tab05_HSC_like_CytoTRACE2_metadata.csv"), row.names = FALSE)
  
  p_hsc_group <- ggplot(hsc_df, aes(x = disease, y = CytoTRACE2_Score, fill = disease)) +
    geom_violin(width = 0.75, trim = TRUE, color = "grey25", linewidth = 0.25, alpha = 0.70) +
    geom_boxplot(width = 0.22, outlier.size = 0.20, linewidth = 0.32, color = "black", fill = "white") +
    scale_fill_manual(values = group_cols, guide = "none") +
    labs(x = NULL, y = "CytoTRACE2 score") +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_blank(), axis.text.x = element_text(face = "bold"))
  save_plot_pub(p_hsc_group, file.path(OUT_DIR, "Fig14_HSC_like_CytoTRACE2_score_by_disease"), width = 4.8, height = 5.2)
  
  hsc_pot_df <- hsc_df %>%
    count(disease, Potency_Category, name = "CellNumber") %>%
    group_by(disease) %>%
    mutate(Percent = CellNumber / sum(CellNumber) * 100) %>%
    ungroup()
  p_hsc_pot <- ggplot(hsc_pot_df, aes(x = disease, y = Percent / 100, fill = Potency_Category)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.16) +
    scale_fill_manual(values = potency_category_colors, drop = FALSE, name = "Potency") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "HSC-like cell proportion") +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_blank(), axis.text.x = element_text(face = "bold"))
  save_plot_pub(p_hsc_pot, file.path(OUT_DIR, "Fig15_HSC_like_potency_composition_by_disease"), width = 5.8, height = 5.2)
}

############################################################
## 8. Save object with CytoTRACE2 metadata
############################################################
#qsave(seurat_obj, file.path(OUT_DIR, "seurat_with_CytoTRACE2_pub_final.qs"))

message("CytoTRACE2 re-plotted figures exported to: ", OUT_DIR)

# END ORIGINAL SCRIPT P15

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P16: 10_SC_enrich.R
# ==============================================================================
# 清空环境
rm(list = ls())
# 加载必要的包
library(Seurat)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(qs)

# 根目录(按需修改)
ANNO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
ENRICH_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\7.enrich_analysis\\new"
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(ENRICH_DIR)) dir.create(ENRICH_DIR, recursive = TRUE)
# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}


setwd(ANNO_DIR)
#选择分辨率进行降维----
scRNA <- qread("2.resolution_0.5_annotation.qs")
setwd(ENRICH_DIR)
# 查看细胞类型分布
print("细胞类型分布：")
print(table(scRNA$celltype))


celltype_markers <- FindAllMarkers(scRNA, 
                                   group.by = "celltype",
                                   logfc.threshold = 0.5, 
                                   min.pct = 0.25,
                                   only.pos = TRUE)

# 保存标记基因结果，便于后续参考
write.csv(celltype_markers, "Tab01_celltype_markers.csv", row.names = FALSE)

# 创建基因-细胞类型分组数据框
gene_groups <- data.frame(gene = celltype_markers$gene, 
                          celltype = celltype_markers$cluster)


# 将基因符号转换为Entrez ID，用于通路分析
gene_id_mapping <- bitr(celltype_markers$gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Hs.eg.db")

# 输出基因映射统计信息
print(paste0("成功映射 ", nrow(gene_id_mapping), " 个基因，共 ", 
             length(unique(celltype_markers$gene)), " 个"))

# 合并基因ID与细胞类型信息
enrichment_data <- merge(gene_id_mapping, gene_groups, 
                         by.x = 'SYMBOL', 
                         by.y = 'gene')



# 按细胞类型进行KEGG通路富集分析
kegg_results <- compareCluster(ENTREZID ~ celltype,
                               data = enrichment_data,
                               fun = "enrichKEGG",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,  # 设置较为宽松的阈值以捕获更多通路
                               qvalueCutoff = 0.1,   
                               organism = "hsa")

# 将结果中的Entrez ID转回基因符号，提高可读性
kegg_results <- setReadable(kegg_results, 
                            OrgDb = "org.Hs.eg.db", 
                            keyType = "ENTREZID")

# 导出完整的富集分析结果
kegg_table <- kegg_results@compareClusterResult
write.csv(kegg_table, "Tab02_celltype_kegg_enrichment.csv", row.names = FALSE)
table(kegg_table$Cluster)

library(dplyr)

# 假设kegg_table已经存在
kegg_table <- kegg_table %>%
  filter(pvalue < 0.05) %>%                             # 筛选p值小于0.05
  filter(!grepl("Human Diseases", category)) %>%        # 去掉Category含有Human Diseases的行
  group_by(Cluster) %>%                                 # 按Cluster分组
  arrange(desc(Count), .by_group = TRUE) %>%            # 每组按Count降序排列
  ungroup()

write.csv(kegg_table, "Tab03_filtered_celltype_kegg_enrichment.csv", row.names = FALSE)


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(stringr)

# 读取数据
data <- kegg_table

# 查看数据结构
# str(data)
# head(data)

# 查看所有可用的细胞类型
# cat("可用的细胞类型：\n")
print(unique(data$celltype))

# 用户选择细胞类型（修改这里来选择您想要的细胞类型）
# 示例：选择前6个细胞类型，用户可以根据需要修改
selected_celltypes <- c("HSC-like", "GMP-like", "CyclingProg", "Progenitor","GranuleProg","MyeloidProg","LymphoidProg")

# 如果您想要使用所有细胞类型，可以用下面这行代替上面的代码：
# selected_celltypes <- unique(data$celltype)

# 过滤数据，只保留选定的细胞类型
filtered_data <- data %>%
  filter(celltype %in% selected_celltypes)

# 数据预处理：为每个celltype选择zScore最高的前15条通路
top_pathways <- filtered_data %>%
  group_by(celltype) %>%
  arrange(desc(zScore)) %>%
  slice_head(n = 15) %>%
  ungroup()

# 为了更好的可视化，截断过长的Description
top_pathways$Description_short <- str_wrap(top_pathways$Description, width = 100)

# 创建颜色调色板
unique_celltypes <- unique(top_pathways$celltype)
n_celltypes <- length(unique_celltypes)

# 使用多样化的颜色
colors <- c("#E8A5A5", "#A5E8A5", "#A5A5E8", "#E8E8A5", "#E8A5E8", 
            "#A5E8E8", "#FFB366", "#66B3FF", "#66FFB3", "#B366FF",
            "#FF66B3", "#B3FF66", "#FF6666", "#66FF66", "#6666FF",
            "#FFA366", "#66FFA3", "#A366FF", "#FF66A3", "#A3FF66")

# 如果细胞类型数量超过预设颜色，使用RColorBrewer生成更多颜色
if(n_celltypes > length(colors)) {
  colors <- rainbow(n_celltypes)
}

# 创建颜色映射
color_mapping <- setNames(colors[1:n_celltypes], unique_celltypes)

# 为每个通路创建唯一的y轴位置
plot_data <- top_pathways %>%
  arrange(celltype, desc(zScore)) %>%
  mutate(
    # 为每个通路创建唯一的y轴标识
    pathway_id = paste0(celltype, "_", row_number()),
    # 创建用于排序的数字
    celltype_num = as.numeric(as.factor(celltype)),
    # 在同一细胞类型内按zScore排序
    within_celltype_rank = ave(zScore, celltype, FUN = function(x) rank(-x))
  ) %>%
  # 按细胞类型和zScore排序
  arrange(celltype_num, desc(zScore))   #可以选择用那一类来升序或降序排列

# 使用ggh4x包来实现不同的strip背景色
# 如果没有安装ggh4x，请先运行：install.packages("ggh4x")
library(ggh4x)

# 创建每个细胞类型对应的背景色列表
strip_colors <- color_mapping[unique_celltypes]
names(strip_colors) <- unique_celltypes

# 创建水平条形图 - 使用ggh4x实现不同背景色
p_modified <- ggplot(plot_data, aes(x = zScore + 1, y = reorder(Description_short, zScore))) +  # 将条形图向右偏移1个单位
  geom_col(aes(fill = celltype), alpha = 0.8, width = 0.8) +  # 增加条形宽度，使字体能被包围在条形里
  # 在颜色条块里添加通路名字，垂直居中，左侧偏移1个单位
  geom_text(aes(label = Description_short), 
            hjust = 0, x = 1.1, size = 6, vjust = 0.5, color = "black") +  # 增大字体到6，让文字更清晰
  scale_fill_manual(values = color_mapping, guide = "none") +
  facet_grid2(celltype ~ ., scales = "free_y", space = "free_y", switch = "y",
              strip = strip_themed(background_y = elem_list_rect(fill = strip_colors))) +  # 使用ggh4x设置不同背景色
  theme_minimal() +
  theme(
    # 隐藏y轴文本，因为文本现在在条形图内部
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),  # 增大x轴字体
    axis.title = element_text(size = 18),  # 增大轴标题字体，不加粗
    plot.title = element_text(size = 22, hjust = 0.5),  # 增大标题字体，不加粗
    # 调整细胞名显示：左对齐，显示在左侧
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 18, 
                                     margin = margin(l = 8, r = 15, t = 8, b = 8),
                                     color = "black"),  # 增大字体，增加边距
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.6, "lines"),  # 增加面板间距，使每个条形图的上下边距增大
    # 调整strip位置到左侧
    strip.placement = "outside",
    # 增加图的整体边距
    plot.margin = margin(30, 30, 30, 30)
  ) +
  labs(
    title = "Enrichment of celltype marker genes",
    x = "Z-Score",
    y = NULL
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  # 调整y轴的扩展，为每个条形图增加上下边距
  scale_y_discrete(expand = expansion(add = c(0.2, 0.2)))  # 为每个条形图增加上下0.2个单位的边距

# 创建一个更大的绘图设备并保存为 PDF
pdf("Fig01_celltype_kegg_enrichment_modified.pdf", width = 15, height = 30)  # 增大图像宽度，让元素更清晰
print(p_modified)
dev.off()




# END ORIGINAL SCRIPT P16

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P17: 11_SC_cellchat.R
# ==============================================================================
rm(list = ls())  # 清空环境变量
library(Seurat)
library(CellChat)
library(future)
library(ggplot2)
library(dplyr)
library(tidyr)
library(qs)


# install.packages("CellChat")
# BiocManager::install("CellChat")
# remotes::install_github("sqjin/CellChat")

# 根目录(按需修改)
ANNO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation\\new"
COMMUNI_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\4.Sc_cell_communication\\final"
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(COMMUNI_DIR)) dir.create(COMMUNI_DIR, recursive = TRUE)
# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

setwd(ANNO_DIR)
# 细胞通讯----
seurat <- qread("2.resolution_0.5_annotation.qs")

DefaultAssay(seurat) <- "RNA"

stopifnot("celltype" %in% colnames(seurat@meta.data))

# 清理 celltype（CellChat 不接受 NA/空）
seurat$celltype <- as.character(seurat$celltype)
seurat <- subset(seurat, subset = !is.na(celltype) & celltype != "")

# NormalizeData 建议在 sequential 下运行（避免 future 大对象传输）
future::plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
Idents(seurat) <- seurat$celltype

cat("Genes x Cells:", nrow(seurat), "x", ncol(seurat), "\n")
cat("Celltype levels:", length(levels(Idents(seurat))), "\n")

# 取 Seurat v5 的 data layer（log-normalized）
data.input <- GetAssayData(seurat, assay = "RNA", layer = "data")

meta <- seurat@meta.data
meta <- meta[match(colnames(data.input), rownames(meta)), , drop = FALSE]
stopifnot(all(colnames(data.input) == rownames(meta)))




setwd(COMMUNI_DIR)
# 创建 CellChat（matrix + meta 接口，兼容 Seurat v5）
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat <- setIdent(cellchat, ident.use = "celltype")
groupSize <- as.numeric(table(cellchat@idents))      #每个细胞类型的细胞数目

#导入配体和受体数据库
CellChatDB <- CellChatDB.human       #如果物种是小鼠改成CellChatDB.mouse

# cellchat <- subsetData(cellchat)
cellchat@DB <- CellChatDB.human  # DB：明确使用 human 全库（不要 subsetDB 空参数）
stopifnot(nrow(cellchat@DB$interaction) > 0)

#查看配体和受体对的类型
pdf(file="Fig01_DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

#对表达数据进行预处理
# subset signaling genes（关键 QC 点）
cellchat <- subsetData(cellchat)
cat("data.signaling dim:", paste(dim(cellchat@data.signaling), collapse = " x "), "\n")
stopifnot(nrow(cellchat@data.signaling) > 0)

options(future.globals.maxSize = 20 * 1024^3)
future::plan("multisession", workers = 10) #设置使用的电脑核心数，并行计算节省时间
cellchat <- identifyOverExpressedGenes(cellchat)      #识别每种细胞中过表达的基因

cellchat <- identifyOverExpressedInteractions(cellchat)      #识别基因的相互作用

#cellchat <- smoothData(cellchat, adj = PPI.human)
#smoothData() 用于根据蛋白互作网络（PPI）平滑基因表达。但是从 CellChat v1.5 起
#该步骤已整合进 subsetData() 或被弃用，因为它在大多数分析中不再需要，且可能造成重复计算。

#计算细胞通讯的概率
cellchat <- computeCommunProb(cellchat, population.size = TRUE)

#过滤掉小于10个细胞的细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#输出细胞间的通讯关系
df.net=subsetCommunication(cellchat)
library(openxlsx)
write.xlsx(df.net, "Tab01_network.xlsx")

#在信号通路的水平进一步推测胞间的通讯, 推断通路水平的互作网络
cellchat <- computeCommunProbPathway(cellchat)

#对计算结果汇总整合，展示整体细胞通讯状态
cellchat <- aggregateNet(cellchat)
#输出细胞通讯的图形(互作数量的图形)
pdf(file="Fig02_cellNetworkCount.pdf", width=10, height=12)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

#输出细胞通讯的图形(互作强度的图形)
pdf(file="Fig03_cellNetworkWeight.pdf", width=10, height=12)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

# 热图展示细胞群体之间的通讯强度
# 生成热图并保存为 PDF 文件
pdf(file = "Fig04_cellCommunicationHeatmap.pdf", width = 10, height = 8)  # 设置PDF文件名及大小
netVisual_heatmap(cellchat, measure = "weight")  # 生成热图
dev.off()  # 关闭 PDF 输出设备
# 提取相关性矩阵并保存为 Excel 文件
correlation_matrix <- cellchat@net$weight  # 提取相关性矩阵
# 将相关性矩阵保存为 Excel 文件
write.xlsx(correlation_matrix, file = "Tab02_cellCommunicationCorrelation.xlsx", row.names = TRUE)


#分细胞类型展示(把单个细胞类型提取出来,观察这个细胞与其他细胞的通讯)
pdf(file="Fig05_singleCell.pdf", width=18, height=15)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=1,title.name=cel)
}
dev.off()


# 计算网络中心度
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# netP 是 CellChat 包中的一个数据槽（slot），它包含了细胞通讯网络的“概率”信息。
# 在 netAnalysis_computeCentrality 函数中，netP 表示细胞通讯网络中每个细胞在通信中的中心性度量。
# 这有助于评估每个细胞在整个通讯网络中的重要性
# 更高的中心度通常意味着该细胞在网络中承担更关键的角色，是一个重要的“信号传播者”

# 生成并保存调整大小后的 PDF 文件
pdf(file = "Fig06_signalingRoleHeatmap.pdf", width = 10, height = 6)  # 设置 PDF 宽度和高度

# 绘制并合并两个热图
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# 合并并显示热图
ht1 + ht2
# 关闭 PDF 输出设备
dev.off()

# outgoing：表示“外向”信号，指的是一个细胞向其他细胞发送信号。
# incoming：表示“内向”信号，指的是其他细胞向该细胞发送信号。
# 在细胞通讯分析中，这两个方向很重要，因为它们有助于理解一个细胞在网络中的“角色”。
# 例如，一个细胞可能是信号的接收者（incoming），也可能是信号的发送者（outgoing）。



# 绘制散点图并添加斜向 45 度虚线
gg1 <- netAnalysis_signalingRole_scatter(cellchat) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 添加斜向 45 度虚线
  labs(x = "Outgoing interaction strength", y = "Incoming interaction strength")  # 添加坐标轴标签

# 保存 gg1 为 PDF 文件
ggsave("Fig07_signalingRoleScatter.pdf", plot = gg1, width = 8, height = 6)




# 评估通讯模式----
# outgoing
library(NMF)
library(ggalluvial)
# 生成图形并保存为 PDF 文件
pdf(file = "Fig08_selectK_outgoing.pdf", width = 10, height = 8)  # 设置 PDF 文件名及大小
# 绘制并显示图形
selectK(cellchat, pattern = "outgoing")
# 关闭 PDF 输出设备
dev.off()
# selectK() 函数的主要作用是根据不同的 k 值（即细胞通讯模式的数量）计算科菲尼系数和轮廓系数，并绘制曲线，
# 帮助选择最佳的 k 值。通过查看这些曲线的变化，你可以看到科菲尼系数和轮廓系数在某些 k 值下出现剧烈的变化，
# 进而推测出最佳的 k 值。

# 当外向模式数量为 6 时，科菲尼系数和轮廓系数都突然下降
nPatterns = 5
# 生成并保存为 PDF 文件
pdf(file = "Fig09_outgoing_communicationPatterns.pdf", width = 10, height = 5)  # 设置 PDF 文件名和大小
# 调用 identifyCommunicationPatterns 函数并绘制图形
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# 关闭 PDF 输出设备
dev.off()

# 在 identifyCommunicationPatterns 函数中，k 代表需要提取的通讯模式的数量，即要识别的“信号模式”的个数。
# nPatterns 是指你希望从细胞通讯网络中识别的模式数目。
# 通常，选择合适的 k 值是通过分析模型性能、科菲尼系数（Cohen's Kappa）和轮廓系数（Silhouette Coefficient）来确定的。
# # 选择合适的 k 值帮助我们识别最具生物学意义的细胞通讯模式。

# 桑葚图
# 生成并保存桑葚图为 PDF 文件
pdf(file = "Fig10_outgoing_riverPlot.pdf", width = 12, height = 8)  # 设置 PDF 文件名和大小
# 绘制桑葚图
netAnalysis_river(cellchat, pattern = "outgoing")
# 关闭 PDF 输出设备
dev.off()

#点图
# 生成并保存点图为 PDF 文件
pdf(file = "Fig11_outgoing_dotPlot.pdf", width = 12, height = 8)  # 设置 PDF 文件名和大小
# 绘制点图
netAnalysis_dot(cellchat, pattern = "outgoing")
# 关闭 PDF 输出设备
dev.off()






# incoming
library(NMF)
library(ggalluvial)
# 生成图形并保存为 PDF 文件
pdf(file = "Fig08_selectK_incoming.pdf", width = 10, height = 8)  # 设置 PDF 文件名及大小
# 绘制并显示图形
selectK(cellchat, pattern = "incoming")
# 关闭 PDF 输出设备
dev.off()

# 科菲尼系数和轮廓系数都突然下降的前一个值
nPatterns = 3
# 生成并保存为 PDF 文件
pdf(file = "Fig09_incoming_communicationPatterns.pdf", width = 10, height = 5)  # 设置 PDF 文件名和大小
# 调用 identifyCommunicationPatterns 函数并绘制图形
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# 关闭 PDF 输出设备
dev.off()

# 桑葚图
# 生成并保存桑葚图为 PDF 文件
pdf(file = "Fig10_incoming_riverPlot.pdf", width = 12, height = 8)  # 设置 PDF 文件名和大小
# 绘制桑葚图
netAnalysis_river(cellchat, pattern = "incoming")
# 关闭 PDF 输出设备
dev.off()

#点图
# 生成并保存点图为 PDF 文件
pdf(file = "Fig11_incoming_dotPlot.pdf", width = 12, height = 8)  # 设置 PDF 文件名和大小
# 绘制点图
netAnalysis_dot(cellchat, pattern = "incoming")
# 关闭 PDF 输出设备
dev.off()


############################################################
# 保存完整 CellChat 对象，供后续统一风格重绘使用
############################################################
qsave(cellchat, file.path(COMMUNI_DIR, "cellchat_object.qs"))


# END ORIGINAL SCRIPT P17

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P18: 11_SC_cellchat/plot_cellchat_pub.R
# ==============================================================================
############################################################
## 03_plot_cellchat_pub_hsclike_final.R
## Publication-style CellChat plotting
## Focus: HSC-like communication and LR/signaling patterns
## Palette strictly follows 00_pub_theme_colors.R
############################################################

suppressPackageStartupMessages({
  library(qs)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(grid)
  library(patchwork)
})

source("00_pub_theme_colors.R")

############################################################
## 0. Output paths and basic parameters
############################################################
OUT_DIR <- file.path(COMMUNI_DIR, "pub_figures_cellchat")
HSC_DIR <- file.path(OUT_DIR, "HSC_like_specific")
safe_mkdir(OUT_DIR)
safe_mkdir(HSC_DIR)

HSC_CELL <- "HSC-like"
TOP_N_PATHWAYS <- 15
TOP_N_LR <- 35
PVALUE_CUTOFF <- 0.05
MIN_PROB_FOR_LR <- 0

## Circle/network plot settings:
## keep label size modest, only make labels bold to avoid GMP-like overlap.
CIRCLE_LABEL_CEX <- 0.86
CIRCLE_LABEL_FONT <- 2
PATHWAY_LABEL_CEX <- 0.68
UNIFORM_VERTEX_WEIGHT <- 8
UNIFORM_VERTEX_SIZE_MAX <- 6

############################################################
## 1. Compatibility helpers
############################################################
`%||%` <- function(a, b) if (!is.null(a)) a else b

save_base_plot <- function(filename_base, width = 7, height = 6, expr, dpi = 450) {
  safe_mkdir(dirname(filename_base))

  ## PDF. Use grDevices::pdf rather than cairo_pdf for CellChat/igraph base plots.
  grDevices::pdf(paste0(filename_base, ".pdf"), width = width, height = height, useDingbats = FALSE)
  tryCatch(
    force(expr),
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Plot failed:", conditionMessage(e)), cex = 0.8)
    },
    finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    }
  )

  ## PNG. Re-open and redraw, because base graphics cannot be copied reliably.
  grDevices::png(paste0(filename_base, ".png"), width = width * dpi, height = height * dpi, res = dpi)
  tryCatch(
    force(expr),
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Plot failed:", conditionMessage(e)), cex = 0.8)
    },
    finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    }
  )
  invisible(filename_base)
}

save_plot_pub_wide <- function(plot, filename_base, width = 20, height = 8, dpi = 450, bg = "white") {
  safe_mkdir(dirname(filename_base))
  ggplot2::ggsave(
    paste0(filename_base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    bg = bg,
    useDingbats = FALSE,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(filename_base, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  invisible(filename_base)
}

make_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x
}

shorten_celltype_name <- function(x) {
  x <- as.character(x)
  out <- x
  if (exists("celltype_map")) {
    mapped <- unname(celltype_map[out])
    out[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }
  out <- dplyr::case_when(
    out %in% c("HSC-like cell", "HSC_like", "HSC-like-cell") ~ "HSC-like",
    out %in% c("GMP-like cell") ~ "GMP-like",
    out %in% c("Progenitor-like cell") ~ "Progenitor",
    out %in% c("Cycling Progenitor") ~ "CyclingProg",
    out %in% c("LMPP-like cell") ~ "LymphoidProg",
    out %in% c("MBP-like cell") ~ "GranuleProg",
    out %in% c("MDP-like cell") ~ "MyeloidProg",
    out %in% c("ProMyelocyte") ~ "Promyelo",
    out %in% c("Naive T cell") ~ "NaiveT",
    out %in% c("NK cell") ~ "NK_T",
    out %in% c("B cell") ~ "Bcell",
    out %in% c("Plasma cell") ~ "Plasma",
    out %in% c("Mast cell") ~ "Mast",
    TRUE ~ out
  )
  out
}

rename_dimnames_short <- function(x, map_fun = shorten_celltype_name) {
  if (is.null(dim(x)) || length(dim(x)) < 2) return(x)
  dn <- dimnames(x)
  if (length(dn) >= 1 && !is.null(dn[[1]])) dn[[1]] <- map_fun(dn[[1]])
  if (length(dn) >= 2 && !is.null(dn[[2]])) dn[[2]] <- map_fun(dn[[2]])
  dimnames(x) <- dn
  x
}

standardize_cellchat_names <- function(cellchat) {
  old_levels <- levels(cellchat@idents)
  if (is.null(old_levels)) old_levels <- unique(as.character(cellchat@idents))
  new_idents <- shorten_celltype_name(as.character(cellchat@idents))
  new_levels <- c(intersect(celltype_levels_short, unique(new_idents)), setdiff(unique(new_idents), celltype_levels_short))
  cellchat@idents <- factor(new_idents, levels = new_levels)

  for (nm in names(cellchat@net)) {
    cellchat@net[[nm]] <- rename_dimnames_short(cellchat@net[[nm]])
  }
  for (nm in names(cellchat@netP)) {
    cellchat@netP[[nm]] <- rename_dimnames_short(cellchat@netP[[nm]])
  }
  cellchat
}

get_ordered_levels <- function(cellchat) {
  lev <- levels(cellchat@idents)
  if (is.null(lev)) lev <- unique(as.character(cellchat@idents))
  if (!is.null(cellchat@net$weight) && !is.null(dimnames(cellchat@net$weight)[[1]])) {
    lev <- dimnames(cellchat@net$weight)[[1]]
  }
  c(intersect(celltype_levels_short, lev), setdiff(lev, celltype_levels_short))
}

get_ordered_matrix <- function(mat, levels_use) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat))) rownames(mat) <- levels_use[seq_len(nrow(mat))]
  if (is.null(colnames(mat))) colnames(mat) <- levels_use[seq_len(ncol(mat))]
  keep <- intersect(levels_use, intersect(rownames(mat), colnames(mat)))
  mat <- mat[keep, keep, drop = FALSE]
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  mat
}

get_cell_cols_for <- function(levels_use) {
  get_celltype_colors(levels_use, use_short = TRUE)
}

solid_size_legend <- function(title_text = "Communication probability") {
  guide_legend(
    title = title_text,
    override.aes = list(shape = 16, color = "black", fill = "black", alpha = 1)
  )
}

extract_lr_label <- function(df) {
  if (all(c("ligand", "receptor") %in% colnames(df))) {
    return(paste(df$ligand, df$receptor, sep = "-"))
  }
  if ("interaction_name_2" %in% colnames(df)) {
    return(gsub(" \\+ ", "+", gsub(" - ", "-", as.character(df$interaction_name_2))))
  }
  if ("interaction_name" %in% colnames(df)) {
    return(gsub("_", "-", as.character(df$interaction_name)))
  }
  paste0("LR_", seq_len(nrow(df)))
}

build_pathway_lr_label_map <- function(comm_df) {
  if (is.null(comm_df) || nrow(comm_df) == 0 || !"pathway_name" %in% colnames(comm_df)) {
    return(character(0))
  }
  df <- comm_df
  if (!"prob" %in% colnames(df)) df$prob <- 1
  df$LR_pair <- extract_lr_label(df)
  df$receptor_part <- ifelse(grepl("-", df$LR_pair), sub("^[^-]+-", "", df$LR_pair), NA_character_)
  df %>%
    filter(!is.na(pathway_name), pathway_name != "") %>%
    arrange(desc(prob)) %>%
    group_by(pathway_name) %>%
    summarise(receptor = dplyr::first(na.omit(receptor_part)), .groups = "drop") %>%
    mutate(
      receptor = gsub("_", "+", receptor),
      label = ifelse(is.na(receptor) | receptor == "", pathway_name, paste0(pathway_name, "-", receptor))
    ) %>%
    {setNames(.$label, .$pathway_name)}
}

############################################################
## 2. Read CellChat object and standardize names/colors
############################################################
candidate_files <- c(
  file.path(COMMUNI_DIR, "cellchat_object.qs"),
  file.path(COMMUNI_DIR, "cellchat.qs"),
  file.path(COMMUNI_DIR, "CellChat_object.qs"),
  file.path(PROJECT_DIR, "4.Sc_cell_communication", "new", "cellchat_object.qs"),
  file.path(PROJECT_DIR, "4.Sc_cell_communication", "final", "cellchat_object.qs"),
  file.path(PROJECT_DIR, "4.Sc_cell_communication", "cellchat_object.qs")
)
cellchat_file <- candidate_files[file.exists(candidate_files)][1]

if (is.na(cellchat_file) || length(cellchat_file) == 0) {
  stop("Cannot find CellChat object. Expected file name: cellchat_object.qs")
}

cellchat <- qread(cellchat_file)
cellchat <- standardize_cellchat_names(cellchat)

idents <- get_ordered_levels(cellchat)
cell_cols <- get_cell_cols_for(idents)

weight_mat <- get_ordered_matrix(cellchat@net$weight, idents)
count_mat  <- get_ordered_matrix(cellchat@net$count, idents)
idents <- rownames(weight_mat)
cell_cols <- get_cell_cols_for(idents)

object_levels <- levels(cellchat@idents)
object_cols <- get_cell_cols_for(object_levels)

hsc_name <- if (HSC_CELL %in% idents) HSC_CELL else grep("HSC", idents, ignore.case = TRUE, value = TRUE)[1]
if (is.na(hsc_name) || length(hsc_name) == 0) {
  stop("Cannot find HSC-like identity. Existing CellChat groups: ", paste(idents, collapse = ", "))
}
hsc_idx <- match(hsc_name, idents)

vertex_uniform <- rep(UNIFORM_VERTEX_WEIGHT, length(idents))
names(vertex_uniform) <- idents

write.csv(data.frame(celltype = idents, color = unname(cell_cols)),
          file.path(OUT_DIR, "Tab_cellchat_celltype_colors.csv"), row.names = FALSE)
write.csv(data.frame(HSC_like_identity = hsc_name, HSC_like_index = hsc_idx),
          file.path(HSC_DIR, "Tab_HSC_like_identity_used.csv"), row.names = FALSE)

############################################################
## 3. Circle plot helper
############################################################
plot_circle_matrix <- function(mat, title_text = NULL, edge_max = NULL) {
  mat <- get_ordered_matrix(mat, idents)
  cols <- get_cell_cols_for(rownames(mat))
  vw <- rep(UNIFORM_VERTEX_WEIGHT, nrow(mat))
  names(vw) <- rownames(mat)
  if (is.null(edge_max)) edge_max <- max(mat, na.rm = TRUE)
  if (!is.finite(edge_max) || edge_max <= 0) edge_max <- 1

  par(mfrow = c(1, 1), mar = c(1.2, 1.2, 4.2, 1.2), xpd = NA)
  
  netVisual_circle(
    mat,
    vertex.weight = vw,
    vertex.weight.max = UNIFORM_VERTEX_WEIGHT,
    vertex.size.max = UNIFORM_VERTEX_SIZE_MAX,
    weight.scale = TRUE,
    edge.weight.max = edge_max,
    edge.width.max = 6,
    label.edge = FALSE,
    color.use = cols,
    vertex.label.cex = CIRCLE_LABEL_CEX,
    vertex.label.color = "black",
    title.name = ""
  )
  
  if (!is.null(title_text) && title_text != "") {
    mtext(title_text, side = 3, line = 1.2, cex = 1.15, font = 2)
  }
}

############################################################
## 4. General CellChat overview plots
############################################################
safe_run("Fig01 overall interaction number circle", {
  save_base_plot(file.path(OUT_DIR, "Fig01_overall_interaction_number_circle"), 9, 8, {
    plot_circle_matrix(count_mat, "Number of interactions", edge_max = max(count_mat, na.rm = TRUE))
  })
})

safe_run("Fig02 overall interaction strength circle", {
  save_base_plot(file.path(OUT_DIR, "Fig02_overall_interaction_strength_circle"), 9, 8, {
    plot_circle_matrix(weight_mat, "Interaction strength", edge_max = max(weight_mat, na.rm = TRUE))
  })
})

safe_run("Fig03 overall interaction strength heatmap", {
  save_base_plot(file.path(OUT_DIR, "Fig03_overall_interaction_strength_heatmap"), 8.5, 7.2, {
    ht <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
    print(ht)
  })
})

safe_run("Fig04 overall signaling role scatter", {
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  p <- netAnalysis_signalingRole_scatter(cellchat) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "black") +
    labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55))
  save_plot_pub(p, file.path(OUT_DIR, "Fig04_overall_signaling_role_scatter"), 7.2, 6.4)
})

safe_run("Fig05 overall signaling role heatmap", {
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  save_base_plot(file.path(OUT_DIR, "Fig05_overall_signaling_role_heatmap_out_in"), 12, 6.5, {
    ht_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.heatmap = "Reds")
    ht_in  <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "Blues")
    print(ht_out + ht_in)
  })
})

############################################################
## 5. Fig06 full source.target LR bubble
## Keep full source.target categories and use very wide canvas.
############################################################
safe_run("Fig06 full source.target LR bubble", {
  comm_all <- subsetCommunication(cellchat)
  if (is.null(comm_all) || nrow(comm_all) == 0) stop("No communication table found")
  if (!"prob" %in% colnames(comm_all)) comm_all$prob <- 1
  if (!"pval" %in% colnames(comm_all)) comm_all$pval <- 1
  if (!"pathway_name" %in% colnames(comm_all)) comm_all$pathway_name <- "Unknown"

  comm_all$source <- shorten_celltype_name(comm_all$source)
  comm_all$target <- shorten_celltype_name(comm_all$target)
  comm_all$LR_pair <- extract_lr_label(comm_all)

  df <- comm_all %>%
    filter(prob > MIN_PROB_FOR_LR) %>%
    filter(is.na(pval) | pval <= PVALUE_CUTOFF) %>%
    mutate(
      source = factor(source, levels = idents),
      target = factor(target, levels = idents),
      source.target = paste(source, target, sep = " -> "),
      LR_pair = gsub("_", "-", LR_pair)
    ) %>%
    arrange(source, target, desc(prob)) %>%
    mutate(
      source.target = factor(source.target, levels = unique(source.target)),
      LR_pair = factor(LR_pair, levels = rev(unique(LR_pair)))
    )

  if (nrow(df) == 0) {
    df <- comm_all %>%
      filter(prob > MIN_PROB_FOR_LR) %>%
      mutate(
        source = factor(shorten_celltype_name(source), levels = idents),
        target = factor(shorten_celltype_name(target), levels = idents),
        source.target = paste(source, target, sep = " -> "),
        LR_pair = factor(gsub("_", "-", LR_pair))
      ) %>%
      arrange(source, target, desc(prob)) %>%
      mutate(source.target = factor(source.target, levels = unique(source.target)))
  }

  if (nrow(df) == 0) stop("No communications available for Fig06")

  df$neg_log10_p <- -log10(pmax(df$pval, 1e-300))
  n_x <- length(levels(df$source.target))
  n_y <- length(levels(df$LR_pair))
  fig_w <- max(36, 0.42 * n_x + 8)
  fig_h <- max(9, min(65, 0.22 * n_y + 4.2))

  p <- ggplot(df, aes(x = source.target, y = LR_pair)) +
    geom_point(aes(size = neg_log10_p, fill = prob),
               shape = 21, color = "black", stroke = 0.16, alpha = 0.92) +
    scale_fill_gradientn(
      colors = c("#2C7BB6", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C"),
      name = "Communication\nprobability"
    ) +
    scale_size_continuous(name = "-log10(p)", range = c(1.0, 4.5)) +
    labs(x = "Source -> target cell type", y = "Ligand-receptor pair") +
    theme_sc_pub(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 5.0, color = "black"),
      axis.text.y = element_text(size = 7.4, color = "black"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.45),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.18),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.30, "cm"),
      legend.text = element_text(size = 8.0)
    )

  save_plot_pub_wide(p, file.path(OUT_DIR, "Fig06_overall_full_source_target_LR_bubble"), fig_w, fig_h)
})

############################################################
## 6. HSC-like related communication scores and pathway tables
############################################################
hsc_out_weight <- weight_mat[hsc_name, , drop = TRUE]
hsc_in_weight  <- weight_mat[, hsc_name, drop = TRUE]
hsc_out_count  <- count_mat[hsc_name, , drop = TRUE]
hsc_in_count   <- count_mat[, hsc_name, drop = TRUE]

hsc_score_df <- bind_rows(
  data.frame(direction = "Outgoing from HSC-like", partner = names(hsc_out_weight), weight = as.numeric(hsc_out_weight), count = as.numeric(hsc_out_count)),
  data.frame(direction = "Incoming to HSC-like", partner = names(hsc_in_weight), weight = as.numeric(hsc_in_weight), count = as.numeric(hsc_in_count))
) %>% filter(partner != hsc_name | weight > 0 | count > 0)
write.csv(hsc_score_df, file.path(HSC_DIR, "Tab_HSC_like_partner_interaction_scores.csv"), row.names = FALSE)

get_hsc_pathway_scores <- function(cellchat, hsc_name) {
  pathways <- cellchat@netP$pathways %||% dimnames(cellchat@netP$prob)[[3]]
  prob <- cellchat@netP$prob
  if (is.null(pathways) || is.null(prob) || length(dim(prob)) < 3) return(NULL)
  dn <- dimnames(prob)
  if (is.null(dn[[1]])) dn[[1]] <- levels(cellchat@idents)
  if (is.null(dn[[2]])) dn[[2]] <- levels(cellchat@idents)
  if (is.null(dn[[3]])) dn[[3]] <- pathways
  dn[[1]] <- shorten_celltype_name(dn[[1]])
  dn[[2]] <- shorten_celltype_name(dn[[2]])
  dimnames(prob) <- dn
  if (!hsc_name %in% dn[[1]] || !hsc_name %in% dn[[2]]) return(NULL)
  out_score <- apply(prob[hsc_name, , , drop = FALSE], 3, sum, na.rm = TRUE)
  in_score  <- apply(prob[, hsc_name, , drop = FALSE], 3, sum, na.rm = TRUE)
  data.frame(
    pathway = names(out_score),
    outgoing_strength = as.numeric(out_score),
    incoming_strength = as.numeric(in_score),
    total_HSC_related_strength = as.numeric(out_score) + as.numeric(in_score),
    stringsAsFactors = FALSE
  )
}

hsc_pathway_df <- get_hsc_pathway_scores(cellchat, hsc_name)
if (!is.null(hsc_pathway_df)) {
  hsc_pathway_df <- hsc_pathway_df %>% arrange(desc(total_HSC_related_strength))
  write.csv(hsc_pathway_df, file.path(HSC_DIR, "Tab_HSC_like_related_pathway_strength.csv"), row.names = FALSE)
  pathways_hsc_out <- hsc_pathway_df %>% filter(outgoing_strength > 0) %>% arrange(desc(outgoing_strength)) %>% pull(pathway)
  pathways_hsc_in  <- hsc_pathway_df %>% filter(incoming_strength > 0) %>% arrange(desc(incoming_strength)) %>% pull(pathway)
  pathways_hsc_any <- hsc_pathway_df %>% filter(total_HSC_related_strength > 0) %>% pull(pathway)
} else {
  pathways_hsc_out <- pathways_hsc_in <- pathways_hsc_any <- character(0)
}

############################################################
## 7. HSC-like circle networks
############################################################
safe_run("HSC01 HSC-like outgoing circle", {
  mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  mat[hsc_name, ] <- weight_mat[hsc_name, ]
  save_base_plot(file.path(HSC_DIR, "HSC01_HSC_like_outgoing_circle"), 8.5, 7.8, {
    plot_circle_matrix(mat, "HSC-like outgoing signaling", edge_max = max(weight_mat, na.rm = TRUE))
  })
})

safe_run("HSC02 HSC-like incoming circle", {
  mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  mat[, hsc_name] <- weight_mat[, hsc_name]
  save_base_plot(file.path(HSC_DIR, "HSC02_HSC_like_incoming_circle"), 8.5, 7.8, {
    plot_circle_matrix(mat, "HSC-like incoming signaling", edge_max = max(weight_mat, na.rm = TRUE))
  })
})

safe_run("HSC03 HSC-like bidirectional circle", {
  mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  mat[hsc_name, ] <- weight_mat[hsc_name, ]
  mat[, hsc_name] <- weight_mat[, hsc_name]
  save_base_plot(file.path(HSC_DIR, "HSC03_HSC_like_bidirectional_circle"), 8.5, 7.8, {
    plot_circle_matrix(mat, "HSC-like bidirectional signaling", edge_max = max(weight_mat, na.rm = TRUE))
  })
})

############################################################
## 8. HSC-like partner scores
############################################################
safe_run("HSC04 HSC-like partner scores", {
  df <- hsc_score_df %>% filter(weight > 0 | count > 0)
  if (nrow(df) == 0) stop("No HSC-related partner interactions")
  p <- ggplot(df, aes(x = fct_reorder(partner, weight), y = weight, fill = direction)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "black", linewidth = 0.18) +
    coord_flip() +
    scale_fill_manual(values = c("Outgoing from HSC-like" = "#1F4E79", "Incoming to HSC-like" = "#8C6D31")) +
    labs(x = NULL, y = "Communication strength", fill = NULL) +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55))
  save_plot_pub(p, file.path(HSC_DIR, "HSC04_HSC_like_partner_communication_strength"), 7.2, 6.2)
})

############################################################
## 9. HSC-like LR bubble plots
############################################################
comm_all <- tryCatch(subsetCommunication(cellchat), error = function(e) NULL)
if (!is.null(comm_all) && nrow(comm_all) > 0) {
  if (!"prob" %in% colnames(comm_all)) comm_all$prob <- 1
  if (!"pval" %in% colnames(comm_all)) comm_all$pval <- 1
  if (!"pathway_name" %in% colnames(comm_all)) comm_all$pathway_name <- "Unknown"
  comm_all$source <- shorten_celltype_name(comm_all$source)
  comm_all$target <- shorten_celltype_name(comm_all$target)
  comm_all$LR_pair <- extract_lr_label(comm_all)
  comm_all$LR_pair <- gsub("_", "-", comm_all$LR_pair)
}

comm_out <- if (!is.null(comm_all)) comm_all %>% filter(source == hsc_name, prob > MIN_PROB_FOR_LR, pval <= PVALUE_CUTOFF) else NULL
comm_in  <- if (!is.null(comm_all)) comm_all %>% filter(target == hsc_name, prob > MIN_PROB_FOR_LR, pval <= PVALUE_CUTOFF) else NULL

plot_lr_bubble <- function(df, direction = c("outgoing", "incoming"), out_prefix) {
  direction <- match.arg(direction)
  if (is.null(df) || nrow(df) == 0) stop("No HSC-like LR communications for ", direction)
  partner_col <- if (direction == "outgoing") "target" else "source"
  df2 <- df %>%
    arrange(desc(prob)) %>%
    slice_head(n = TOP_N_LR) %>%
    mutate(
      partner = factor(.data[[partner_col]], levels = idents),
      LR_pair = fct_reorder(LR_pair, prob)
    )
  p <- ggplot(df2, aes(x = partner, y = LR_pair)) +
    geom_point(aes(size = prob, fill = pathway_name), shape = 21, color = "black", stroke = 0.20, alpha = 0.92) +
    scale_size_continuous(name = "Communication probability", range = c(2.0, 7.5)) +
    guides(size = solid_size_legend("Communication probability"), fill = "none", color = "none", colour = "none") +
    labs(x = ifelse(direction == "outgoing", "Target cell type", "Source cell type"),
         y = "Ligand-receptor pair") +
    theme_sc_pub(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.50),
      legend.key.size = unit(0.32, "cm"),
      legend.text = element_text(size = 8.5)
    )
  save_plot_pub(p, out_prefix, width = 9.0, height = max(5.5, 0.20 * nrow(df2) + 2.8))
}

safe_run("HSC05 HSC-like outgoing LR bubble", {
  plot_lr_bubble(comm_out, "outgoing", file.path(HSC_DIR, "HSC05_HSC_like_outgoing_LR_bubble"))
})

safe_run("HSC06 HSC-like incoming LR bubble", {
  plot_lr_bubble(comm_in, "incoming", file.path(HSC_DIR, "HSC06_HSC_like_incoming_LR_bubble"))
})

############################################################
## 10. HSC-like pathway strength barplots
############################################################
safe_run("HSC07 HSC-like outgoing pathway strength", {
  if (is.null(hsc_pathway_df) || nrow(hsc_pathway_df) == 0) stop("No HSC pathway score table")
  df <- hsc_pathway_df %>% filter(outgoing_strength > 0) %>% arrange(desc(outgoing_strength)) %>% slice_head(n = TOP_N_PATHWAYS) %>% mutate(pathway = fct_reorder(pathway, outgoing_strength))
  p <- ggplot(df, aes(x = outgoing_strength, y = pathway)) +
    geom_col(width = 0.72, fill = "#1F4E79", color = "black", linewidth = 0.16) +
    labs(x = "Outgoing pathway strength from HSC-like", y = NULL) +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55))
  save_plot_pub(p, file.path(HSC_DIR, "HSC07_HSC_like_outgoing_pathway_strength"), 6.5, max(4.8, 0.26 * nrow(df) + 1.6))
})

safe_run("HSC08 HSC-like incoming pathway strength", {
  if (is.null(hsc_pathway_df) || nrow(hsc_pathway_df) == 0) stop("No HSC pathway score table")
  df <- hsc_pathway_df %>% filter(incoming_strength > 0) %>% arrange(desc(incoming_strength)) %>% slice_head(n = TOP_N_PATHWAYS) %>% mutate(pathway = fct_reorder(pathway, incoming_strength))
  p <- ggplot(df, aes(x = incoming_strength, y = pathway)) +
    geom_col(width = 0.72, fill = "#8C6D31", color = "black", linewidth = 0.16) +
    labs(x = "Incoming pathway strength to HSC-like", y = NULL) +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55))
  save_plot_pub(p, file.path(HSC_DIR, "HSC08_HSC_like_incoming_pathway_strength"), 6.5, max(4.8, 0.26 * nrow(df) + 1.6))
})

############################################################
## 11. CellChat bubble plots for HSC-like pathways
############################################################
safe_run("HSC09 CellChat bubble HSC-like outgoing", {
  pathways_use <- head(pathways_hsc_out, TOP_N_PATHWAYS)
  if (length(pathways_use) == 0) stop("No HSC outgoing pathways")
  p <- netVisual_bubble(cellchat, sources.use = hsc_name, signaling = pathways_use, remove.isolate = FALSE) +
    theme_sc_pub(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.45))
  save_plot_pub(p, file.path(HSC_DIR, "HSC09_CellChat_bubble_HSC_like_outgoing_pathways"), 11, 7.5)
})

safe_run("HSC10 CellChat bubble HSC-like incoming", {
  pathways_use <- head(pathways_hsc_in, TOP_N_PATHWAYS)
  if (length(pathways_use) == 0) stop("No HSC incoming pathways")
  p <- netVisual_bubble(cellchat, targets.use = hsc_name, signaling = pathways_use, remove.isolate = FALSE) +
    theme_sc_pub(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.45))
  save_plot_pub(p, file.path(HSC_DIR, "HSC10_CellChat_bubble_HSC_like_incoming_pathways"), 11, 7.5)
})

############################################################
## 12. HSC11: valid ggplot heatmap of HSC-like pathway direction strength
############################################################
safe_run("HSC11 HSC-like pathway direction heatmap", {
  if (is.null(hsc_pathway_df) || nrow(hsc_pathway_df) == 0) stop("No HSC pathway score table")
  df <- hsc_pathway_df %>%
    filter(total_HSC_related_strength > 0) %>%
    arrange(desc(total_HSC_related_strength)) %>%
    slice_head(n = TOP_N_PATHWAYS) %>%
    select(pathway, outgoing_strength, incoming_strength) %>%
    pivot_longer(cols = c(outgoing_strength, incoming_strength), names_to = "direction", values_to = "strength") %>%
    mutate(
      direction = recode(direction, outgoing_strength = "Outgoing", incoming_strength = "Incoming"),
      pathway = factor(pathway, levels = rev(unique(pathway)))
    )
  p <- ggplot(df, aes(x = direction, y = pathway, fill = strength)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradientn(colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5"), name = "Strength") +
    labs(x = NULL, y = NULL) +
    theme_sc_pub(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.55),
          axis.text.x = element_text(face = "bold"))
  save_plot_pub(p, file.path(HSC_DIR, "HSC11_HSC_like_pathway_direction_strength_heatmap"), 5.6, max(5.0, 0.30 * length(unique(df$pathway)) + 1.8))
})

############################################################
## 12B. HSC11A/HSC11B: compact HSC-like outgoing/incoming signaling patterns
## Purpose:
##   1. Avoid duplicated default titles from netAnalysis_signalingRole_heatmap()
##   2. Reduce blank area
##   3. Keep complete pathway names and cell-type labels
##   4. Fit the lower-left/remaining space of Fig.2
############################################################

safe_run("HSC11A/HSC11B compact signaling role heatmaps", {
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for compact HSC11A/HSC11B heatmaps.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for compact HSC11A/HSC11B heatmaps.")
  }
  
  prob_arr <- cellchat@netP$prob
  if (is.null(prob_arr) || length(dim(prob_arr)) < 3) {
    stop("cellchat@netP$prob is missing or not a 3D array.")
  }
  
  dn <- dimnames(prob_arr)
  if (is.null(dn[[1]])) dn[[1]] <- levels(cellchat@idents)
  if (is.null(dn[[2]])) dn[[2]] <- levels(cellchat@idents)
  if (is.null(dn[[3]])) dn[[3]] <- cellchat@netP$pathways
  
  dn[[1]] <- shorten_celltype_name(dn[[1]])
  dn[[2]] <- shorten_celltype_name(dn[[2]])
  dimnames(prob_arr) <- dn
  
  cell_levels_use <- c(
    intersect(celltype_levels_short, intersect(dn[[1]], dn[[2]])),
    setdiff(intersect(dn[[1]], dn[[2]]), celltype_levels_short)
  )
  
  cell_levels_use <- intersect(cell_levels_use, idents)
  
  if (length(cell_levels_use) == 0) {
    cell_levels_use <- intersect(dn[[1]], dn[[2]])
  }
  
  cell_cols_use <- get_celltype_colors(cell_levels_use, use_short = TRUE)
  names(cell_cols_use) <- cell_levels_use
  
  build_signaling_role_mat <- function(pathways_use,
                                       pattern = c("outgoing", "incoming")) {
    pattern <- match.arg(pattern)
    
    pathways_use <- intersect(pathways_use, dimnames(prob_arr)[[3]])
    if (length(pathways_use) == 0) {
      stop("No valid pathways found for ", pattern)
    }
    
    mat_abs <- matrix(
      0,
      nrow = length(pathways_use),
      ncol = length(cell_levels_use),
      dimnames = list(pathways_use, cell_levels_use)
    )
    
    for (pw in pathways_use) {
      pmat <- prob_arr[, , pw, drop = TRUE]
      pmat <- pmat[
        intersect(cell_levels_use, rownames(pmat)),
        intersect(cell_levels_use, colnames(pmat)),
        drop = FALSE
      ]
      
      if (pattern == "outgoing") {
        v <- rowSums(pmat, na.rm = TRUE)
      } else {
        v <- colSums(pmat, na.rm = TRUE)
      }
      
      mat_abs[pw, names(v)] <- v
    }
    
    row_max <- apply(mat_abs, 1, max, na.rm = TRUE)
    row_max[!is.finite(row_max) | row_max <= 0] <- 1
    
    mat_rel <- mat_abs / row_max
    
    list(
      abs = mat_abs,
      rel = mat_rel
    )
  }
  
  save_compact_role_heatmap <- function(mat_list,
                                        filename_base,
                                        title_text,
                                        color_mode = c("Reds", "Blues"),
                                        width = 6.35,
                                        height = 4.45,
                                        dpi = 450) {
    color_mode <- match.arg(color_mode)
    
    mat_rel <- mat_list$rel
    mat_abs <- mat_list$abs
    
    if (nrow(mat_rel) == 0 || ncol(mat_rel) == 0) {
      stop("Empty heatmap matrix for ", title_text)
    }
    
    if (color_mode == "Reds") {
      col_fun <- circlize::colorRamp2(
        c(0, 0.5, 1),
        c("#FFFFFF", "#FCAE91", "#A50F15")
      )
    } else {
      col_fun <- circlize::colorRamp2(
        c(0, 0.5, 1),
        c("#FFFFFF", "#9ECAE1", "#084594")
      )
    }
    
    col_total <- colSums(mat_abs, na.rm = TRUE)
    row_total <- rowSums(mat_abs, na.rm = TRUE)
    
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      Strength = ComplexHeatmap::anno_barplot(
        col_total,
        height = grid::unit(7, "mm"),
        border = FALSE,
        gp = grid::gpar(fill = cell_cols_use[colnames(mat_rel)], col = NA),
        axis_param = list(
          gp = grid::gpar(fontsize = 6),
          side = "left"
        )
      ),
      Cell = ComplexHeatmap::anno_simple(
        colnames(mat_rel),
        col = cell_cols_use,
        height = grid::unit(2.2, "mm"),
        border = FALSE
      ),
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold")
    )
    
    right_anno <- ComplexHeatmap::rowAnnotation(
      Strength = ComplexHeatmap::anno_barplot(
        row_total,
        width = grid::unit(8, "mm"),
        border = TRUE,
        gp = grid::gpar(fill = "grey55", col = "black", lwd = 0.4),
        axis_param = list(
          gp = grid::gpar(fontsize = 6)
        )
      ),
      annotation_name_side = "bottom",
      annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold")
    )
    
    ht <- ComplexHeatmap::Heatmap(
      mat_rel,
      name = "Relative\nstrength",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      top_annotation = top_anno,
      right_annotation = right_anno,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 8.5, fontface = "plain"),
      column_names_gp = grid::gpar(fontsize = 7.0, fontface = "plain"),
      column_names_rot = 90,
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 7.5, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 7),
        legend_height = grid::unit(22, "mm")
      ),
      border = FALSE,
      rect_gp = grid::gpar(col = NA),
      column_title = title_text,
      column_title_gp = grid::gpar(fontsize = 10.5, fontface = "bold"),
      column_title_side = "top"
    )
    
    grDevices::pdf(
      paste0(filename_base, ".pdf"),
      width = width,
      height = height,
      useDingbats = FALSE
    )
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(1.5, 1.5, 1.5, 1.5), "mm")
    )
    grDevices::dev.off()
    
    grDevices::png(
      paste0(filename_base, ".png"),
      width = width * dpi,
      height = height * dpi,
      res = dpi
    )
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(1.5, 1.5, 1.5, 1.5), "mm")
    )
    grDevices::dev.off()
    
    invisible(filename_base)
  }
  
  if (is.null(hsc_pathway_df) || nrow(hsc_pathway_df) == 0) {
    stop("No HSC pathway score table found.")
  }
  
  pathways_out_compact <- hsc_pathway_df %>%
    dplyr::filter(outgoing_strength > 0) %>%
    dplyr::arrange(dplyr::desc(outgoing_strength)) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::pull(pathway)
  
  pathways_in_compact <- hsc_pathway_df %>%
    dplyr::filter(incoming_strength > 0) %>%
    dplyr::arrange(dplyr::desc(incoming_strength)) %>%
    dplyr::slice_head(n = 9) %>%
    dplyr::pull(pathway)
  
  mat_out <- build_signaling_role_mat(
    pathways_use = pathways_out_compact,
    pattern = "outgoing"
  )
  
  mat_in <- build_signaling_role_mat(
    pathways_use = pathways_in_compact,
    pattern = "incoming"
  )
  
  save_compact_role_heatmap(
    mat_list = mat_out,
    filename_base = file.path(
      HSC_DIR,
      "HSC11A_HSC_like_related_outgoing_signaling_role_heatmap"
    ),
    title_text = "HSC-like outgoing signaling patterns",
    color_mode = "Reds",
    width = 6.35,
    height = 4.35
  )
  
  save_compact_role_heatmap(
    mat_list = mat_in,
    filename_base = file.path(
      HSC_DIR,
      "HSC11B_HSC_like_related_incoming_signaling_role_heatmap"
    ),
    title_text = "HSC-like incoming signaling patterns",
    color_mode = "Blues",
    width = 6.35,
    height = 4.55
  )
})
############################################################
## 13. HSC12/HSC13 pathway-specific networks
############################################################
save_one_hsc_pathway_network <- function(pathway, direction = c("outgoing", "incoming"), idx = 1) {
  direction <- match.arg(direction)
  prefix <- ifelse(direction == "outgoing", "HSC12", "HSC13")
  subdir <- file.path(HSC_DIR, paste0(prefix, "_individual_pathway_networks"))
  safe_mkdir(subdir)
  out_base <- file.path(subdir, paste0(prefix, "_", sprintf("%02d", idx), "_", make_safe_filename(pathway), "_", direction))
  save_base_plot(out_base, 9.5, 8.5, {
    par(mfrow = c(1, 1), xpd = TRUE, mar = c(1.2, 1.2, 3.2, 1.2))
    if (direction == "outgoing") {
      netVisual_aggregate(
        cellchat,
        signaling = pathway,
        sources.use = hsc_name,
        layout = "circle",
        remove.isolate = TRUE,
        color.use = object_cols,
        vertex.weight = rep(UNIFORM_VERTEX_WEIGHT, length(object_levels)),
        vertex.weight.max = UNIFORM_VERTEX_WEIGHT,
        vertex.size.max = UNIFORM_VERTEX_SIZE_MAX,
        vertex.label.cex = PATHWAY_LABEL_CEX
      )
      title(main = paste0(pathway, " | HSC-like outgoing"), cex.main = 1.0, font.main = 2)
    } else {
      netVisual_aggregate(
        cellchat,
        signaling = pathway,
        targets.use = hsc_name,
        layout = "circle",
        remove.isolate = TRUE,
        color.use = object_cols,
        vertex.weight = rep(UNIFORM_VERTEX_WEIGHT, length(object_levels)),
        vertex.weight.max = UNIFORM_VERTEX_WEIGHT,
        vertex.size.max = UNIFORM_VERTEX_SIZE_MAX,
        vertex.label.cex = PATHWAY_LABEL_CEX
      )
      title(main = paste0(pathway, " | HSC-like incoming"), cex.main = 1.0, font.main = 2)
    }
  })
  out_base
}

safe_run("HSC12 top outgoing pathway networks", {
  pathways_use <- head(pathways_hsc_out, min(8, length(pathways_hsc_out)))
  if (length(pathways_use) == 0) stop("No HSC outgoing pathways")
  files <- character(0)
  for (i in seq_along(pathways_use)) {
    files <- c(files, tryCatch(save_one_hsc_pathway_network(pathways_use[i], "outgoing", i), error = function(e) NA_character_))
  }
  write.csv(data.frame(pathway = pathways_use, file_base = files),
            file.path(HSC_DIR, "HSC12_top_outgoing_pathway_networks_from_HSC_like_index.csv"), row.names = FALSE)
})

safe_run("HSC13 top incoming pathway networks", {
  pathways_use <- head(pathways_hsc_in, min(8, length(pathways_hsc_in)))
  if (length(pathways_use) == 0) stop("No HSC incoming pathways")
  files <- character(0)
  for (i in seq_along(pathways_use)) {
    files <- c(files, tryCatch(save_one_hsc_pathway_network(pathways_use[i], "incoming", i), error = function(e) NA_character_))
  }
  write.csv(data.frame(pathway = pathways_use, file_base = files),
            file.path(HSC_DIR, "HSC13_top_incoming_pathway_networks_to_HSC_like_index.csv"), row.names = FALSE)
})

############################################################
## 14. HSC14 communication pattern dot plots
## Restore original netAnalysis_dot style.
## Remove cell/pathway color legends; keep solid black Contribution legend.
############################################################
safe_run("HSC14 global outgoing and incoming communication pattern dot plots", {
  pathway_lr_labels <- build_pathway_lr_label_map(comm_all)

  relabel_pathways <- function(x) {
    out <- pathway_lr_labels[x]
    out[is.na(out)] <- x[is.na(out)]
    unname(out)
  }

  format_pattern_dot <- function(p) {
    p +
      scale_x_discrete(labels = relabel_pathways) +
      guides(
        size = solid_size_legend("Contribution"),
        color = "none",
        colour = "none",
        fill = "none"
      ) +
      labs(x = "Ligand-receptor pathway", y = NULL, color = NULL, fill = NULL) +
      theme_sc_pub(base_size = 10.5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.45),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.20),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 10, face = "bold")
      )
  }

  p_out <- format_pattern_dot(netAnalysis_dot(cellchat, pattern = "outgoing")) +
    ggtitle("Outgoing communication patterns of secreting cells")
  p_in <- format_pattern_dot(netAnalysis_dot(cellchat, pattern = "incoming")) +
    ggtitle("Incoming communication patterns of target cells")

  save_plot_pub_wide(
    p_out + p_in + plot_layout(widths = c(1, 1)),
    file.path(HSC_DIR, "HSC14_outgoing_incoming_communication_pattern_dotplot"),
    width = 16,
    height = 8.5
  )
})

message("CellChat publication figures exported to: ", OUT_DIR)
message("HSC-like specific figures exported to: ", HSC_DIR)

# END ORIGINAL SCRIPT P18

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P19: 12_SC_hdWGCNA.R
# ==============================================================================
# 初始化 ----
rm(list = ls())

# 加载包 ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(parallel)
library(scCustomize)
library(doParallel)
library(qs)
library(openxlsx)
library(future)

# 建议关闭多线程，减少结果波动
disableWGCNAThreads()

# 根目录(按需修改)
ANNO_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\2.Sc_annotation"
hdWGCNA_DIR <- "D:\\0.DataAnalysis\\AML_SC_RNA\\GSE116256\\1.16AML+5NORMAL\\6.hdWGCNA\\0.HSC-like_30hub_new"
if (!dir.exists(ANNO_DIR)) dir.create(ANNO_DIR, recursive = TRUE)
if (!dir.exists(hdWGCNA_DIR)) dir.create(hdWGCNA_DIR, recursive = TRUE)

# 图像风格辅助：随机色板
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}
randomColors <- replicate(100, randomColor())

# 统一图标题的辅助函数（自动编号）
.fig_id <- 0L
fig_id <- function() { .fig_id <<- .fig_id + 1L; .fig_id }
fig_title <- function(title) sprintf("Fig.%d %s", fig_id(), title)

# PDF 保存辅助
save_pdf <- function(file, expr, w = 8, h = 6) {
  pdf(file, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  print(expr)
}

# ---------------------------
# 1. 读入对象
# ---------------------------
setwd(ANNO_DIR)
seurat <- qread("2.resolution_0.5_annotation.qs")
seurat_obj <- seurat
setwd(hdWGCNA_DIR)

seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)

# 确保使用对象中原始、完整的 celltype
# 如果对象的 Idents 不是完整名称，就不要覆盖 celltype
# seurat_obj$celltype = Idents(seurat_obj)   # 这一句删掉，避免把长名覆盖成短名

# ---------------------------
# 2. 新增简写列：仅用于输出和绘图
# ---------------------------
celltype_map <- c(
  "B cell" = "Bcell",
  "cDC" = "cDC",
  "cMono" = "cMono",
  "Cycling Progenitor" = "CyclingProg",
  "Erythrocyte" = "Erythrocyte",
  "GMP-like cell" = "GMP-like",
  "HSC-like cell" = "HSC-like",
  "intMono" = "intMono",
  "LMPP-like cell" = "LymphoidProg",
  "Macrophage" = "Macrophage",
  "Mast cell" = "Mast",
  "MBP-like cell" = "GranuleProg",
  "MDP-like cell" = "MyeloidProg",
  "moDC" = "moDC",
  "Naive T cell" = "NaiveT",
  "ncMono" = "ncMono",
  "NK cell" = "NK_T",
  "pDC" = "pDC",
  "Plasma cell" = "Plasma",
  "Progenitor-like cell" = "Progenitor",
  "ProMyelocyte" = "Promyelo"
)

seurat_obj$celltype <- as.character(seurat_obj$celltype)
seurat_obj$celltype_short <- unname(celltype_map[seurat_obj$celltype])

# 检查是否有未映射成功的名称
unmapped_types <- unique(seurat_obj$celltype[is.na(seurat_obj$celltype_short)])
if (length(unmapped_types) > 0) {
  stop(
    paste0(
      "以下 celltype 没有在 celltype_map 中定义，请补充映射：\n",
      paste(unmapped_types, collapse = ", ")
    )
  )
}

# 为了后续绘图顺序稳定，把 short 设为因子，顺序按你的目标顺序
short_levels <- c(
  "Bcell", "cDC", "cMono", "CyclingProg", "Erythrocyte",
  "GMP-like", "GranuleProg", "HSC-like", "intMono", "LymphoidProg",
  "Macrophage", "Mast", "moDC", "MyeloidProg", "NaiveT",
  "ncMono", "NK_T", "pDC", "Plasma", "Progenitor", "Promyelo"
)
seurat_obj$celltype_short <- factor(seurat_obj$celltype_short, levels = short_levels)

# type 保持为 AML / Normal
seurat_obj$type <- as.character(seurat_obj$type)

# 输出检查
print(table(seurat_obj$type))
print(table(seurat_obj$celltype))
print(table(seurat_obj$celltype_short, useNA = "ifany"))
print(head(seurat_obj@meta.data[, c("orig.ident", "type", "celltype", "celltype_short")]))

# ---------------------------
# 3. hdWGCNA 预处理
# ---------------------------
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "hdWGCNA"
)

# 构造 metacells：分析阶段仍然用完整 celltype
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "type"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "type",
  min_cells = 100
)

# 对 Metacells 表达矩阵进行矫正
seurat_obj <- NormalizeMetacells(seurat_obj)

# ---------------------------
# 4. 设置分析目标细胞群：用完整名称
# ---------------------------
cellName <- "HSC-like cell"      # 分析阶段用完整名称
cellName_short <- "HSC-like"     # 输出阶段用简写名称

# 创建表达矩阵：仍然用完整 celltype
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cellName,
  group.by = "celltype",
  assay = "RNA",
  layer = "data"
)

# ---------------------------
# 5. 选择最佳 soft power
# ---------------------------
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed"
)

plot_list <- PlotSoftPowers(seurat_obj)
pdf(file = "Fig01_power.pdf", width = 8, height = 7)
wrap_plots(plot_list, ncol = 2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file = "Tab00_powerTable.csv", row.names = FALSE)

# !!! 按 Fig01_power.pdf 决定
selectPower <- 8

# ---------------------------
# 6. 构建网络
# ---------------------------
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = selectPower,
  tom_name = cellName_short   # 输出文件名用简写
)

pdf(file = "Fig02_cluster.pdf", width = 8, height = 6)
PlotDendrogram(seurat_obj, main = paste0(cellName_short, " hdWGCNA Dendrogram"))
dev.off()

# 获取 TOM
TOM <- GetTOM(seurat_obj)

# ---------------------------
# 7. 模块特征值与连接度
# ---------------------------
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# 这里仍然按完整 celltype 计算，保证分析一致
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "celltype"
)

hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized = FALSE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "celltype",
  group_name = cellName
)

# 模块名重命名：输出名用简写前缀
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(cellName_short, "-M")
)

# 重新取模块名
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized = FALSE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]

# ---------------------------
# 8. kME 图
# ---------------------------
pdf(file = "Fig03_kME.pdf", width = 15, height = 15)
PlotKMEs(seurat_obj, ncol = 5, n_hubs = 30)
dev.off()

# 模块基因表
modules <- GetModules(seurat_obj) %>% subset(module != "grey")

# 如果表里有 cell_type/group_name 等字段，可加简写列
if ("module" %in% colnames(modules)) {
  modules$target_celltype <- cellName_short
}

write.csv(modules, file = "Tab01_moduleGene.csv", row.names = FALSE)

# hub 基因表
hub_df <- GetHubGenes(seurat_obj, n_hubs = 30)
if (nrow(hub_df) > 0) {
  hub_df$target_celltype <- cellName_short
}
write.csv(hub_df, file = "Tab02_hubGene.csv", row.names = FALSE)

# ---------------------------
# 9. 模块评分
# ---------------------------
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 30,
  method = "Seurat"
)

# ---------------------------
# 10. FeaturePlot
# ---------------------------
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "tsne",
  features = "hMEs",
  order = TRUE
)

Ncol <- 5
Height <- ifelse(length(mods) / 5 * 3 > 6, length(mods) / 5 * 3, 6)

pdf(file = "Fig04_FeaturePlot.pdf", width = 15, height = Height)
wrap_plots(plot_list, ncol = Ncol)
dev.off()

# ---------------------------
# 11. 雷达图
# ---------------------------
# 这里 cluster 原代码是从 orig.ident 按空格拆分，通常没必要
# 直接用 type 或 orig.ident 更稳
seurat_obj$cluster <- seurat_obj$type

pdf(file = "Fig05_ModuleRadarPlot.pdf", width = 8, height = 8)
ModuleRadarPlot(
  seurat_obj,
  group.by = "cluster",
  barcodes = seurat_obj@meta.data %>%
    subset(celltype == cellName) %>%
    rownames(),
  axis.label.size = 4,
  grid.label.size = 4
)
dev.off()

# ---------------------------
# 12. 模块相关性图
# ---------------------------
pdf(file = "Fig06_ModuleCorrelogram.pdf", width = 8, height = 7)
ModuleCorrelogram(seurat_obj)
dev.off()

# ---------------------------
# 13. 模块-细胞类型相关性图与导表
# 输出阶段改用简写 celltype_short
# ---------------------------
MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p6 <- DotPlot(seurat_obj, features = mods, group.by = "celltype_short")
dotplot_data <- p6$data

# 导出宽表时也用简写名
dotplot_wide_avg <- dotplot_data %>%
  dplyr::select(features.plot, id, avg.exp) %>%
  tidyr::pivot_wider(names_from = id, values_from = avg.exp)

dotplot_wide_pct <- dotplot_data %>%
  dplyr::select(features.plot, id, pct.exp) %>%
  tidyr::pivot_wider(names_from = id, values_from = pct.exp)

write.xlsx(
  list(
    "AverageExpression" = dotplot_wide_avg,
    "PercentExpression" = dotplot_wide_pct
  ),
  file = "Tab03_cellModuleCor_wide.xlsx",
  rowNames = FALSE
)

p6 <- p6 +
  RotatedAxis() +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")

pdf(file = "Fig07_cellModuleCor.pdf", width = 9, height = 6)
print(p6)
dev.off()

# ---------------------------
# 14. 小提琴图
# 输出阶段改用简写 celltype_short
# ---------------------------
p7 <- VlnPlot(
  seurat_obj,
  features = mods,
  group.by = "celltype_short",
  ncol = Ncol,
  sort = "increasing",
  pt.size = 0
)
p7 <- p7 + geom_boxplot(width = 0.25, fill = "white")
p7 <- p7 + xlab("") + ylab("hME") + NoLegend()

pdf(file = "Fig08_VlnPlot.pdf", width = 16, height = Height)
print(p7)
dev.off()

# ---------------------------
# 15. 模块网络图
# ---------------------------
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "all")

# ---------------------------
# 16. hub 网络图
# ---------------------------
options(future.globals.maxSize = 20 * 1024^3)

pdf(file = "Fig09_hubNetwork.pdf", width = 7, height = 7)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3,
  n_other = 5,
  edge_prop = 0.75,
  mods = "all"
)
dev.off()

# ---------------------------
# 17. 基因 UMAP
# ---------------------------
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- GetModuleUMAP(seurat_obj)
write.csv(umap_df, file = "Tab04_geneUMAP.csv", row.names = FALSE)

pdf(file = "Fig10_geneUmap.pdf", width = 9, height = 7)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = FALSE
)
dev.off()

# ---------------------------
# 18. 保存对象
# ---------------------------
qsave(seurat_obj, file = "hdWGCNA_object.qs")






# END ORIGINAL SCRIPT P19

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P20: 13_single_gene_GSEA.R
# ==============================================================================
rm(list = ls())

############################################################
## Single-gene GSEA
## Target-gene correlation ranking -> KEGG GSEA
## Final version
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(openxlsx)
  library(stringr)
  library(cowplot)
  library(grid)
})

options(stringsAsFactors = FALSE)

############################################################
## 0) Global style
############################################################

COL_RED  <- "#CB5148"
COL_BLUE <- "#6A9BCB"

## single-pathway GSEA classic red-blue
COL_GSEA_HIGH <- "#D73027"
COL_GSEA_MID  <- "#F7F7F7"
COL_GSEA_LOW  <- "#2C7FB8"

COL_GREY_LINE <- "#D0D0D0"
COL_GRID <- "#EAEAEA"

theme_clean_gsea <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        color = "black"
      ),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

format_pvalue <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      "<0.001",
      sprintf("%.3g", x)
    )
  )
}

fmt_axis_label <- function(x) {
  x[abs(x) < 1e-10] <- 0
  out <- formatC(x, format = "f", digits = 2)
  out <- sub("0$", "", out)
  out <- sub("\\.0$", "", out)
  out
}

make_safe_range <- function(x, pad = 1e-8) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(0, 1))
  }
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) {
    r <- r + c(-pad, pad)
  }
  r
}

save_plot_highres <- function(plot_obj,
                              filename_prefix,
                              width = 8,
                              height = 6,
                              dpi = 300) {
  
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  ggsave(
    filename = paste0(filename_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = pdf_device
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "tiff",
    compression = "lzw"
  )
  
  ggsave(
    filename = paste0(filename_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "png",
    bg = "white"
  )
}

############################################################
## 1) User parameters
############################################################

setwd("D:/0.DataAnalysis/AML_BULK_RNA/小论文/GSE30029_new/6.single_gene_GSEA")

expr_file   <- "raw_data.xlsx"
target_gene <- "GNL3L"

## "All" or "AML_only"
analysis_scope <- "All"

n_healthy <- 31
n_aml     <- 46

## 手动指定 4 条单通路图展示内容
## 支持 KEGG ID 或通路名称关键词
single_pathway_terms <- c(
  "Cell cycle",
  "DNA replication",
  "Oxidative phosphorylation",
  "TNF signaling pathway"
)

## 如果想自动选 4 条，改成：
# single_pathway_terms <- character(0)

target_gene <- toupper(target_gene)

############################################################
## 2) Read expression matrix
############################################################

expr_df <- read_xlsx(expr_file)
colnames(expr_df)[1] <- "gene_symbol"

expr_df$gene_symbol <- as.character(expr_df$gene_symbol)
expr_df$gene_symbol <- trimws(expr_df$gene_symbol)
expr_df$gene_symbol <- toupper(expr_df$gene_symbol)

expr_num <- expr_df[, -1, drop = FALSE]
expr_num[] <- lapply(expr_num, function(x) as.numeric(as.character(x)))

expr_df2 <- data.frame(
  gene_symbol = expr_df$gene_symbol,
  expr_num,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_df2 <- expr_df2 %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

expr_mat <- as.matrix(expr_df2[, -1, drop = FALSE])
rownames(expr_mat) <- expr_df2$gene_symbol
storage.mode(expr_mat) <- "numeric"

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

if (!(target_gene %in% rownames(expr_mat))) {
  stop(paste0("Target gene ", target_gene, " not found in expression matrix."))
}

############################################################
## 3) Sample metadata
############################################################

expected_n <- n_healthy + n_aml

if (ncol(expr_mat) != expected_n) {
  stop(
    paste0(
      "Sample number mismatch: expression matrix has ",
      ncol(expr_mat),
      " samples, but n_healthy + n_aml = ",
      expected_n,
      ". Please check n_healthy and n_aml."
    )
  )
}

sample_info <- data.frame(
  sample_id = colnames(expr_mat),
  group = factor(
    c(rep("Healthy", n_healthy), rep("AML", n_aml)),
    levels = c("Healthy", "AML")
  ),
  stringsAsFactors = FALSE
)

if (analysis_scope == "AML_only") {
  keep_samples <- sample_info$group == "AML"
  expr_use <- expr_mat[, keep_samples, drop = FALSE]
  sample_info_use <- sample_info[keep_samples, , drop = FALSE]
} else if (analysis_scope == "All") {
  expr_use <- expr_mat
  sample_info_use <- sample_info
} else {
  stop("analysis_scope must be either 'All' or 'AML_only'.")
}

cat("Analysis scope:", analysis_scope, "\n")
cat("Samples used:", ncol(expr_use), "\n")

############################################################
## 4) Compute target-gene correlation with all genes
############################################################

target_expr <- as.numeric(expr_use[target_gene, ])
names(target_expr) <- colnames(expr_use)

if (sd(target_expr, na.rm = TRUE) == 0) {
  stop(paste0(target_gene, " has zero variance in selected samples."))
}

cor_vec <- apply(expr_use, 1, function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  cor(
    x,
    target_expr,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
})

cor_vec <- cor_vec[!is.na(cor_vec)]

gene_cor_df <- data.frame(
  SYMBOL = names(cor_vec),
  spearman_r = as.numeric(cor_vec),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(spearman_r))

write.xlsx(
  gene_cor_df,
  file = paste0("01_", target_gene, "_correlation_ranking_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

############################################################
## 5) SYMBOL -> ENTREZID
############################################################

id_map <- bitr(
  geneID = gene_cor_df$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_df <- gene_cor_df %>%
  left_join(id_map, by = "SYMBOL") %>%
  filter(!is.na(ENTREZID))

gene_df <- gene_df %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs(spearman_r),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

gene_list <- gene_df$spearman_r
names(gene_list) <- gene_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Genes used for GSEA:", length(gene_list), "\n")

if (length(gene_list) < 100) {
  warning("The ranked gene list contains fewer than 100 genes. GSEA results may be unstable.")
}

############################################################
## 6) KEGG GSEA
############################################################

gsea_result <- gseKEGG(
  geneList      = gene_list,
  organism      = "hsa",
  keyType       = "kegg",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  eps           = 0
)

gsea_df <- as.data.frame(gsea_result)

write.xlsx(
  gsea_df,
  file = paste0("02_", target_gene, "_single_gene_GSEA_KEGG_all_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

sig_df <- gsea_df %>%
  filter(p.adjust < 0.05)

write.xlsx(
  sig_df,
  file = paste0("03_", target_gene, "_single_gene_GSEA_KEGG_sig_FDR0.05_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

gsea_pos <- sig_df %>%
  filter(NES > 0)

gsea_neg <- sig_df %>%
  filter(NES < 0)

write.xlsx(
  gsea_pos,
  file = paste0("04_", target_gene, "_single_gene_GSEA_KEGG_positive_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

write.xlsx(
  gsea_neg,
  file = paste0("05_", target_gene, "_single_gene_GSEA_KEGG_negative_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

if (nrow(sig_df) == 0) {
  stop("No significant KEGG pathways found at FDR < 0.05.")
}

############################################################
## 7) Ridgeplot
############################################################

p_ridge <- ridgeplot(
  gsea_result,
  showCategory = min(20, nrow(gsea_df)),
  core_enrichment = TRUE,
  fill = "p.adjust"
) +
  scale_fill_gradient(
    low = COL_RED,
    high = COL_BLUE,
    name = "p.adjust"
  ) +
  labs(
    title = paste0(target_gene, " single-gene GSEA: KEGG ridgeplot"),
    x = "Running enrichment distribution"
  ) +
  theme_clean_gsea(12)

save_plot_highres(
  p_ridge,
  paste0("06_", target_gene, "_GSEA_KEGG_ridgeplot_", analysis_scope),
  width = 9.5,
  height = 12.5,
  dpi = 300
)

############################################################
## 8) Overall KEGG lollipop plot
############################################################

n_show_each <- 10

## 固定 NES 坐标范围
X_LIMIT <- 3.5
X_BREAKS <- c(-3.5, -1.75, 0, 1.75, 3.5)

overall_pos <- sig_df %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = n_show_each) %>%
  mutate(
    KEGG_pathway_label = str_wrap(Description, width = 42)
  )

overall_neg <- sig_df %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = n_show_each) %>%
  mutate(
    KEGG_pathway_label = str_wrap(Description, width = 42)
  )

overall_plot_df <- bind_rows(overall_pos, overall_neg) %>%
  distinct(ID, .keep_all = TRUE)

if (nrow(overall_plot_df) == 0) {
  stop("No positive or negative significant pathways for overall lollipop plot.")
}

if (max(abs(overall_plot_df$NES), na.rm = TRUE) > X_LIMIT) {
  warning(
    "Some NES values are outside [-3.5, 3.5] and will be clipped in the plot."
  )
}

## 上半部分 positive enrichment，下半部分 negative enrichment
pathway_order_top_to_bottom <- c(
  overall_pos %>%
    arrange(desc(NES)) %>%
    pull(KEGG_pathway_label),
  overall_neg %>%
    arrange(NES) %>%
    pull(KEGG_pathway_label)
)

pathway_order_top_to_bottom <- unique(pathway_order_top_to_bottom)

n_pathway <- length(pathway_order_top_to_bottom)

y_label_df <- tibble(
  KEGG_pathway_label = pathway_order_top_to_bottom,
  y_num = n_pathway - seq_along(pathway_order_top_to_bottom) + 1
)

overall_plot_df <- overall_plot_df %>%
  left_join(y_label_df, by = "KEGG_pathway_label") %>%
  mutate(
    NES_plot = pmax(pmin(NES, X_LIMIT), -X_LIMIT)
  )

count_range <- make_safe_range(overall_plot_df$setSize)
pval_range  <- make_safe_range(overall_plot_df$p.adjust)

STRIP_BOTTOM <- n_pathway + 0.52
STRIP_TOP    <- n_pathway + 1.32
STRIP_MID    <- (STRIP_BOTTOM + STRIP_TOP) / 2

strip_df <- tibble(
  xmin = c(-X_LIMIT, 0),
  xmax = c(0, X_LIMIT),
  ymin = STRIP_BOTTOM,
  ymax = STRIP_TOP,
  label = c(
    paste0("Negative enrichment\n(", target_gene, "-high)"),
    paste0("Positive enrichment\n(", target_gene, "-high)")
  ),
  xmid = c(-X_LIMIT / 2, X_LIMIT / 2)
)

p_overall <- ggplot() +
  geom_rect(
    data = strip_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "#E5E5E5",
    color = "black",
    linewidth = 0.65
  ) +
  geom_text(
    data = strip_df,
    aes(
      x = xmid,
      y = STRIP_MID,
      label = label
    ),
    fontface = "bold",
    size = 3.65,
    color = "black",
    lineheight = 0.95
  ) +
  geom_segment(
    data = overall_plot_df,
    aes(
      x = 0,
      xend = NES_plot,
      y = y_num,
      yend = y_num
    ),
    color = COL_GREY_LINE,
    linewidth = 1.10,
    lineend = "round"
  ) +
  geom_point(
    data = overall_plot_df,
    aes(
      x = NES_plot,
      y = y_num,
      size = setSize,
      color = p.adjust
    ),
    shape = 16
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.55,
    color = "grey35"
  ) +
  scale_x_continuous(
    limits = c(-X_LIMIT, X_LIMIT),
    breaks = X_BREAKS,
    labels = fmt_axis_label,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_label_df$y_num,
    labels = y_label_df$KEGG_pathway_label,
    limits = c(0.5, STRIP_TOP),
    expand = c(0, 0)
  ) +
  scale_color_gradient(
    low = COL_RED,
    high = COL_BLUE,
    limits = pval_range,
    name = "p.adjust"
  ) +
  scale_size_area(
    max_size = 6.8,
    limits = count_range,
    name = "Count"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        colour = "black",
        fill = "black"
      )
    ),
    color = guide_colorbar(
      order = 2,
      barheight = unit(42, "pt"),
      barwidth = unit(8, "pt")
    )
  ) +
  labs(
    x = "NES",
    y = NULL
  ) +
  theme_clean_gsea(12.5) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.65
    ),
    panel.grid.major = element_line(
      color = "#EAEAEA",
      linewidth = 0.38
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 11.5,
      face = "bold",
      color = "black",
      lineheight = 0.96,
      margin = margin(r = 4)
    ),
    axis.text.x = element_text(
      size = 10.8,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 11.5,
      color = "black",
      margin = margin(t = 4)
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.40
    ),
    legend.title = element_text(
      size = 10.5,
      color = "black"
    ),
    legend.text = element_text(
      size = 9.8,
      color = "black"
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(4, "pt"),
    plot.margin = margin(5, 6, 6, 6)
  )

overall_height <- max(9.6, 0.54 * nrow(overall_plot_df) + 2.4)

save_plot_highres(
  p_overall,
  paste0("07_", target_gene, "_GSEA_KEGG_overall_lollipop_", analysis_scope),
  width = 7.2,
  height = overall_height,
  dpi = 300
)

write.xlsx(
  overall_plot_df,
  paste0("07_", target_gene, "_GSEA_KEGG_overall_lollipop_source_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

############################################################
## 9) Single-pathway GSEA template plots
############################################################

calculate_running_es <- function(gene_list,
                                 gene_set,
                                 exponent = 1) {
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_names <- names(gene_list)
  
  gene_set <- intersect(as.character(gene_set), gene_names)
  
  N <- length(gene_list)
  hit_ind <- gene_names %in% gene_set
  Nh <- sum(hit_ind)
  
  if (Nh == 0) {
    stop("No genes in this gene set were found in the ranked gene list.")
  }
  
  if (Nh == N) {
    stop("All ranked genes are included in this gene set; running ES cannot be calculated.")
  }
  
  weights <- abs(gene_list)^exponent
  
  hit_score <- ifelse(
    hit_ind,
    weights / sum(weights[hit_ind]),
    0
  )
  
  miss_score <- ifelse(
    !hit_ind,
    1 / (N - Nh),
    0
  )
  
  running_ES <- cumsum(hit_score - miss_score)
  
  data.frame(
    rank = seq_len(N),
    gene = gene_names,
    rank_metric = as.numeric(gene_list),
    hit = hit_ind,
    running_ES = running_ES
  )
}

plot_single_gsea_template <- function(gsea_obj,
                                      gene_list,
                                      geneSetID,
                                      pvalX = 0.73,
                                      pvalY = 0.83,
                                      title_width = 42) {
  
  result_df <- as.data.frame(gsea_obj)
  
  one_result <- result_df %>%
    filter(ID == geneSetID)
  
  if (nrow(one_result) == 0) {
    stop(paste0("geneSetID not found in GSEA result: ", geneSetID))
  }
  
  gene_set <- gsea_obj@geneSets[[geneSetID]]
  
  es_df <- calculate_running_es(
    gene_list = gene_list,
    gene_set = gene_set,
    exponent = 1
  )
  
  hit_df <- es_df %>%
    filter(hit)
  
  N <- nrow(es_df)
  
  y_min <- min(es_df$running_ES, na.rm = TRUE)
  y_max <- max(es_df$running_ES, na.rm = TRUE)
  y_span <- y_max - y_min
  
  x_text <- N * pvalX
  y_text <- y_min + y_span * pvalY
  
  pathway_title <- one_result$Description[1]
  pathway_title <- str_wrap(pathway_title, width = title_width)
  
  label_text <- paste0(
    "NES: ", sprintf("%.2f", one_result$NES[1]), "\n",
    "Pvalue: ", format_pvalue(one_result$pvalue[1])
  )
  
  p_curve <- ggplot(
    es_df,
    aes(x = rank, y = running_ES)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.50,
      color = "grey35"
    ) +
    geom_segment(
      data = hit_df,
      aes(
        x = rank,
        xend = rank,
        y = 0,
        yend = running_ES,
        color = rank_metric
      ),
      linewidth = 0.30,
      alpha = 0.78
    ) +
    geom_line(
      aes(color = rank_metric),
      linewidth = 1.15,
      alpha = 0.95
    ) +
    geom_point(
      aes(color = rank_metric),
      size = 0.28,
      alpha = 0.35
    ) +
    annotate(
      "text",
      x = x_text,
      y = y_text,
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 3.5,
      color = "black",
      fontface = "italic"
    ) +
    scale_color_gradient2(
      low = COL_GSEA_LOW,
      mid = COL_GSEA_MID,
      high = COL_GSEA_HIGH,
      midpoint = 0,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(1, N),
      expand = c(0, 0)
    ) +
    labs(
      title = pathway_title,
      x = NULL,
      y = "Running Enrichment Score"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 14.5,
        lineheight = 0.95
      ),
      axis.title.y = element_text(
        color = "black",
        size = 11.5
      ),
      axis.text.y = element_text(
        color = "black",
        size = 9.5
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_line(
        color = "#EEEEEE",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.60
      ),
      plot.margin = margin(4, 4, 0, 4)
    )
  
  n_bins <- 12
  
  bin_df <- es_df %>%
    mutate(
      bin = cut(
        rank,
        breaks = n_bins,
        labels = FALSE
      )
    ) %>%
    group_by(bin) %>%
    summarise(
      xmin = min(rank),
      xmax = max(rank),
      rank_metric = mean(rank_metric, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_heat <- ggplot(bin_df) +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = 1,
        fill = rank_metric
      ),
      color = NA
    ) +
    scale_fill_gradient2(
      low = COL_GSEA_LOW,
      mid = COL_GSEA_MID,
      high = COL_GSEA_HIGH,
      midpoint = 0,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(1, N),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    theme_void() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.60
      ),
      plot.margin = margin(2, 4, 4, 4)
    )
  
  cowplot::plot_grid(
    p_curve,
    p_heat,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(1, 0.12)
  )
}

############################################################
## 9.1) Select single pathways for 4-panel plot
############################################################

select_auto_single_terms <- function(df, n_total = 4) {
  
  n_pos <- ceiling(n_total / 2)
  n_neg <- floor(n_total / 2)
  
  pos_ids <- df %>%
    filter(NES > 0) %>%
    arrange(p.adjust, pvalue, desc(NES)) %>%
    slice_head(n = n_pos) %>%
    pull(ID)
  
  neg_ids <- df %>%
    filter(NES < 0) %>%
    arrange(p.adjust, pvalue, NES) %>%
    slice_head(n = n_neg) %>%
    pull(ID)
  
  selected_ids <- c(pos_ids, neg_ids)
  
  if (length(selected_ids) < n_total) {
    extra_ids <- df %>%
      filter(!ID %in% selected_ids) %>%
      arrange(p.adjust, desc(abs(NES))) %>%
      slice_head(n = n_total - length(selected_ids)) %>%
      pull(ID)
    
    selected_ids <- c(selected_ids, extra_ids)
  }
  
  unique(selected_ids)
}

resolve_single_terms <- function(terms, result_df, sig_result_df) {
  
  if (length(terms) == 0) {
    return(select_auto_single_terms(sig_result_df, n_total = 4))
  }
  
  resolved_ids <- c()
  
  for (one_term in terms) {
    
    one_term_trim <- str_trim(one_term)
    
    if (one_term_trim %in% result_df$ID) {
      resolved_ids <- c(resolved_ids, one_term_trim)
      next
    }
    
    exact_desc <- result_df %>%
      filter(tolower(Description) == tolower(one_term_trim)) %>%
      pull(ID)
    
    if (length(exact_desc) > 0) {
      resolved_ids <- c(resolved_ids, exact_desc[1])
      next
    }
    
    fuzzy_desc <- result_df %>%
      filter(str_detect(
        tolower(Description),
        fixed(tolower(one_term_trim))
      )) %>%
      arrange(p.adjust, pvalue) %>%
      pull(ID)
    
    if (length(fuzzy_desc) > 0) {
      resolved_ids <- c(resolved_ids, fuzzy_desc[1])
      next
    }
    
    warning(paste0("No matched pathway was found for: ", one_term_trim))
  }
  
  unique(resolved_ids)
}

single_pathway_ids <- resolve_single_terms(
  terms = single_pathway_terms,
  result_df = gsea_df,
  sig_result_df = sig_df
)

if (length(single_pathway_ids) == 0) {
  stop("No single pathway was selected. Please check single_pathway_terms.")
}

single_pathway_info <- gsea_df %>%
  filter(ID %in% single_pathway_ids) %>%
  arrange(match(ID, single_pathway_ids))

print(
  single_pathway_info[, c("ID", "Description", "NES", "pvalue", "p.adjust")]
)

write.xlsx(
  single_pathway_info,
  paste0("08_", target_gene, "_GSEA_single_pathway_selected_terms_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

############################################################
## 9.2) Generate 4-panel single-pathway plots
############################################################

gsea_single_list <- lapply(
  single_pathway_ids,
  function(x) {
    plot_single_gsea_template(
      gsea_obj = gsea_result,
      gene_list = gene_list,
      geneSetID = x,
      pvalX = 0.73,
      pvalY = 0.83,
      title_width = 42
    )
  }
)

p_single_combined <- cowplot::plot_grid(
  plotlist = gsea_single_list,
  ncol = 2,
  align = "hv"
)

save_plot_highres(
  p_single_combined,
  paste0("09_", target_gene, "_GSEA_single_pathway_template_style_combined_", analysis_scope),
  width = 10.5,
  height = 8.4,
  dpi = 300
)

############################################################
## 9.3) Save each single pathway separately
############################################################

for (i in seq_along(gsea_single_list)) {
  
  one_id <- single_pathway_ids[i]
  
  one_name <- gsea_df$Description[
    match(one_id, gsea_df$ID)
  ]
  
  one_name_clean <- one_name %>%
    str_replace_all("[/\\\\:*?\"<>|]", "_") %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_sub(1, 80)
  
  save_plot_highres(
    gsea_single_list[[i]],
    paste0(
      "10_",
      target_gene,
      "_GSEA_single_pathway_",
      i,
      "_",
      one_id,
      "_",
      one_name_clean,
      "_",
      analysis_scope
    ),
    width = 5.5,
    height = 4.9,
    dpi = 300
  )
}

############################################################
## 10) Top pathways summary lollipop plot
############################################################

sig_df2 <- sig_df %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      paste0(target_gene, "-positive"),
      paste0(target_gene, "-negative")
    ),
    Direction = factor(
      Direction,
      levels = c(
        paste0(target_gene, "-negative"),
        paste0(target_gene, "-positive")
      )
    ),
    absNES = abs(NES),
    KEGG_pathway_label = str_wrap(Description, width = 40)
  )

top_n <- min(25, nrow(sig_df2))

top25_paths <- sig_df2 %>%
  arrange(desc(absNES)) %>%
  slice_head(n = top_n) %>%
  mutate(
    KEGG_pathway_label = factor(
      KEGG_pathway_label,
      levels = rev(unique(KEGG_pathway_label))
    )
  )

p_summary <- ggplot(
  top25_paths,
  aes(x = NES, y = KEGG_pathway_label)
) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = KEGG_pathway_label,
      yend = KEGG_pathway_label
    ),
    color = COL_GREY_LINE,
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    aes(
      size = setSize,
      color = p.adjust
    ),
    shape = 16
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    color = "grey40"
  ) +
  scale_color_gradient(
    low = COL_RED,
    high = COL_BLUE,
    name = "p.adjust"
  ) +
  scale_size_area(
    max_size = 6.8,
    name = "Count"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        colour = "black",
        fill = "black"
      )
    ),
    color = guide_colorbar(order = 2)
  ) +
  labs(
    title = paste0("Top ", top_n, " KEGG pathways ranked by |NES|"),
    subtitle = paste0(target_gene, " single-gene GSEA (", analysis_scope, ")"),
    x = "NES",
    y = NULL
  ) +
  theme_clean_gsea(12) +
  theme(
    axis.text.y = element_text(size = 9.8, color = "black"),
    axis.text.x = element_text(size = 10.5, color = "black"),
    legend.position = "right"
  )

save_plot_highres(
  p_summary,
  paste0("11_", target_gene, "_GSEA_summary_top25_absNES_", analysis_scope),
  width = 9.5,
  height = 10,
  dpi = 300
)

write.xlsx(
  top25_paths,
  file = paste0("12_", target_gene, "_GSEA_top25_absNES_", analysis_scope, ".xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

############################################################
## 11) Bar-style summary plot
############################################################

direction_cols <- setNames(
  object = c(COL_BLUE, COL_RED),
  nm = c(
    paste0(target_gene, "-negative"),
    paste0(target_gene, "-positive")
  )
)

p_bar <- ggplot(
  top25_paths,
  aes(
    x = reorder(KEGG_pathway_label, NES),
    y = NES,
    fill = Direction
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = direction_cols,
    name = "Direction"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.6,
    color = "black"
  ) +
  labs(
    title = paste0("Top ", top_n, " KEGG pathways ranked by |NES|"),
    subtitle = paste0(target_gene, " single-gene GSEA (", analysis_scope, ")"),
    x = NULL,
    y = "NES"
  ) +
  theme_clean_gsea(12) +
  theme(
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(size = 10.5, color = "black")
  )

save_plot_highres(
  p_bar,
  paste0("13_", target_gene, "_GSEA_barplot_top25_absNES_", analysis_scope),
  width = 9.5,
  height = 10,
  dpi = 300
)

############################################################
## 12) Output information
############################################################

cat("DONE: revised single-gene GSEA finished for ", target_gene, ".\n", sep = "")
cat("Main outputs:\n")
cat("01_", target_gene, "_correlation_ranking_", analysis_scope, ".xlsx\n", sep = "")
cat("02_", target_gene, "_single_gene_GSEA_KEGG_all_", analysis_scope, ".xlsx\n", sep = "")
cat("03_", target_gene, "_single_gene_GSEA_KEGG_sig_FDR0.05_", analysis_scope, ".xlsx\n", sep = "")
cat("06_", target_gene, "_GSEA_KEGG_ridgeplot_", analysis_scope, ".pdf/.tiff/.png\n", sep = "")
cat("07_", target_gene, "_GSEA_KEGG_overall_lollipop_", analysis_scope, ".pdf/.tiff/.png\n", sep = "")
cat("09_", target_gene, "_GSEA_single_pathway_template_style_combined_", analysis_scope, ".pdf/.tiff/.png\n", sep = "")
cat("11_", target_gene, "_GSEA_summary_top25_absNES_", analysis_scope, ".pdf/.tiff/.png\n", sep = "")
cat("13_", target_gene, "_GSEA_barplot_top25_absNES_", analysis_scope, ".pdf/.tiff/.png\n", sep = "")


# END ORIGINAL SCRIPT P20

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P21: 13_BULK_SC_gene_intersection_AND_single_gene_GSEA/BULK_SC_gene_upset_inputs/upset_evidence_plots.R
# ==============================================================================
############################################################
## Fig3 UpSet / Venn / Evidence matrix plotting FINAL_FIXED
## Project: AML小论文_GNL3L_TAM
##
## Input root:
##   D:/Users/lab01/Desktop/AML小论文_GNL3L_TAM/0.Fig3_upset_inputs
##
## Output:
##   Fig3_upset_outputs/
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################################
## 0. Packages
############################################################

need_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "stringr", "readr",
  "pheatmap", "RColorBrewer", "grid", "VennDiagram",
  "patchwork"
)

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[INFO] Installing missing package: ", pkg)
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(need_pkgs)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(VennDiagram)
  library(patchwork)
})

############################################################
## 1. Paths
############################################################

BASE_DIR <- "D:/Users/lab01/Desktop/AML小论文_GNL3L_TAM/0.Fig3_upset_inputs"

OUT_DIR <- file.path(BASE_DIR, "Fig3_upset_outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_DIR <- file.path(BASE_DIR, "01_gene_integration_upset")
TF_DIR   <- file.path(BASE_DIR, "02_TF_venn_or_upset")
DRUG_DIR <- file.path(BASE_DIR, "03_drug_target_evidence")

############################################################
## 2. Publication-style colors
############################################################

fig3_source_cols <- c(
  "Bulk_DEG"      = "#D55E5E",
  "Bulk_WGCNA"    = "#3F5FA8",
  "Sc_marker"     = "#2E8B57",
  "Sc_hdWGCNA"    = "#9A9A00",
  
  "CistromeDB_TF"     = "#2E8B57",
  "Motif_TF"          = "#D9C400",
  "HSC_likeMarkerTF"  = "#3F5FA8",
  
  "TF_DGIdb"      = "#4DBBD5",
  "TF_PanDrugs"   = "#B77DCC",
  "BULK_EpiMed"   = "#86C166",
  "BULK_L2S2"     = "#F28E8C"
)

fig3_main_cols <- c(
  "gene_upset" = "#2F4F4F",
  "tf_upset"   = "#4A4A4A",
  "drug_upset" = "#3A3A3A"
)

drug_source_cols <- c(
  "TF_DGIdb"    = "#4DBBD5",
  "TF_PanDrugs" = "#B77DCC",
  "BULK_EpiMed" = "#86C166",
  "BULK_L2S2"   = "#F28E8C"
)

evidence_count_cols <- colorRampPalette(
  c("#EAF4F4", "#9BD7C5", "#2E8B57")
)(100)

############################################################
## 3. Helper functions
############################################################

find_file_recursive <- function(root_dir, filename) {
  hits <- list.files(
    root_dir,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hits) == 0) {
    stop("Cannot find file: ", filename, " under ", root_dir)
  }
  
  hits[1]
}

read_one_column_txt <- function(file,
                                to_upper = TRUE,
                                remove_header_words = TRUE) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  
  if (length(x) == 0) {
    x <- readLines(file, warn = FALSE)
  }
  
  x <- gsub("\ufeff", "", x)
  x <- trimws(x)
  x <- x[x != ""]
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'|'$", "", x)
  x <- trimws(x)
  
  x <- sapply(strsplit(x, "\t|,"), `[`, 1)
  x <- trimws(x)
  
  if (remove_header_words) {
    bad <- tolower(x) %in% c(
      "gene", "genes", "symbol", "genesymbol", "gene_symbol",
      "tf", "transcription_factor",
      "drug", "drugs",
      "candidate", "name", "item", "x"
    )
    x <- x[!bad]
  }
  
  x <- x[!is.na(x) & x != ""]
  
  if (to_upper) {
    x <- toupper(x)
  }
  
  unique(x)
}

make_membership_df <- function(set_list) {
  set_list <- lapply(set_list, function(x) {
    unique(x[!is.na(x) & x != ""])
  })
  
  all_items <- sort(unique(unlist(set_list)))
  
  df <- data.frame(
    item = all_items,
    stringsAsFactors = FALSE
  )
  
  for (nm in names(set_list)) {
    df[[nm]] <- as.integer(df$item %in% set_list[[nm]])
  }
  
  df$Evidence_count <- rowSums(df[names(set_list)], na.rm = TRUE)
  
  df
}

write_membership_outputs <- function(df, set_names, prefix) {
  write.csv(
    df,
    paste0(prefix, "_membership_matrix.csv"),
    row.names = FALSE
  )
  
  common_all <- df %>%
    dplyr::filter(Evidence_count == length(set_names)) %>%
    dplyr::arrange(item)
  
  write.csv(
    common_all,
    paste0(prefix, "_all_sets_intersection.csv"),
    row.names = FALSE
  )
  
  multi_evidence <- df %>%
    dplyr::filter(Evidence_count >= 2) %>%
    dplyr::arrange(dplyr::desc(Evidence_count), item)
  
  write.csv(
    multi_evidence,
    paste0(prefix, "_multi_evidence_candidates.csv"),
    row.names = FALSE
  )
  
  invisible(
    list(
      common_all = common_all,
      multi_evidence = multi_evidence
    )
  )
}

############################################################
## 4. Custom publication-style UpSet plot
## Revised:
##   1. remove y-axis title in the upper-right bar plot
##   2. move the left panels rightward and make the layout compact
##   3. draw inactive grey dots first, then solid connection lines,
##      then active black dots, so the vertical lines are continuous
############################################################

save_upset_plot <- function(df,
                            set_names,
                            filename_base,
                            mainbar_y = "Intersecting genes",
                            sets_x = "Genes per source",
                            width = 7.2,
                            height = 4.8,
                            nintersects = 15,
                            main_bar_color = "#2F4F4F",
                            set_bar_colors = NULL,
                            base_size = 11) {
  
  message("[INFO] Drawing custom UpSet plot: ", basename(filename_base))
  
  missing_sets <- setdiff(set_names, colnames(df))
  if (length(missing_sets) > 0) {
    stop("Missing set columns in df: ", paste(missing_sets, collapse = ", "))
  }
  
  if (is.null(set_bar_colors)) {
    set_bar_colors <- rep("#4D4D4D", length(set_names))
    names(set_bar_colors) <- set_names
  } else {
    tmp_cols <- rep("#4D4D4D", length(set_names))
    names(tmp_cols) <- set_names
    tmp_cols[names(set_bar_colors)] <- set_bar_colors[names(set_bar_colors)]
    set_bar_colors <- tmp_cols
  }
  
  ##########################################################
  ## 1. Prepare membership matrix
  ##########################################################
  
  plot_df <- df[, set_names, drop = FALSE]
  plot_df[] <- lapply(plot_df, function(x) as.integer(as.numeric(x) > 0))
  
  if (nrow(plot_df) == 0) {
    stop("No rows available for UpSet plot: ", filename_base)
  }
  
  set_sizes <- colSums(plot_df[, set_names, drop = FALSE], na.rm = TRUE)
  if (all(set_sizes == 0)) {
    stop("All sets are empty for UpSet plot: ", filename_base)
  }
  
  ##########################################################
  ## 2. Build intersection table
  ##########################################################
  
  signatures <- apply(plot_df[, set_names, drop = FALSE], 1, paste0, collapse = "")
  
  inter_df <- as.data.frame(table(signatures), stringsAsFactors = FALSE)
  colnames(inter_df) <- c("signature", "count")
  inter_df$count <- as.numeric(inter_df$count)
  
  zero_sig <- paste(rep(0, length(set_names)), collapse = "")
  
  inter_df <- inter_df %>%
    dplyr::filter(signature != zero_sig, count > 0)
  
  if (nrow(inter_df) == 0) {
    stop("No non-empty intersections found for: ", filename_base)
  }
  
  inter_df$degree <- sapply(
    strsplit(inter_df$signature, ""),
    function(x) sum(x == "1")
  )
  
  inter_df <- inter_df %>%
    dplyr::arrange(dplyr::desc(count), dplyr::desc(degree), signature) %>%
    dplyr::slice_head(n = nintersects) %>%
    dplyr::mutate(inter_idx = seq_len(dplyr::n()))
  
  ##########################################################
  ## 3. Matrix data
  ##########################################################
  
  set_levels <- rev(set_names)
  
  mat_list <- lapply(seq_len(nrow(inter_df)), function(i) {
    sig_vec <- as.integer(strsplit(inter_df$signature[i], "")[[1]])
    data.frame(
      inter_idx = inter_df$inter_idx[i],
      set_name = set_names,
      active = sig_vec,
      stringsAsFactors = FALSE
    )
  })
  
  mat_df <- dplyr::bind_rows(mat_list)
  mat_df$set_name <- factor(mat_df$set_name, levels = set_levels)
  mat_df$y <- as.numeric(mat_df$set_name)
  
  seg_df <- mat_df %>%
    dplyr::filter(active == 1) %>%
    dplyr::group_by(inter_idx) %>%
    dplyr::summarise(
      ymin = min(y),
      ymax = max(y),
      .groups = "drop"
    ) %>%
    dplyr::filter(ymax > ymin)
  
  ##########################################################
  ## 4. Set-size data
  ##########################################################
  
  set_df <- data.frame(
    set_name = set_names,
    size = as.numeric(set_sizes[set_names]),
    stringsAsFactors = FALSE
  )
  
  set_df$set_name <- factor(set_df$set_name, levels = set_levels)
  set_df$y <- as.numeric(set_df$set_name)
  set_df$fill_col <- set_bar_colors[as.character(set_df$set_name)]
  
  max_set_size <- max(set_df$size, na.rm = TRUE)
  axis_max <- ceiling(max_set_size / 100) * 100
  if (!is.finite(axis_max) || axis_max <= 0) {
    axis_max <- max_set_size
  }
  
  bar_breaks <- pretty(c(0, axis_max), n = 4)
  bar_breaks <- bar_breaks[bar_breaks >= 0 & bar_breaks <= axis_max]
  bar_breaks <- sort(unique(c(bar_breaks, 0)))
  
  set_df <- set_df %>%
    dplyr::mutate(
      xmin = -size,
      xmax = 0,
      ymin = y - 0.18,
      ymax = y + 0.18
    )
  
  ##########################################################
  ## 5. Left-top legend panel
  ## Move rightward and reduce empty space
  ##########################################################
  
  legend_df <- data.frame(
    set_name = set_names,
    y = rev(seq_along(set_names)),
    x = 0.20,
    stringsAsFactors = FALSE
  )
  
  legend_df$fill_col <- set_bar_colors[legend_df$set_name]
  
  p_legend <- ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = set_name),
      shape = 15,
      size = 4.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.38, label = set_name),
      hjust = 0,
      size = 4.15,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_color_manual(values = set_bar_colors, guide = "none") +
    ggplot2::coord_cartesian(
      xlim = c(0.05, 2.85),
      ylim = c(0.5, length(set_names) + 0.5),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 0, 2, 2)
    )
  
  ##########################################################
  ## 6. Right-top intersection bar plot
  ## Remove y-axis title
  ##########################################################
  
  p_top <- ggplot2::ggplot(inter_df, ggplot2::aes(x = inter_idx, y = count)) +
    ggplot2::geom_col(
      fill = main_bar_color,
      width = 0.62
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = count),
      vjust = -0.25,
      size = 3.2,
      color = "#333333"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.5, max(inter_df$inter_idx) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.10))
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "#333333"),
      plot.margin = ggplot2::margin(8, 8, 4, 2)
    )
  
  ##########################################################
  ## 7. Right-bottom matrix plot
  ## Key fix:
  ##   draw grey points first -> draw black solid lines -> draw black points
  ##########################################################
  
  shade_df <- data.frame(
    y = seq_along(set_levels),
    ymin = seq_along(set_levels) - 0.5,
    ymax = seq_along(set_levels) + 0.5,
    shade = rep(c(TRUE, FALSE), length.out = length(set_levels))
  )
  
  p_matrix <- ggplot2::ggplot() +
    
    ggplot2::geom_rect(
      data = shade_df[shade_df$shade, ],
      ggplot2::aes(
        xmin = 0.5,
        xmax = max(inter_df$inter_idx) + 0.5,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "#F5F5F5",
      color = NA
    ) +
    
    ## 1. inactive grey dots first
    ggplot2::geom_point(
      data = mat_df,
      ggplot2::aes(x = inter_idx, y = y),
      size = 2.6,
      color = "#DDDDDD"
    ) +
    
    ## 2. solid lines second; they pass through intermediate grey dots
    ggplot2::geom_segment(
      data = seg_df,
      ggplot2::aes(
        x = inter_idx,
        xend = inter_idx,
        y = ymin,
        yend = ymax
      ),
      linewidth = 0.75,
      color = "#222222",
      lineend = "round"
    ) +
    
    ## 3. active black dots last
    ggplot2::geom_point(
      data = mat_df %>% dplyr::filter(active == 1),
      ggplot2::aes(x = inter_idx, y = y),
      size = 2.9,
      color = "#222222"
    ) +
    
    ggplot2::scale_x_continuous(
      limits = c(0.5, max(inter_df$inter_idx) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(set_levels),
      labels = rep("", length(set_levels)),
      limits = c(0.5, length(set_levels) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 8, 8, 2)
    )
  
  ##########################################################
  ## 8. Left-bottom set-size colored bars
  ## Move rightward, make compact, keep colored modules visible
  ##########################################################
  
  p_setsize <- ggplot2::ggplot() +
    
    ggplot2::geom_rect(
      data = set_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = set_name
      ),
      color = NA
    ) +
    
    ggplot2::scale_fill_manual(values = set_bar_colors, guide = "none") +
    
    ggplot2::scale_y_continuous(
      breaks = seq_along(set_levels),
      labels = rep("", length(set_levels)),
      limits = c(0.5, length(set_levels) + 0.5),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_x_continuous(
      limits = c(-axis_max * 1.03, axis_max * 0.08),
      breaks = -bar_breaks,
      labels = bar_breaks,
      expand = c(0, 0)
    ) +
    
    ggplot2::coord_cartesian(clip = "off") +
    
    ggplot2::labs(x = sets_x, y = NULL) +
    
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        face = "plain",
        size = base_size - 0.5,
        margin = ggplot2::margin(t = 2)
      ),
      axis.text.x = ggplot2::element_text(
        color = "#333333",
        size = base_size - 2,
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 0, 12, 2)
    )
  
  ##########################################################
  ## 9. Assemble with patchwork
  ## Make left column narrower and closer to the right
  ##########################################################
  
  left_col <- p_legend / p_setsize +
    patchwork::plot_layout(heights = c(2.0, 1.7))
  
  right_col <- p_top / p_matrix +
    patchwork::plot_layout(heights = c(2.0, 1.7))
  
  p_final <- (left_col | right_col) +
    patchwork::plot_layout(widths = c(1.45, 4.95))
  
  ##########################################################
  ## 10. Save
  ##########################################################
  
  ggplot2::ggsave(
    paste0(filename_base, ".pdf"),
    plot = p_final,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    dpi = 450,
    bg = "white"
  )
  
  ggplot2::ggsave(
    paste0(filename_base, ".png"),
    plot = p_final,
    width = width,
    height = height,
    dpi = 450,
    bg = "white"
  )
  
  invisible(p_final)
}
############################################################
## 5. Read input files
############################################################

gene_files <- c(
  Bulk_DEG    = "Bulk_DEG_genes.txt",
  Bulk_WGCNA  = "Bulk_WGCNA_genes.txt",
  Sc_marker   = "Sc_marker_genes.txt",
  Sc_hdWGCNA  = "Sc_hdWGCNA_genes.txt"
)

gene_sets <- lapply(gene_files, function(fn) {
  read_one_column_txt(find_file_recursive(GENE_DIR, fn), to_upper = TRUE)
})
names(gene_sets) <- names(gene_files)

tf_files <- c(
  CistromeDB_TF     = "CistromeDB_TF.txt",
  Motif_TF          = "Motif_TF.txt",
  HSC_likeMarkerTF  = "HSC_likeMarkerTF.txt"
)

tf_sets <- lapply(tf_files, function(fn) {
  read_one_column_txt(find_file_recursive(TF_DIR, fn), to_upper = TRUE)
})
names(tf_sets) <- names(tf_files)

drug_files <- c(
  TF_DGIdb     = "TF_DGIdb.txt",
  TF_PanDrugs  = "TF_PanDrugs.txt",
  BULK_EpiMed  = "BULK_EpiMed.txt",
  BULK_L2S2    = "BULK_L2S2.txt"
)

drug_sets <- lapply(drug_files, function(fn) {
  read_one_column_txt(find_file_recursive(DRUG_DIR, fn), to_upper = TRUE)
})
names(drug_sets) <- names(drug_files)

set_summary <- bind_rows(
  data.frame(
    category = "Gene integration",
    source = names(gene_sets),
    n = sapply(gene_sets, length)
  ),
  data.frame(
    category = "TF evidence",
    source = names(tf_sets),
    n = sapply(tf_sets, length)
  ),
  data.frame(
    category = "Drug evidence",
    source = names(drug_sets),
    n = sapply(drug_sets, length)
  )
)

write.csv(
  set_summary,
  file.path(OUT_DIR, "00_input_set_summary.csv"),
  row.names = FALSE
)

message("[INFO] Input set summary:")
print(set_summary)

############################################################
## 6. Fig3E: Bulk + single-cell integrated gene UpSet
############################################################

gene_df <- make_membership_df(gene_sets)

gene_out_prefix <- file.path(
  OUT_DIR,
  "Fig3E_Bulk_sc_integrated_gene"
)

write_membership_outputs(
  gene_df,
  names(gene_sets),
  gene_out_prefix
)

save_upset_plot(
  df = gene_df,
  set_names = c("Bulk_DEG", "Bulk_WGCNA", "Sc_marker", "Sc_hdWGCNA"),
  filename_base = file.path(OUT_DIR, "Fig3E_Bulk_sc_integrated_gene_UpSet"),
  mainbar_y = "Intersecting genes",
  sets_x = "Genes per source",
  width = 7.2,
  height = 4.8,
  nintersects = 12,
  main_bar_color = "#355C5C",
  set_bar_colors = fig3_source_cols[c("Bulk_DEG", "Bulk_WGCNA", "Sc_marker", "Sc_hdWGCNA")]
)

############################################################
## 7. Fig3F: TF evidence Venn + optional UpSet
############################################################

tf_df <- make_membership_df(tf_sets)

tf_out_prefix <- file.path(
  OUT_DIR,
  "Fig3F_TF_evidence"
)

write_membership_outputs(
  tf_df,
  names(tf_sets),
  tf_out_prefix
)

## 7A. Three-set Venn
tf_venn_cols <- c("#7080C0", "#F0E866", "#7AC36A")

pdf(
  file.path(OUT_DIR, "Fig3F_TF_evidence_3set_Venn.pdf"),
  width = 5.2,
  height = 4.8,
  useDingbats = FALSE
)
grid.newpage()
venn_obj <- VennDiagram::venn.diagram(
  x = tf_sets,
  filename = NULL,
  fill = tf_venn_cols,
  alpha = 0.58,
  col = "grey35",
  lwd = 1.0,
  cex = 1.20,
  fontface = "bold",
  cat.cex = 1.05,
  cat.fontface = "bold",
  cat.col = c("#253494", "#8C8600", "#1B7837"),
  margin = 0.06
)
grid.draw(venn_obj)
dev.off()

png(
  file.path(OUT_DIR, "Fig3F_TF_evidence_3set_Venn.png"),
  width = 5.2 * 450,
  height = 4.8 * 450,
  res = 450
)
grid.newpage()
venn_obj <- VennDiagram::venn.diagram(
  x = tf_sets,
  filename = NULL,
  fill = tf_venn_cols,
  alpha = 0.58,
  col = "grey35",
  lwd = 1.0,
  cex = 1.20,
  fontface = "bold",
  cat.cex = 1.05,
  cat.fontface = "bold",
  cat.col = c("#253494", "#8C8600", "#1B7837"),
  margin = 0.06
)
grid.draw(venn_obj)
dev.off()

## 7B. Optional TF UpSet
save_upset_plot(
  df = tf_df,
  set_names = c("CistromeDB_TF", "Motif_TF", "HSC_likeMarkerTF"),
  filename_base = file.path(OUT_DIR, "Fig3F_TF_evidence_UpSet_optional"),
  mainbar_y = "Intersecting TFs",
  sets_x = "TFs per source",
  width = 6.8,
  height = 4.4,
  nintersects = 8,
  main_bar_color = "#4A4A4A",
  set_bar_colors = fig3_source_cols[c("CistromeDB_TF", "Motif_TF", "HSC_likeMarkerTF")]
)
############################################################
## 8. Fig3G: Drug evidence UpSet
############################################################

drug_df <- make_membership_df(drug_sets)

drug_out_prefix <- file.path(
  OUT_DIR,
  "Fig3G_drug_evidence"
)

write_membership_outputs(
  drug_df,
  names(drug_sets),
  drug_out_prefix
)

save_upset_plot(
  df = drug_df,
  set_names = c("TF_DGIdb", "TF_PanDrugs", "BULK_EpiMed", "BULK_L2S2"),
  filename_base = file.path(OUT_DIR, "Fig3G_drug_evidence_UpSet"),
  mainbar_y = "Intersecting candidates",
  sets_x = "Candidates per source",
  width = 7.2,
  height = 4.8,
  nintersects = 14,
  main_bar_color = "#404040",
  set_bar_colors = fig3_source_cols[c("TF_DGIdb", "TF_PanDrugs", "BULK_EpiMed", "BULK_L2S2")]
)

############################################################
## 9. Fig3H: Drug evidence matrix
############################################################

plot_evidence_matrix <- function(df,
                                 set_names,
                                 filename_base,
                                 top_n = 30,
                                 min_evidence = 2,
                                 show_title = FALSE) {
  
  mat_df <- df %>%
    dplyr::filter(Evidence_count >= min_evidence) %>%
    dplyr::arrange(dplyr::desc(Evidence_count), item)
  
  if (nrow(mat_df) == 0) {
    mat_df <- df %>%
      dplyr::arrange(dplyr::desc(Evidence_count), item) %>%
      dplyr::slice_head(n = top_n)
  } else {
    mat_df <- mat_df %>%
      dplyr::slice_head(n = top_n)
  }
  
  mat <- as.matrix(mat_df[, set_names, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- mat_df$item
  
  ann_row <- data.frame(
    Evidence_count = mat_df$Evidence_count
  )
  rownames(ann_row) <- mat_df$item
  
  ann_col <- data.frame(
    Source = factor(set_names, levels = set_names)
  )
  rownames(ann_col) <- set_names
  
  source_cols_use <- drug_source_cols[set_names]
  source_cols_use[is.na(source_cols_use)] <- "#BDBDBD"
  
  ann_colors <- list(
    Source = source_cols_use,
    Evidence_count = evidence_count_cols
  )
  
  fig_h <- max(4.8, min(8.0, 0.13 * nrow(mat) + 1.8))
  fig_w <- 5.4
  
  plot_title <- if (show_title) "Drug / target evidence matrix" else ""
  
  pheatmap::pheatmap(
    mat,
    color = c("#FFFFFF", "#252525"),
    breaks = c(-0.1, 0.5, 1.1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8.2,
    fontsize_row = ifelse(nrow(mat) > 35, 5.6, 6.4),
    fontsize_col = 8.8,
    angle_col = 45,
    border_color = "#D9D9D9",
    legend = FALSE,
    main = plot_title,
    filename = paste0(filename_base, ".pdf"),
    width = fig_w,
    height = fig_h
  )
  
  pheatmap::pheatmap(
    mat,
    color = c("#FFFFFF", "#252525"),
    breaks = c(-0.1, 0.5, 1.1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8.2,
    fontsize_row = ifelse(nrow(mat) > 35, 5.6, 6.4),
    fontsize_col = 8.8,
    angle_col = 45,
    border_color = "#D9D9D9",
    legend = FALSE,
    main = plot_title,
    filename = paste0(filename_base, ".png"),
    width = fig_w,
    height = fig_h
  )
  
  write.csv(
    mat_df,
    paste0(filename_base, "_plotted_candidates.csv"),
    row.names = FALSE
  )
  
  invisible(mat_df)
}

plot_evidence_matrix(
  df = drug_df,
  set_names = names(drug_sets),
  filename_base = file.path(
    OUT_DIR,
    "Fig3H_drug_evidence_matrix_top_candidates"
  ),
  top_n = 30,
  min_evidence = 2,
  show_title = FALSE
)

############################################################
## 10. Highlight important candidates
############################################################

important_names <- c(
  "GNL3L", "GNL3", "MYC", "TP53", "SPI1", "IRF1", "STAT3", "MAX",
  "RBPJ", "BACH1", "GATA1",
  "VENETOCLAX", "AZACITIDINE", "DECITABINE", "CYTARABINE",
  "BORTEZOMIB", "SIROLIMUS", "VORINOSTAT", "PALBOCICLIB",
  "FENOFIBRATE"
)

highlight_gene <- gene_df %>%
  dplyr::filter(item %in% important_names) %>%
  dplyr::mutate(category = "Gene integration")

highlight_tf <- tf_df %>%
  dplyr::filter(item %in% important_names) %>%
  dplyr::mutate(category = "TF evidence")

highlight_drug <- drug_df %>%
  dplyr::filter(item %in% important_names) %>%
  dplyr::mutate(category = "Drug evidence")

write.csv(
  dplyr::bind_rows(highlight_gene, highlight_tf, highlight_drug),
  file.path(OUT_DIR, "09_highlighted_candidates_across_sources.csv"),
  row.names = FALSE
)

############################################################
## 11. Finished
############################################################

message("\n[DONE] Fig3 overlap/evidence plots have been exported to:")
message(OUT_DIR)

message("\n[KEY OUTPUTS]")
message("1) Fig3E_Bulk_sc_integrated_gene_UpSet.pdf/png")
message("2) Fig3F_TF_evidence_3set_Venn.pdf/png")
message("3) Fig3F_TF_evidence_UpSet_optional.pdf/png")
message("4) Fig3G_drug_evidence_UpSet.pdf/png")
message("5) Fig3H_drug_evidence_matrix_top_candidates.pdf/png")
message("6) CSV files for membership matrices and multi-evidence candidates")
# END ORIGINAL SCRIPT P21

# ==============================================================================
# BEGIN ORIGINAL SCRIPT P22: 14_TF_CistromeDB_decoupleR.R
# ==============================================================================
################################################################################
##
##  GNL3L TF regulation analysis - figure optimized full version (v2_no_motif)
##  Based on the user's currently runnable full pipeline.
##
##  Core analysis logic preserved except motif scanning removed:
##    1) CistromeDB candidate TFs
##    2) decoupleR ULM TF activity using CollecTRI
##       - expression matrix from RNA layer=data
##       - TF activity matrix built from score column
##    3) FindMarkers(HSC-like vs Others) on tfsulm assay
##    4) Intersection of CistromeDB candidate TFs and differential TF activity
##       to define putative regulators
##
##  Figure updates in this version:
##    - Volcano plots: Scheme A / A2 / A3 all exported
##    - Motif scanning and motif-related plots removed
##
################################################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(qs)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(decoupleR)
  library(RColorBrewer)
  library(pheatmap)
  library(scales)
  library(grid)
  library(gridExtra)
  library(Matrix)
})

# =============================================================================
# 0. Global settings
# =============================================================================
WORK_DIR <- "D:/0.DataAnalysis/AML_SC_RNA/GSE116256/1.16AML+5NORMAL"
setwd(WORK_DIR)

OUTPUT_DIR <- file.path(WORK_DIR, "11.GNL3L_TF分析/0.new")
STEP1_DIR <- file.path(OUTPUT_DIR, "Step1_CistromeDB")
STEP2_DIR <- file.path(OUTPUT_DIR, "Step2_TF_Activity")
STEP3_DIR <- file.path(OUTPUT_DIR, "Step3_Differential_TF_Activity")
STEP4_DIR <- file.path(OUTPUT_DIR, "Step4_Putative_Regulators")
FIG_DIR   <- file.path(OUTPUT_DIR, "Publication_Ready_Figures")
for (d in c(OUTPUT_DIR, STEP1_DIR, STEP2_DIR, STEP3_DIR, STEP4_DIR, FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

CISTROME_FILE <- "9.transcription_factor_decoupleR/CistromeDB_GNL3L_TF.xlsx"
SEURAT_FILE   <- "2.Sc_annotation/new/2.resolution_0.5_annotation.qs"
REGULON_CACHE_FILE <- file.path(OUTPUT_DIR, "CollecTRI_human_splitFALSE.rds")


CELLTYPE_COL <- "celltype"
HSC_LABEL    <- "HSC-like"
PVAL_THRESH  <- 0.05
LOGFC_THRESH <- 0.5
set.seed(42)

COL_HSC_FILL   <- "#CB5148"
COL_HSC_EDGE   <- "#A13E37"
COL_OTHER_FILL <- "#6A9BCB"
COL_OTHER_EDGE <- "#4F7EAA"
COL_NS_FILL    <- "#B3B3B3"
COL_NS_EDGE    <- "#8E8E8E"
BG_COL         <- "white"
COL_CISTROME   <- "#7A5195"
COL_CORE       <- "#CB5148"

# =============================================================================
# 1. Helpers
# =============================================================================
stop_if_missing <- function(path_string) {
  if (!file.exists(path_string)) stop("Required file not found: ", path_string)
}

check_required_columns <- function(df, cols, df_name = "data.frame") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop(df_name, " is missing required columns: ", paste(miss, collapse = ", "))
}

safe_ggsave <- function(filename, plot_obj, width, height, dpi = 320, bg = "white") {
  ext <- tolower(tools::file_ext(filename))
  if (ext == "pdf") {
    if (capabilities("cairo")) {
      ggsave(filename, plot = plot_obj, width = width, height = height, device = cairo_pdf, bg = bg)
    } else {
      ggsave(filename, plot = plot_obj, width = width, height = height, device = grDevices::pdf, bg = bg)
    }
  } else {
    ggsave(filename, plot = plot_obj, width = width, height = height, dpi = dpi, bg = bg)
  }
}

save_plot_pdf_png <- function(plot_obj, pdf_file, png_file, width = 8, height = 6, dpi = 320, bg = "white") {
  safe_ggsave(pdf_file, plot_obj, width, height, dpi, bg)
  safe_ggsave(png_file, plot_obj, width, height, dpi, bg)
}

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0.5, color = "gray35", size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(8, 10, 8, 10)
    )
}

prepare_collectri_original_semantics <- function(cache_file) {
  message("Preparing CollecTRI regulons ...")
  cache_file <- normalizePath(cache_file, winslash = "/", mustWork = FALSE)
  message("Cache path: ", cache_file)
  
  if (file.exists(cache_file)) {
    message("Using local cached CollecTRI regulons.")
    regs <- readRDS(cache_file)
  } else {
    message("Local cache not found. Running decoupleR::get_collectri() ...")
    regs <- suppressWarnings(
      decoupleR::get_collectri(
        organism = "human",
        split_complexes = FALSE
      )
    )
    check_required_columns(regs, c("source", "target", "mor"), "collectri_net")
    saveRDS(regs, cache_file)
    message("CollecTRI cache saved to: ", cache_file)
  }
  
  check_required_columns(regs, c("source", "target", "mor"), "collectri_net")
  regs
}

# -------------------------------
# Volcano: scheme A (existing style)
# -------------------------------
make_tf_activity_volcano_A <- function(diff_df,
                                       fc_cut = 0.5,
                                       p_cut = 0.05,
                                       title = "HSC-like vs Others",
                                       left_label = "Higher in Others",
                                       right_label = "Higher in HSC-like",
                                       label_top_n = 8,
                                       x_display_limit = 6) {
  check_required_columns(diff_df, c("TF", "avg_log2FC", "p_val"), "diff_df")
  
  vol_df <- diff_df %>%
    mutate(
      Log2FC_raw = as.numeric(avg_log2FC),
      pvalue = as.numeric(p_val),
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      pvalue = ifelse(pvalue == 0, 1e-300, pvalue),
      neglogp = -log10(pvalue),
      regulated = case_when(
        pvalue < p_cut & Log2FC_raw >  fc_cut ~ "Higher in HSC-like",
        pvalue < p_cut & Log2FC_raw < -fc_cut ~ "Higher in Others",
        TRUE ~ "Not significant"
      ),
      Log2FC_plot = pmax(pmin(Log2FC_raw, x_display_limit), -x_display_limit),
      is_clipped = abs(Log2FC_raw) > x_display_limit
    )
  
  total_n <- nrow(vol_df)
  right_n <- sum(vol_df$regulated == "Higher in HSC-like")
  left_n  <- sum(vol_df$regulated == "Higher in Others")
  y_max <- ceiling(max(vol_df$neglogp, na.rm = TRUE))
  y_max <- max(y_max, 10)
  
  label_df <- bind_rows(
    vol_df %>% filter(regulated == "Higher in HSC-like") %>% arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>% slice_head(n = label_top_n),
    vol_df %>% filter(regulated == "Higher in Others") %>% arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>% slice_head(n = label_top_n)
  ) %>% distinct(TF, .keep_all = TRUE)
  
  ggplot(vol_df, aes(x = Log2FC_plot, y = neglogp)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = BG_COL) +
    geom_hline(yintercept = -log10(p_cut), color = "grey65", linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), color = "grey55", linetype = 2, linewidth = 0.9) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
    geom_point(data = subset(vol_df, regulated == "Not significant"), shape = 21, size = 2.2, stroke = 0.2, color = COL_NS_EDGE, fill = COL_NS_FILL, alpha = 0.75) +
    geom_point(data = subset(vol_df, regulated == "Higher in Others"), shape = 21, size = 2.8, stroke = 0.25, color = COL_OTHER_EDGE, fill = COL_OTHER_FILL, alpha = 0.95) +
    geom_point(data = subset(vol_df, regulated == "Higher in HSC-like"), shape = 21, size = 2.8, stroke = 0.25, color = COL_HSC_EDGE, fill = COL_HSC_FILL, alpha = 0.95) +
    # clipped triangles as visual cue
    geom_point(data = subset(vol_df, is_clipped & Log2FC_raw < -x_display_limit), aes(x = -x_display_limit), shape = 24, size = 2.8, stroke = 0.25, color = COL_OTHER_EDGE, fill = COL_OTHER_FILL) +
    geom_point(data = subset(vol_df, is_clipped & Log2FC_raw > x_display_limit), aes(x = x_display_limit), shape = 25, size = 2.8, stroke = 0.25, color = COL_HSC_EDGE, fill = COL_HSC_FILL) +
    ggrepel::geom_text_repel(data = label_df, aes(label = TF), size = 3.3, fontface = "bold", box.padding = 0.35, point.padding = 0.20, segment.color = "grey45", segment.size = 0.30, min.segment.length = 0, max.overlaps = Inf, seed = 42) +
    annotate("segment", x = -0.15, xend = -x_display_limit, y = y_max * 0.90, yend = y_max * 0.90, colour = COL_OTHER_FILL, linewidth = 1.0, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("segment", x = 0.15, xend = x_display_limit, y = y_max * 0.90, yend = y_max * 0.90, colour = COL_HSC_FILL, linewidth = 1.0, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("text", x = -x_display_limit * 0.62, y = y_max * 0.92, label = paste0(round(left_n / total_n * 100, 1), " %  ", left_label), colour = COL_OTHER_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = x_display_limit * 0.62, y = y_max * 0.92, label = paste0(right_label, "  ", round(right_n / total_n * 100, 1), " %"), colour = COL_HSC_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = x_display_limit, y = y_max * 0.04, label = paste0("|x| > ", x_display_limit, " truncated for display"), hjust = 1, size = 3.2, color = "grey35") +
    scale_x_continuous(limits = c(-x_display_limit, x_display_limit), breaks = seq(-x_display_limit, x_display_limit, by = 1)) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = expression(Log[2]~foldchange), y = expression(-Log[10]~italic(P)~"-value")) +
    theme_classic(base_size = 13) +
    theme(plot.background = element_rect(fill = BG_COL, color = NA), panel.background = element_rect(fill = BG_COL, color = NA), axis.line = element_line(linewidth = 1.0, color = "black"), axis.ticks = element_line(linewidth = 0.9, color = "black"), axis.title = element_text(face = "bold", size = 17), axis.text = element_text(size = 11), plot.title = element_text(face = "bold", hjust = 0.5, size = 14), plot.margin = margin(10, 16, 10, 16))
}

# -------------------------------
# Volcano: scheme A2 (clipped circles, no triangles)
# -------------------------------
make_tf_activity_volcano_A2 <- function(diff_df,
                                     fc_cut = 0.5,
                                     p_cut = 0.05,
                                     title = "HSC-like vs Others",
                                     left_label = "Higher in Others",
                                     right_label = "Higher in HSC-like",
                                     label_top_n = 8,
                                     x_display_limit = 6) {
  check_required_columns(diff_df, c("TF", "avg_log2FC", "p_val"), "diff_df")

  vol_df <- diff_df %>%
    dplyr::mutate(
      Log2FC_raw = as.numeric(avg_log2FC),
      pvalue = as.numeric(p_val),
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      pvalue = ifelse(pvalue == 0, 1e-300, pvalue),
      neglogp = -log10(pvalue),
      regulated = dplyr::case_when(
        pvalue < p_cut & Log2FC_raw >  fc_cut ~ "Higher in HSC-like",
        pvalue < p_cut & Log2FC_raw < -fc_cut ~ "Higher in Others",
        TRUE ~ "Not significant"
      ),
      # 仅用于显示，不改变真实筛选结果
      Log2FC_plot = pmax(pmin(Log2FC_raw, x_display_limit), -x_display_limit)
    )

  total_n <- nrow(vol_df)
  right_n <- sum(vol_df$regulated == "Higher in HSC-like")
  left_n  <- sum(vol_df$regulated == "Higher in Others")

  y_max <- ceiling(max(vol_df$neglogp, na.rm = TRUE))
  y_max <- max(y_max, 10)

  label_df <- dplyr::bind_rows(
    vol_df %>%
      dplyr::filter(regulated == "Higher in HSC-like") %>%
      dplyr::arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>%
      dplyr::slice_head(n = label_top_n),
    vol_df %>%
      dplyr::filter(regulated == "Higher in Others") %>%
      dplyr::arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>%
      dplyr::slice_head(n = label_top_n)
  ) %>%
    dplyr::distinct(TF, .keep_all = TRUE)

  p <- ggplot(vol_df, aes(x = Log2FC_plot, y = neglogp)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = BG_COL) +
    geom_hline(yintercept = -log10(p_cut), color = "grey65", linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), color = "grey55", linetype = 2, linewidth = 1.0) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +

    geom_point(
      data = subset(vol_df, regulated == "Not significant"),
      shape = 21, size = 2.2, stroke = 0.20,
      color = COL_NS_EDGE, fill = COL_NS_FILL, alpha = 0.75
    ) +
    geom_point(
      data = subset(vol_df, regulated == "Higher in Others"),
      shape = 21, size = 2.8, stroke = 0.25,
      color = COL_OTHER_EDGE, fill = COL_OTHER_FILL, alpha = 0.95
    ) +
    geom_point(
      data = subset(vol_df, regulated == "Higher in HSC-like"),
      shape = 21, size = 2.8, stroke = 0.25,
      color = COL_HSC_EDGE, fill = COL_HSC_FILL, alpha = 0.95
    ) +

    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = TF),
      size = 3.3,
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.20,
      segment.color = "grey45",
      segment.size = 0.30,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 42
    ) +

    annotate("segment",
             x = -0.15, xend = -x_display_limit,
             y = y_max * 0.90, yend = y_max * 0.90,
             colour = COL_OTHER_FILL, linewidth = 1.0,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("segment",
             x = 0.15, xend = x_display_limit,
             y = y_max * 0.90, yend = y_max * 0.90,
             colour = COL_HSC_FILL, linewidth = 1.0,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +

    annotate("text",
             x = -x_display_limit * 0.62, y = y_max * 0.92,
             label = paste0(round(left_n / total_n * 100, 1), " %  ", left_label),
             colour = COL_OTHER_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text",
             x = x_display_limit * 0.62, y = y_max * 0.92,
             label = paste0(right_label, "  ", round(right_n / total_n * 100, 1), " %"),
             colour = COL_HSC_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +

    scale_x_continuous(
      limits = c(-x_display_limit, x_display_limit),
      breaks = seq(-x_display_limit, x_display_limit, by = 1)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = expression(Log[2]~foldchange),
      y = expression(-Log[10]~italic(P)~"-value")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = BG_COL, color = NA),
      panel.background = element_rect(fill = BG_COL, color = NA),
      axis.line = element_line(linewidth = 1.0, color = "black"),
      axis.ticks = element_line(linewidth = 0.9, color = "black"),
      axis.title = element_text(face = "bold", size = 17),
      axis.text = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.margin = margin(10, 16, 10, 16)
    )

  return(p)
}
# -------------------------------
# Volcano: scheme A3 (quantile-based clipping)
# -------------------------------
make_tf_activity_volcano_A3 <- function(diff_df,
                                        fc_cut = 0.5,
                                        p_cut = 0.05,
                                        title = "HSC-like vs Others",
                                        left_label = "Higher in Others",
                                        right_label = "Higher in HSC-like",
                                        label_top_n = 6,
                                        lower_quantile = 0.01,
                                        upper_quantile = 0.99) {
  check_required_columns(diff_df, c("TF", "avg_log2FC", "p_val"), "diff_df")
  
  raw_fc <- as.numeric(diff_df$avg_log2FC)
  q_low  <- as.numeric(stats::quantile(raw_fc, probs = lower_quantile, na.rm = TRUE))
  q_high <- as.numeric(stats::quantile(raw_fc, probs = upper_quantile, na.rm = TRUE))
  x_display_limit <- max(abs(q_low), abs(q_high))
  x_display_limit <- max(2, ceiling(x_display_limit))
  
  vol_df <- diff_df %>%
    mutate(
      Log2FC_raw = as.numeric(avg_log2FC),
      pvalue = as.numeric(p_val),
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      pvalue = ifelse(pvalue == 0, 1e-300, pvalue),
      neglogp = -log10(pvalue),
      regulated = case_when(
        pvalue < p_cut & Log2FC_raw >  fc_cut ~ "Higher in HSC-like",
        pvalue < p_cut & Log2FC_raw < -fc_cut ~ "Higher in Others",
        TRUE ~ "Not significant"
      ),
      Log2FC_plot = pmax(pmin(Log2FC_raw, x_display_limit), -x_display_limit)
    )
  
  total_n <- nrow(vol_df)
  right_n <- sum(vol_df$regulated == "Higher in HSC-like")
  left_n  <- sum(vol_df$regulated == "Higher in Others")
  y_max <- ceiling(max(vol_df$neglogp, na.rm = TRUE))
  y_max <- max(y_max, 10)
  
  label_df <- bind_rows(
    vol_df %>% filter(regulated == "Higher in HSC-like") %>% arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>% slice_head(n = label_top_n),
    vol_df %>% filter(regulated == "Higher in Others") %>% arrange(desc(neglogp), desc(abs(Log2FC_raw))) %>% slice_head(n = label_top_n)
  ) %>% distinct(TF, .keep_all = TRUE)
  
  ggplot(vol_df, aes(x = Log2FC_plot, y = neglogp)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = BG_COL) +
    geom_hline(yintercept = -log10(p_cut), color = "grey65", linetype = 3, linewidth = 0.8) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), color = "grey45", linetype = 2, linewidth = 1.0) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
    geom_point(data = subset(vol_df, regulated == "Not significant"), shape = 21, size = 2.2, stroke = 0.2, color = COL_NS_EDGE, fill = COL_NS_FILL, alpha = 0.75) +
    geom_point(data = subset(vol_df, regulated == "Higher in Others"), shape = 21, size = 2.8, stroke = 0.25, color = COL_OTHER_EDGE, fill = COL_OTHER_FILL, alpha = 0.95) +
    geom_point(data = subset(vol_df, regulated == "Higher in HSC-like"), shape = 21, size = 2.8, stroke = 0.25, color = COL_HSC_EDGE, fill = COL_HSC_FILL, alpha = 0.95) +
    ggrepel::geom_text_repel(data = label_df, aes(label = TF), size = 3.2, fontface = "bold", box.padding = 0.32, point.padding = 0.18, segment.color = "grey45", segment.size = 0.30, min.segment.length = 0, max.overlaps = Inf, seed = 42) +
    annotate("segment", x = -0.15, xend = -x_display_limit, y = y_max * 0.90, yend = y_max * 0.90, colour = COL_OTHER_FILL, linewidth = 1.0, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("segment", x = 0.15, xend = x_display_limit, y = y_max * 0.90, yend = y_max * 0.90, colour = COL_HSC_FILL, linewidth = 1.0, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("text", x = -x_display_limit * 0.62, y = y_max * 0.92, label = paste0(round(left_n / total_n * 100, 1), " %  ", left_label), colour = COL_OTHER_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = x_display_limit * 0.62, y = y_max * 0.92, label = paste0(right_label, "  ", round(right_n / total_n * 100, 1), " %"), colour = COL_HSC_FILL, fontface = "bold", size = 5.0, hjust = 0.5) +
    annotate("text", x = x_display_limit, y = y_max * 0.04, label = paste0("display range based on ", lower_quantile*100, "%–", upper_quantile*100, "% quantiles"), hjust = 1, size = 3.0, color = "grey35") +
    scale_x_continuous(limits = c(-x_display_limit, x_display_limit), breaks = pretty(c(-x_display_limit, x_display_limit), n = 7)) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = expression(Log[2]~foldchange), y = expression(-Log[10]~italic(P)~"-value")) +
    theme_classic(base_size = 13) +
    theme(plot.background = element_rect(fill = BG_COL, color = NA), panel.background = element_rect(fill = BG_COL, color = NA), axis.line = element_line(linewidth = 1.0, color = "black"), axis.ticks = element_line(linewidth = 0.9, color = "black"), axis.title = element_text(face = "bold", size = 17), axis.text = element_text(size = 11), plot.title = element_text(face = "bold", hjust = 0.5, size = 14), plot.margin = margin(10, 16, 10, 16))
}

make_tf_workflow_schematic <- function(cistrome_n, diff_n, core_n) {
  flow_df <- data.frame(
    step = factor(c("CistromeDB\\nquery", "Differential TF\\nactivity", "Final putative\\nregulators"),
                  levels = c("CistromeDB\\nquery", "Differential TF\\nactivity", "Final putative\\nregulators")),
    x = c(1, 3, 5),
    y = 1,
    fill = c("#EEE7F7", "#FDEAE8", "#F7D9D6"),
    lab = c(
      paste0("CistromeDB\\n", cistrome_n, " candidate TFs"),
      paste0("FindMarkers on TF activity\\n", diff_n, " differential TFs\\n(p < 0.05, |log2FC| > 0.5)"),
      paste0("Intersection\\n", core_n, " putative regulators")
    )
  )
  arrow_df <- data.frame(x = c(1.45, 3.45), xend = c(2.55, 4.55), y = 1, yend = 1)

  ggplot() +
    geom_segment(data = arrow_df, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                 linewidth = 0.9, color = "gray35") +
    geom_label(data = flow_df, aes(x = x, y = y, label = lab, fill = step),
               label.size = 0.35, color = "black", fontface = "bold", lineheight = 1.0,
               size = 4.0, label.r = unit(0.18, "lines"), label.padding = unit(0.35, "lines")) +
    scale_fill_manual(values = setNames(flow_df$fill, levels(flow_df$step))) +
    annotate("text", x = c(2, 4), y = c(1.26, 1.26),
             label = c("public TF-gene\\nevidence", "cell-state-specific\\nactivity support"),
             color = "gray35", size = 3.3) +
    labs(title = "Integrative screening pipeline for upstream TFs of GNL3L") +
    coord_cartesian(xlim = c(0.45, 5.55), ylim = c(0.78, 1.34), clip = "off") +
    theme_void(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5, size = 14, margin = margin(b = 10)))
}

# =============================================================================
# 2. Step 1. CistromeDB
# =============================================================================
cat("Step 1: CistromeDB\n")
stop_if_missing(CISTROME_FILE)

cistrome_raw <- read_excel(CISTROME_FILE)
colnames(cistrome_raw) <- gsub("_", "", trimws(colnames(cistrome_raw)))
bad_idx <- is.na(colnames(cistrome_raw)) | colnames(cistrome_raw) == ""
if (any(bad_idx)) colnames(cistrome_raw)[bad_idx] <- paste0("X", which(bad_idx))
if ("Factor" %in% colnames(cistrome_raw)) colnames(cistrome_raw)[colnames(cistrome_raw) == "Factor"] <- "TF"
if ("RPscore" %in% colnames(cistrome_raw)) colnames(cistrome_raw)[colnames(cistrome_raw) == "RPscore"] <- "RP_score"
check_required_columns(cistrome_raw, c("TF", "RP_score"), "cistrome_raw")

cistrome_summary <- cistrome_raw %>%
  mutate(TF = toupper(TF)) %>%
  group_by(TF) %>%
  summarise(
    N_Peaks = n(),
    Max_RP = max(RP_score, na.rm = TRUE),
    Mean_RP = mean(RP_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_RP))

CistromeDB_TFs <- cistrome_summary$TF
write.csv(cistrome_raw, file.path(STEP1_DIR, "CistromeDB_Raw_Data.csv"), row.names = FALSE)
write.csv(cistrome_summary, file.path(STEP1_DIR, "CistromeDB_TF_Summary.csv"), row.names = FALSE)
writeLines(CistromeDB_TFs, file.path(STEP1_DIR, "CistromeDB_TF_List.txt"))

# =============================================================================
# 3. Step 2. TF activity with decoupleR (original semantics)
# =============================================================================
cat("Step 2: TF activity calculation\n")
stop_if_missing(SEURAT_FILE)

seurat_obj <- qread(SEURAT_FILE)
cat("【数据概览】\n")
cat("  总细胞数:", ncol(seurat_obj), "\n")
cat("  细胞类型分布:\n")
print(table(seurat_obj@meta.data[[CELLTYPE_COL]]))
cat("\n")

expr_mat <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data"))

cat("【decoupleR参数】\n")
cat("  方法: ULM\n")
cat("  网络: CollecTRI\n")
cat("  minsize: 5 (默认)\n\n")

collectri_net <- prepare_collectri_original_semantics(REGULON_CACHE_FILE)
cat("  CollecTRI调控关系数:", nrow(collectri_net), "\n")
cat("  运行decoupleR...\n")

tf_activities <- run_ulm(
  mat = expr_mat,
  net = collectri_net,
  .source = "source",
  .target = "target",
  .mor = "mor"
)

check_required_columns(tf_activities, c("source", "condition", "score", "statistic"), "tf_activities")

tf_activity_mat <- tf_activities %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

mode(tf_activity_mat) <- "numeric"

cat("\n【Step2 结果】\n")
cat("  计算得到TF数:", nrow(tf_activity_mat), "\n")

CistromeDB_TFs_with_activity <- intersect(CistromeDB_TFs, rownames(tf_activity_mat))
cat("  CistromeDB_TFs中可计算活性的:", length(CistromeDB_TFs_with_activity), "/", length(CistromeDB_TFs), "\n\n")

if (!all(colnames(seurat_obj) %in% colnames(tf_activity_mat))) {
  miss <- setdiff(colnames(seurat_obj), colnames(tf_activity_mat))
  stop("Some Seurat cells are missing from tf_activity_mat: ", paste(head(miss, 10), collapse = ", "))
}

tf_activity_mat_ordered <- tf_activity_mat[, colnames(seurat_obj), drop = FALSE]
mode(tf_activity_mat_ordered) <- "numeric"

seurat_obj[["tfsulm"]] <- CreateAssayObject(data = tf_activity_mat_ordered)

saveRDS(tf_activity_mat, file.path(STEP2_DIR, "TF_Activity_Matrix.rds"))
qsave(seurat_obj, file.path(STEP2_DIR, "Seurat_with_TF_Activity.qs"))
qsave(seurat_obj, file.path(OUTPUT_DIR, "Seurat_with_TF_Activity.qs"))

celltypes <- unique(seurat_obj@meta.data[[CELLTYPE_COL]])
mean_activity_by_ct <- data.frame()
for (ct in celltypes) {
  cells <- colnames(seurat_obj)[seurat_obj@meta.data[[CELLTYPE_COL]] == ct]
  ct_activity <- tf_activity_mat[, cells, drop = FALSE]
  
  mean_act <- data.frame(
    TF = rownames(ct_activity),
    CellType = ct,
    Mean_Activity = rowMeans(ct_activity, na.rm = TRUE),
    SD_Activity = apply(ct_activity, 1, sd, na.rm = TRUE),
    N_Cells = length(cells)
  )
  mean_activity_by_ct <- rbind(mean_activity_by_ct, mean_act)
}

mean_activity_wide <- mean_activity_by_ct %>%
  select(TF, CellType, Mean_Activity) %>%
  pivot_wider(names_from = CellType, values_from = Mean_Activity) %>%
  column_to_rownames("TF") %>%
  as.matrix()

write.csv(mean_activity_by_ct, file.path(STEP2_DIR, "TF_Activity_by_CellType.csv"), row.names = FALSE)
write.csv(data.frame(TF = CistromeDB_TFs_with_activity), file.path(STEP2_DIR, "CistromeDB_TFs_with_activity.csv"), row.names = FALSE)

if (length(CistromeDB_TFs_with_activity) > 0) {
  hm_data <- mean_activity_wide[CistromeDB_TFs_with_activity, , drop = FALSE]
  hm_data_scaled <- t(scale(t(hm_data)))
  pdf(file.path(STEP2_DIR, "CistromeDB_TF_Activity_Heatmap.pdf"), width = 10, height = max(5, 0.3 * nrow(hm_data_scaled) + 2))
  pheatmap::pheatmap(hm_data_scaled,
                     color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     fontsize_row = 8,
                     fontsize_col = 10,
                     border_color = NA,
                     main = "CistromeDB candidate TF activity across cell types")
  dev.off()
}

# =============================================================================
# 4. Step 3. Differential TF activity (original semantics)
# =============================================================================
cat("Step 3: Differential TF activity\n")
DefaultAssay(seurat_obj) <- "tfsulm"
Idents(seurat_obj) <- seurat_obj@meta.data[[CELLTYPE_COL]]

cat("【差异TF活性分析参数】\n")
cat("  比较: HSC-like vs Others\n")
cat("  方法: Wilcoxon rank-sum test\n")
cat("  筛选标准: p<", PVAL_THRESH, " & |log2FC|>", LOGFC_THRESH, "\n\n", sep = "")

diff_tf_activity <- FindMarkers(
  seurat_obj,
  ident.1 = HSC_LABEL,
  assay = "tfsulm"
)

required_step3_cols <- c("p_val", "avg_log2FC")
missing_step4_cols <- setdiff(required_step3_cols, colnames(diff_tf_activity))
if (length(missing_step4_cols) > 0) {
  stop(
    paste0(
      "FindMarkers output does not contain the original expected columns: ",
      paste(missing_step4_cols, collapse = ", "),
      ".\nCurrent columns are: ",
      paste(colnames(diff_tf_activity), collapse = ", "),
      ".\nThis means the current Step2 object state does not match the original TF.R behavior."
    )
  )
}

diff_tf_activity <- diff_tf_activity %>%
  rownames_to_column("TF") %>%
  mutate(
    In_CistromeDB = TF %in% CistromeDB_TFs,
    Category = case_when(
      p_val < PVAL_THRESH & avg_log2FC > LOGFC_THRESH ~ "Upregulated",
      p_val < PVAL_THRESH & avg_log2FC < -LOGFC_THRESH ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  arrange(p_val)

cat("【Step3 结果】\n")
cat("  总差异TF数:", nrow(diff_tf_activity), "\n")

diff_up <- diff_tf_activity %>% filter(Category == "Upregulated")
diff_down <- diff_tf_activity %>% filter(Category == "Downregulated")

cat("  上调TF活性 (p<0.05 & log2FC>0.5):", nrow(diff_up), "\n")
cat("  下调TF活性 (p<0.05 & log2FC<-0.5):", nrow(diff_down), "\n\n")

Diff_TFs_up <- diff_up$TF
Diff_TFs_down <- diff_down$TF

write.csv(diff_tf_activity, file.path(STEP3_DIR, "Differential_TF_Activity_All.csv"), row.names = FALSE)
write.csv(diff_up, file.path(STEP3_DIR, "Upregulated_TF_Activity.csv"), row.names = FALSE)
write.csv(diff_down, file.path(STEP3_DIR, "Downregulated_TF_Activity.csv"), row.names = FALSE)
if (length(Diff_TFs_up) > 0) writeLines(Diff_TFs_up, file.path(STEP3_DIR, "Upregulated_TF_List.txt"))
if (length(Diff_TFs_down) > 0) writeLines(Diff_TFs_down, file.path(STEP3_DIR, "Downregulated_TF_List.txt"))

pdf(file.path(STEP3_DIR, "Differential_TF_Activity_Volcano.pdf"), width = 10, height = 8)
p_step4 <- ggplot(diff_tf_activity, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(aes(color = Category), alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = -log10(PVAL_THRESH), linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = "dashed", color = "gray50", linewidth = 0.6) +
  scale_color_manual(values = c("Upregulated" = COL_HSC_FILL,
                                "Downregulated" = COL_OTHER_FILL,
                                "Not Significant" = COL_NS_FILL)) +
  labs(title = "Differential TF Activity in HSC-like vs Others",
       x = "Average log2 fold change",
       y = "-log10(p-value)",
       color = NULL) +
  theme_publication(13)
print(p_step4)
dev.off()

# Volcano exports: A / A2 / A3
p_volcano_A <- make_tf_activity_volcano_A(diff_tf_activity, fc_cut = LOGFC_THRESH, p_cut = PVAL_THRESH, title = "HSC-like vs Others", left_label = "Higher in Others", right_label = "Higher in HSC-like", label_top_n = 8, x_display_limit = 6)
save_plot_pdf_png(p_volcano_A,
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA.pdf"),
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA.png"),
                  width = 9.8, height = 6.7, dpi = 320, bg = BG_COL)

p_volcano_A2 <- make_tf_activity_volcano_A2(diff_tf_activity, fc_cut = LOGFC_THRESH, p_cut = PVAL_THRESH, title = "HSC-like vs Others", left_label = "Higher in Others", right_label = "Higher in HSC-like", label_top_n = 6, x_display_limit = 6)
save_plot_pdf_png(p_volcano_A2,
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA2.pdf"),
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA2.png"),
                  width = 9.8, height = 6.7, dpi = 320, bg = BG_COL)

p_volcano_A3 <- make_tf_activity_volcano_A3(diff_tf_activity, fc_cut = LOGFC_THRESH, p_cut = PVAL_THRESH, title = "HSC-like vs Others", left_label = "Higher in Others", right_label = "Higher in HSC-like", label_top_n = 6, lower_quantile = 0.01, upper_quantile = 0.99)
save_plot_pdf_png(p_volcano_A3,
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA3.pdf"),
                  file.path(FIG_DIR, "Fig_TFActivity_Volcano_HSC_like_vs_Others_schemeA3.png"),
                  width = 9.8, height = 6.7, dpi = 320, bg = BG_COL)

# =============================================================================
# 5. Step 4. Final putative regulators
# =============================================================================
cat("Step 4: Putative regulators\n")
positive_regulators <- intersect(CistromeDB_TFs, Diff_TFs_up)
negative_regulators <- intersect(CistromeDB_TFs, Diff_TFs_down)
all_regulators <- c(positive_regulators, negative_regulators)

positive_info <- NULL
negative_info <- NULL

if (length(positive_regulators) > 0) {
  positive_info <- diff_up %>%
    filter(TF %in% positive_regulators) %>%
    left_join(cistrome_summary %>% select(TF, N_Peaks, Max_RP), by = "TF") %>%
    mutate(Regulation = "Putative Positive Regulator") %>%
    arrange(p_val)
  write.csv(positive_info, file.path(STEP4_DIR, "Putative_Positive_Regulators.csv"), row.names = FALSE)
  writeLines(positive_regulators, file.path(STEP4_DIR, "Putative_Positive_Regulator_List.txt"))
}

if (length(negative_regulators) > 0) {
  negative_info <- diff_down %>%
    filter(TF %in% negative_regulators) %>%
    left_join(cistrome_summary %>% select(TF, N_Peaks, Max_RP), by = "TF") %>%
    mutate(Regulation = "Putative Negative Regulator") %>%
    arrange(p_val)
  write.csv(negative_info, file.path(STEP4_DIR, "Putative_Negative_Regulators.csv"), row.names = FALSE)
  writeLines(negative_regulators, file.path(STEP4_DIR, "Putative_Negative_Regulator_List.txt"))
}

all_info <- bind_rows(positive_info, negative_info)
if (nrow(all_info) > 0) write.csv(all_info, file.path(STEP4_DIR, "All_Putative_Regulators.csv"), row.names = FALSE)

filtering_summary <- data.frame(
  Step = c("CistromeDB Query", "Differential TF Activity (Up)", "Differential TF Activity (Down)", "Putative Positive Regulators", "Putative Negative Regulators"),
  N_TFs = c(length(CistromeDB_TFs), length(Diff_TFs_up), length(Diff_TFs_down), length(positive_regulators), length(negative_regulators))
)
write.csv(filtering_summary, file.path(STEP4_DIR, "Filtering_Summary.csv"), row.names = FALSE)

# =============================================================================
# 6. Publication-ready figures
# =============================================================================
cat("Publication-ready figures\n")
p_workflow <- make_tf_workflow_schematic(length(CistromeDB_TFs), nrow(diff_up) + nrow(diff_down), length(all_regulators))
save_plot_pdf_png(p_workflow,
                  file.path(FIG_DIR, "Fig_Workflow_GNL3L_TF_screening.pdf"),
                  file.path(FIG_DIR, "Fig_Workflow_GNL3L_TF_screening.png"),
                  width = 11, height = 3.6, dpi = 320, bg = "white")

cat("
Done. Output directory:\n", OUTPUT_DIR, "\n")
cat("Original-semantics TF activity object saved to:\n", file.path(STEP2_DIR, "Seurat_with_TF_Activity.qs"), "\n")
cat("Additional figure outputs saved in:\n", FIG_DIR, "\n")

# END ORIGINAL SCRIPT P22
