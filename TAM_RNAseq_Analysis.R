# Original R scripts: verbatim compilation
# Only section boundary comments were added; source-script bytes are unchanged.
# Read README.md before use. This file is not a one-command pipeline.
# Original paths, source() calls, historical alternatives and environment resets remain.


# ==============================================================================
# BEGIN ORIGINAL SCRIPT T01: 00_config.R
# ==============================================================================
## =============================================================================
## 06_analysis_revised / 00_config.R
## 基于 03_analysis(final clean workflow v1.0) 修订：
##   - 从原始 count 重新读取、完整重跑（不复用 03 的 RDS 结果）
##   - 01 数据准备/DEG、02 共识：代码逻辑沿用 03，不作方法学改动
##   - 新增基于 DEG 的 ORA（GO-BP/CC/MF + KEGG + Hallmark）
##   - GSEA：Hallmark 沿用 03；KEGG 用 msigdbr 2026.1 的 CP:KEGG_LEGACY 修复补全
##   - signature：保留 HSC-like + 5 个 SC hdWGCNA 模块 hub，新增 BULK WGCNA
##     turquoise core(114)；删除 AML_UP/AML_DOWN/GNL3L_POS/GNL3L_NEG；不做 TF 网络
## 绘图风格/配色严格沿用 03（theme_journal、Arial、cairo_pdf）。
## =============================================================================
options(stringsAsFactors = FALSE)

PROJECT <- "D:/Users/lab01/Desktop/GNL3L_TAM小论文"
RAW <- file.path(PROJECT, "03_药物细胞实验与测序数据/02_药物干预后转录组测序/01_rawdata")
OUT <- Sys.getenv("TAM_ANALYSIS_ROOT", unset = normalizePath(getwd(), winslash = "/"))
PUB <- file.path(PROJECT, "02_公共数据分析")

DIR <- list(
  r = file.path(OUT, "R"), metadata = file.path(OUT, "00_metadata"),
  intermediate = file.path(OUT, "01_intermediate"), deg = file.path(OUT, "02_DEG"),
  ora = file.path(OUT, "03_ORA"), gsea = file.path(OUT, "04_GSEA"),
  sig = file.path(OUT, "05_signature"), tables = file.path(OUT, "06_tables"),
  main = file.path(OUT, "07_main"), supp = file.path(OUT, "08_supplementary"),
  logs = file.path(OUT, "09_logs")
)
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))

PKG <- c("data.table", "edgeR", "limma", "ggplot2", "dplyr", "tidyr", "tibble",
         "stringr", "ggrepel", "ComplexHeatmap", "circlize", "patchwork",
         "readxl", "fgsea", "msigdbr", "scales", "matrixStats",
         "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
missing_pkg <- PKG[!vapply(PKG, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) stop("Missing packages: ", paste(missing_pkg, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(PKG, library, character.only = TRUE)))

## ---- 阈值（沿用 03 + 本次 ORA 决策）----
COL <- c(Control = "#4D4D4D", TAM5 = "#E69F00", TAM12 = "#D55E00",
         Up = "#D73027", Down = "#4575B4", NS = "#BDBDBD")
CELL_ORDER <- c("MV4-11", "THP-1", "MOLM-13", "SKM-1")
CONTRAST_ORDER <- as.vector(t(outer(CELL_ORDER, c("TAM5", "TAM12"), paste, sep = "_")))
DEG_LFC <- 1
DEG_P <- 0.05            # 作者锁定：DEG 展示阈值 |log2FC|>1 且 nominal P<0.05（沿用 03）
FDR_GSEA <- 0.05         # GSEA FDR（沿用 03 口径）
GSEA_MIN <- 10; GSEA_MAX <- 500          # Hallmark/KEGG GSEA 基因集大小（沿用 03）
SIG_MIN <- 10; SIG_MAX <- 1000           # signature GSEA 基因集大小（沿用 03）
ORA_P <- 0.05            # ORA 显著性：按决策使用原始 pvalue<0.05（非 p.adjust）
ORA_TOPN <- 10           # 每个数据库展示 top10
FOCUS_GENES <- c("GNL3L", "MYC", "BAX", "BCL2", "CASP3", "PARP1", "TP53", "SPI1")

theme_journal <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold"),
      strip.background = element_blank(), strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(3, "mm")
    )
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
  ggplot2::ggsave(filename, plot, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf, bg = "white")
}

read_tab_matrix <- function(path) {
  x <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  names(x)[1] <- "ensembl_id"
  x$ensembl_id <- sub("\\..*$", "", x$ensembl_id)
  x
}

annotation_map <- function(paths) {
  ans <- lapply(paths, function(p) {
    x <- data.table::fread(p, sep = "\t", select = c("Gene id", "Gene name"),
                          data.table = FALSE, check.names = FALSE)
    names(x) <- c("ensembl_id", "gene_symbol")
    x$ensembl_id <- sub("\\..*$", "", x$ensembl_id)
    x$gene_symbol[is.na(x$gene_symbol)] <- ""
    x
  })
  dplyr::bind_rows(ans) |>
    dplyr::filter(ensembl_id != "") |>
    dplyr::arrange(ensembl_id, dplyr::desc(gene_symbol != "")) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)
}

make_metadata <- function() {
  add <- function(batch, cell, dose, reps = 1:3, source_dose = dose) {
    src_cell <- gsub("-", "_", cell)
    data.frame(
      batch = batch, cell_line = cell, dose = dose, replicate = reps,
      source_sample_id = paste(src_cell, source_dose, reps, sep = "_"),
      analysis_sample_id = paste(batch, src_cell, dose, reps, sep = "__")
    )
  }
  b1 <- dplyr::bind_rows(lapply(CELL_ORDER, function(z)
    dplyr::bind_rows(add("B1_202604", z, "Control", source_dose = "CON"),
                     add("B1_202604", z, "TAM5"))))
  b2 <- dplyr::bind_rows(lapply(CELL_ORDER, function(z) {
    # 厂商将 THP-1 的 TAM12 列误标为 TAM5；保留 source id，在 metadata 层纠正 dose。
    src <- if (z == "THP-1") "TAM5" else "TAM12"
    dplyr::bind_rows(add("B2_202605", z, "Control", source_dose = "CON"),
                     add("B2_202605", z, "TAM12", source_dose = src))
  }))
  dplyr::bind_rows(b1, b2) |>
    dplyr::mutate(cell_line = factor(cell_line, levels = CELL_ORDER),
                  dose = factor(dose, levels = c("Control", "TAM5", "TAM12")),
                  contrast_id = ifelse(dose == "Control", NA_character_,
                                       paste(cell_line, dose, sep = "_")))
}

META <- make_metadata()

PATH <- list(
  b1_count = file.path(RAW, "202604第一批测序/gene.count.matrix.xls"),
  b1_annot = file.path(RAW, "202604第一批测序/gene.TPM.matrix.annot.xls"),
  b2_skm_count = file.path(RAW, "202605第二批测序/1.先送6个样本_SKM-1与TAM12处理后/ExpressMatrix/gene.count.matrix.xls"),
  b2_skm_annot = file.path(RAW, "202605第二批测序/1.先送6个样本_SKM-1与TAM12处理后/ExpressMatrix/gene.TPM.matrix.annot.xls"),
  b2_other_count = file.path(RAW, "202605第二批测序/2.后送18个样本_MV4-11_THP-1_MOLM-13与TAM12处理后/ExpressMatrix/gene.count.matrix.xls"),
  b2_other_annot = file.path(RAW, "202605第二批测序/2.后送18个样本_MV4-11_THP-1_MOLM-13与TAM12处理后/ExpressMatrix/gene.TPM.matrix.annot.xls"),
  ## 修订后 signature 来源（全部来自前期公共数据分析）
  hsc_marker = file.path(PUB, "10_SC_enrich/Tab01_celltype_markers.csv"),
  hd_hub     = file.path(PUB, "12_SC_hdWGCNA/Tab02_hubGene.csv"),
  bulk_wgcna_core = file.path(PUB, "05_BULK_WGCNA/WGCNA_GSE30029_turquoise_core_gene_list.txt")
)

write.csv(META, file.path(DIR$metadata, "sample_metadata_48.csv"), row.names = FALSE)

# END ORIGINAL SCRIPT T01

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T02: 01_prepare_and_DEG.R
# ==============================================================================
## =============================================================================
## 01_prepare_and_DEG.R  —— 沿用 03_analysis，逻辑不做方法学改动，从原始 count 重跑
## 读三批 count -> 三批交集基因 -> 48 样本矩阵 -> 8 组批次内匹配 limma-voom DEG
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

stopifnot(all(file.exists(unlist(PATH[1:6]))))
b1 <- read_tab_matrix(PATH$b1_count)
b2s <- read_tab_matrix(PATH$b2_skm_count)
b2o <- read_tab_matrix(PATH$b2_other_count)
annot <- annotation_map(c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot))
write.csv(annot, file.path(DIR$metadata, "ensembl_gene_symbol_map.csv"), row.names = FALSE)

select_batch <- function(x, meta_batch) {
  ids <- as.character(meta_batch$source_sample_id)
  miss <- setdiff(ids, names(x))
  if (length(miss)) stop("Samples absent from source matrix: ", paste(miss, collapse = ", "))
  y <- x[, c("ensembl_id", ids), drop = FALSE]
  names(y)[-1] <- as.character(meta_batch$analysis_sample_id)
  y
}

x1 <- select_batch(b1, META |> dplyr::filter(batch == "B1_202604"))
x2s <- select_batch(b2s, META |> dplyr::filter(batch == "B2_202605", cell_line == "SKM-1"))
x2o <- select_batch(b2o, META |> dplyr::filter(batch == "B2_202605", cell_line != "SKM-1"))

# 三批取交集，避免把厂商过滤掉的缺失基因在其它批次误当 0
common_ids <- Reduce(intersect, list(x1$ensembl_id, x2s$ensembl_id, x2o$ensembl_id))
writeLines(common_ids, file.path(DIR$metadata, "three_delivery_common_ensembl_ids.txt"))
merge_common <- function(x) x[match(common_ids, x$ensembl_id), , drop = FALSE]
all_counts_df <- Reduce(function(a, b) merge(a, b, by = "ensembl_id", all = FALSE),
                        lapply(list(x1, x2s, x2o), merge_common))
all_counts <- as.matrix(all_counts_df[, -1, drop = FALSE]); mode(all_counts) <- "numeric"
rownames(all_counts) <- all_counts_df$ensembl_id
all_counts <- all_counts[, as.character(META$analysis_sample_id), drop = FALSE]
saveRDS(all_counts, file.path(DIR$intermediate, "counts_48_common_genes.rds"))

# 全局 logCPM 仅用于 QC/热图可视化，不用于跨批次 raw-dose 检验
dge_all <- edgeR::DGEList(all_counts)
keep_all <- edgeR::filterByExpr(dge_all, group = interaction(META$cell_line, META$dose, META$batch))
dge_all <- edgeR::calcNormFactors(dge_all[keep_all, , keep.lib.sizes = FALSE])
logcpm_all <- edgeR::cpm(dge_all, log = TRUE, prior.count = 2)
saveRDS(logcpm_all, file.path(DIR$intermediate, "logCPM_48_common_genes.rds"))

get_source <- function(cell, dose) {
  if (dose == "TAM5") return(x1)
  if (cell == "SKM-1") x2s else x2o
}

run_contrast <- function(cell_i, dose_i) {
  m <- META |> dplyr::filter(.data$cell_line == .env$cell_i,
                            .data$batch == ifelse(.env$dose_i == "TAM5", "B1_202604", "B2_202605"),
                            .data$dose %in% c("Control", .env$dose_i)) |>
    dplyr::arrange(dose, replicate)
  x <- get_source(cell_i, dose_i)
  ids <- as.character(m$analysis_sample_id)
  missing_ids <- setdiff(ids, names(x))
  if (length(missing_ids)) {
    stop("Standardized samples absent from selected matrix for ", cell_i, "_", dose_i,
         ": ", paste(missing_ids, collapse = ", "))
  }
  mat <- as.matrix(x[, ids, drop = FALSE]); mode(mat) <- "numeric"
  rownames(mat) <- x$ensembl_id; colnames(mat) <- ids
  grp <- factor(as.character(m$dose), levels = c("Control", dose_i))
  dge <- edgeR::DGEList(mat)
  keep <- edgeR::filterByExpr(dge, group = grp)
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
  v <- limma::voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  cm <- limma::makeContrasts(contrasts = paste0(dose_i, "-Control"), levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, cm), robust = TRUE)
  tt <- limma::topTable(fit, number = Inf, sort.by = "P") |>
    tibble::rownames_to_column("ensembl_id") |>
    dplyr::left_join(annot, by = "ensembl_id") |>
    dplyr::mutate(display_gene = ifelse(is.na(gene_symbol) | gene_symbol == "", ensembl_id, gene_symbol),
                  direction = dplyr::case_when(
                    logFC > DEG_LFC & P.Value < DEG_P ~ "Up",
                    logFC < -DEG_LFC & P.Value < DEG_P ~ "Down", TRUE ~ "NS"),
                  cell_line = cell_i, dose = dose_i, contrast_id = paste(cell_i, dose_i, sep = "_"))
  od <- file.path(DIR$deg, paste(cell_i, dose_i, sep = "_")); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  write.csv(tt, file.path(od, "DEG_all.csv"), row.names = FALSE)
  write.csv(tt |> dplyr::filter(direction != "NS"), file.path(od, "DEG_nominalP0.05_absLFC1.csv"), row.names = FALSE)
  saveRDS(list(voom = v, fit = fit, metadata = m, keep = keep),
          file.path(od, "limma_voom_object.rds"))
  tt
}

deg_list <- list()
for (cell in CELL_ORDER) for (dose in c("TAM5", "TAM12")) {
  id <- paste(cell, dose, sep = "_"); message("Running ", id)
  deg_list[[id]] <- run_contrast(cell, dose)
}
saveRDS(deg_list, file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))
writeLines(capture.output(sessionInfo()), file.path(DIR$logs, "sessionInfo_01.txt"))
message("01 done: ", length(deg_list), " contrasts.")

# END ORIGINAL SCRIPT T02

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T03: 02_consensus_analysis.R
# ==============================================================================
## =============================================================================
## 02_consensus_analysis.R —— 沿用 03_analysis，不改方法（只分析，不出图）
## 规则法：同方向 >=3/4 株 且 DEG 状态 >=2/4 株
## 探索性 Stouffer meta-z 整合：meta_p<0.05 且同方向 >=3/4 株，组内 BH
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))
long <- dplyr::bind_rows(deg) |>
  dplyr::select(ensembl_id, gene_symbol, display_gene, cell_line, dose,
                contrast_id, logFC, t, P.Value, adj.P.Val)

safe_symbol <- function(x) {
  z <- x[!is.na(x) & x != ""]
  if (length(z)) z[1] else NA_character_
}

## 锁定的规则法定义
consensus_rule <- long |>
  dplyr::group_by(dose, ensembl_id) |>
  dplyr::summarise(
    gene_symbol = safe_symbol(gene_symbol),
    display_gene = dplyr::first(display_gene),
    n_pos = sum(logFC > 0), n_neg = sum(logFC < 0),
    n_sig = sum(abs(logFC) > DEG_LFC & P.Value < DEG_P),
    median_logFC = stats::median(logFC),
    mean_abs_logFC = mean(abs(logFC)), .groups = "drop"
  ) |>
  dplyr::filter(pmax(n_pos, n_neg) >= 3, n_sig >= 2)

## 探索性 Stouffer 整合（与既往实现保持一致）
meta_tbl <- long |>
  dplyr::mutate(z = sign(logFC) * stats::qnorm(pmax(P.Value, 1e-300) / 2,
                                               lower.tail = FALSE)) |>
  dplyr::group_by(dose, ensembl_id) |>
  dplyr::summarise(
    gene_symbol = safe_symbol(gene_symbol),
    display_gene = dplyr::first(display_gene),
    meta_z = sum(z) / sqrt(dplyr::n()),
    meta_p = 2 * stats::pnorm(-abs(meta_z)),
    n_pos = sum(logFC > 0), n_neg = sum(logFC < 0),
    median_logFC = stats::median(logFC), .groups = "drop"
  ) |>
  dplyr::group_by(dose) |>
  dplyr::mutate(meta_fdr = stats::p.adjust(meta_p, "BH")) |>
  dplyr::ungroup() |>
  dplyr::filter(meta_p < 0.05, pmax(n_pos, n_neg) >= 3)

write.csv(consensus_rule,
          file.path(DIR$tables, "Fig6D_consensus_rule_3of4_and_2significant.csv"),
          row.names = FALSE)
write.csv(meta_tbl, file.path(DIR$tables, "Fig6D_consensus_stouffer_meta.csv"),
          row.names = FALSE)
message("02 done: rule=", nrow(consensus_rule), " stouffer=", nrow(meta_tbl))

# END ORIGINAL SCRIPT T03

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T04: 03_ORA_GO_KEGG_Hallmark.R
# ==============================================================================
## =============================================================================
## 03_ORA_GO_KEGG_Hallmark.R  —— 【本次新增】基于 DEG 的过表达分析(ORA)
## 决策固化：
##   - 每个比较组 DEG = Up+Down 合并（不分方向），阈值沿用 |log2FC|>1 & nominal P<0.05
##   - 背景 universe = 该比较组 filterByExpr 后进入检验的全部(有 symbol)基因
##   - GO-BP/CC/MF：clusterProfiler::enrichGO(org.Hs.eg.db, ENTREZID, BH)
##   - KEGG / Hallmark：msigdbr 基因集(CP:KEGG_LEGACY / H) + clusterProfiler::enricher
##                     （离线超几何，与 GSEA 同源，不依赖 KEGG 在线接口）
##   - 显著性按【原始 pvalue < 0.05】（非 p.adjust）；每库按 pvalue 升序取 top10 展示
## 本脚本只出表；ORA 图在 06 出图脚本统一绘制（沿用 03 的分析/出图分离架构）
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))

## msigdbr term2gene（symbol），兼容新版 collection/subcollection 参数
msig_t2g <- function(collection, subcollection = NULL) {
  args <- list(species = "Homo sapiens")
  f <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% f) args$collection <- collection else args$category <- collection
  if (!is.null(subcollection)) {
    if ("subcollection" %in% f) args$subcollection <- subcollection else args$subcategory <- subcollection
  }
  do.call(msigdbr::msigdbr, args) |>
    dplyr::select(term = gs_name, gene = gene_symbol)
}
hall_t2g <- msig_t2g("H")
kegg_t2g <- msig_t2g("C2", "CP:KEGG_LEGACY")
message("Hallmark gene sets: ", dplyr::n_distinct(hall_t2g$term),
        " | KEGG_LEGACY gene sets: ", dplyr::n_distinct(kegg_t2g$term))

valid_sym <- function(x) unique(x[!is.na(x) & x != "" & !grepl("^ENSG", x)])

as_df <- function(eobj) {
  if (is.null(eobj)) return(data.frame())
  z <- as.data.frame(eobj)
  if (is.null(z) || !nrow(z)) return(data.frame())
  z
}

run_ora_one <- function(tt, id) {
  bits <- strsplit(id, "_", fixed = TRUE)[[1]]
  dose_i <- tail(bits, 1); cell_i <- sub(paste0("_", dose_i, "$"), "", id)

  uni_sym <- valid_sym(tt$gene_symbol)                         # 全部进入检验的基因(背景)
  deg_df <- tt |> dplyr::filter(direction != "NS")
  deg_sym <- valid_sym(deg_df$gene_symbol)                     # 合并上下调后的 DEG
  n_up <- sum(deg_df$direction == "Up"); n_dn <- sum(deg_df$direction == "Down")

  ## SYMBOL -> ENTREZID（仅 GO 需要）
  bm <- clusterProfiler::bitr(uni_sym, "SYMBOL", "ENTREZID", org.Hs.eg.db, drop = TRUE)
  deg_ent <- unique(bm$ENTREZID[bm$SYMBOL %in% deg_sym])
  uni_ent <- unique(bm$ENTREZID)

  out <- list()

  ## GO: BP / CC / MF
  for (ont in c("BP", "CC", "MF")) {
    ego <- tryCatch(suppressMessages(
      clusterProfiler::enrichGO(deg_ent, universe = uni_ent, OrgDb = org.Hs.eg.db,
                                keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH",
                                pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE)),
      error = function(e) NULL)
    ego <- tryCatch(clusterProfiler::setReadable(ego, OrgDb = org.Hs.eg.db,
                                                 keyType = "ENTREZID"), error = function(e) ego)
    z <- as_df(ego)
    if (nrow(z)) { z$collection <- paste0("GO_", ont); out[[paste0("GO_", ont)]] <- z }
  }

  ## KEGG_LEGACY / Hallmark：symbol 级 enricher
  for (nm in c("KEGG", "Hallmark")) {
    t2g <- if (nm == "KEGG") kegg_t2g else hall_t2g
    ee <- tryCatch(clusterProfiler::enricher(
      deg_sym, universe = uni_sym, TERM2GENE = t2g,
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), error = function(e) NULL)
    z <- as_df(ee)
    if (nrow(z)) { z$collection <- nm; out[[nm]] <- z }
  }

  res <- dplyr::bind_rows(out)
  if (!nrow(res)) {
    return(data.frame(contrast_id = id, cell_line = cell_i, dose = dose_i,
                      n_deg = length(deg_sym), n_up = n_up, n_down = n_dn,
                      n_universe = length(uni_sym)))
  }
  res |>
    dplyr::mutate(contrast_id = id, cell_line = cell_i, dose = dose_i,
                  n_deg = length(deg_sym), n_up = n_up, n_down = n_dn,
                  n_universe = length(uni_sym),
                  collection = factor(collection,
                    levels = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Hallmark"))) |>
    dplyr::select(contrast_id, cell_line, dose, collection,
                  ID, Description, GeneRatio, BgRatio, pvalue, p.adjust,
                  qvalue, Count, geneID, n_deg, n_up, n_down, n_universe)
}

ora_list <- lapply(CONTRAST_ORDER, function(id) {
  message("ORA ", id)
  run_ora_one(deg[[id]], id)
})
ora_all <- dplyr::bind_rows(ora_list) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                cell_line = factor(cell_line, levels = CELL_ORDER))

## GeneRatio 解析为数值，便于排序/绘图
parse_ratio <- function(x) sapply(strsplit(as.character(x), "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
ora_all$GeneRatioNum <- parse_ratio(ora_all$GeneRatio)

write.csv(ora_all, file.path(DIR$ora, "ORA_all_long.csv"), row.names = FALSE)

## 显著（原始 pvalue<0.05）
ora_sig <- ora_all |> dplyr::filter(pvalue < ORA_P)
write.csv(ora_sig, file.path(DIR$tables, "ORA_sig_p0.05.csv"), row.names = FALSE)

## 每组每库按 pvalue 升序取 top10（且 p<0.05）用于展示
ora_top <- ora_sig |>
  dplyr::arrange(contrast_id, collection, pvalue) |>
  dplyr::group_by(contrast_id, collection) |>
  dplyr::slice_head(n = ORA_TOPN) |>
  dplyr::ungroup()
write.csv(ora_top, file.path(DIR$tables, "ORA_sig_p0.05_top10.csv"), row.names = FALSE)
saveRDS(ora_all, file.path(DIR$intermediate, "ORA_all_long.rds"))

## 留痕：每组 DEG/universe 规模 + 各库显著 term 数
trace <- ora_all |>
  dplyr::group_by(contrast_id, collection) |>
  dplyr::summarise(n_sig_p05 = sum(pvalue < ORA_P), .groups = "drop")
size_trace <- ora_all |>
  dplyr::distinct(contrast_id, n_deg, n_up, n_down, n_universe)
write.csv(size_trace, file.path(DIR$logs, "ORA_01_input_size_trace.csv"), row.names = FALSE)
write.csv(trace, file.path(DIR$logs, "ORA_02_sig_term_count_trace.csv"), row.names = FALSE)
message("03 ORA done. total terms=", nrow(ora_all), " sig(p<0.05)=", nrow(ora_sig),
        " display top10 rows=", nrow(ora_top))

# END ORIGINAL SCRIPT T04

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T05: 04_GSEA_Hallmark_KEGG.R
# ==============================================================================
## =============================================================================
## 04_GSEA_Hallmark_KEGG.R —— GSEA（沿用 03 方法，修复并补全 KEGG）
##   - 排名统计量：limma moderated t（同 symbol 取 |t| 最大一条），降序
##   - Hallmark：msigdbr H（与 03 一致）
##   - KEGG：msigdbr C2 / CP:KEGG_LEGACY（修复 03 中 CP:KEGG 取不到导致 KEGG 缺失）
##   - GO 不做 GSEA（按决策）
##   - fgseaMultilevel，minSize=10 / maxSize=500（沿用 03）
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))
deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))

get_msig <- function(collection, subcollection = NULL) {
  f <- names(formals(msigdbr::msigdbr))
  args <- list(species = "Homo sapiens")
  if ("collection" %in% f) args$collection <- collection else args$category <- collection
  if (!is.null(subcollection)) {
    if ("subcollection" %in% f) args$subcollection <- subcollection else args$subcategory <- subcollection
  }
  ms <- do.call(msigdbr::msigdbr, args) |> dplyr::select(gs_name, gene_symbol)
  split(ms$gene_symbol, ms$gs_name)
}

rank_symbols <- function(x) {
  y <- x |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "", is.finite(t)) |>
    dplyr::group_by(gene_symbol) |>
    dplyr::slice_max(abs(t), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  sort(stats::setNames(y$t, y$gene_symbol), decreasing = TRUE)
}

run_sets <- function(pathways, label) {
  out <- lapply(names(deg), function(id) {
    fgsea::fgseaMultilevel(pathways = pathways, stats = rank_symbols(deg[[id]]),
                           minSize = GSEA_MIN, maxSize = GSEA_MAX) |>
      as.data.frame() |>
      dplyr::mutate(
        leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
        contrast_id = id, collection = label)
  }) |> dplyr::bind_rows()
  write.csv(out, file.path(DIR$gsea, paste0("GSEA_", label, "_all.csv")), row.names = FALSE)
  write.csv(out, file.path(DIR$tables, paste0("GSEA_", label, "_all.csv")), row.names = FALSE)
  out
}

hallmark <- get_msig("H")
kegg_legacy <- get_msig("C2", "CP:KEGG_LEGACY")
message("Hallmark sets=", length(hallmark), "  KEGG_LEGACY sets=", length(kegg_legacy))

gh <- run_sets(hallmark, "Hallmark")
gk <- run_sets(kegg_legacy, "KEGG_LEGACY")
message("04 GSEA done. Hallmark rows=", nrow(gh), " KEGG rows=", nrow(gk))

# END ORIGINAL SCRIPT T05

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T06: 05_signature_GSEA_revised.R
# ==============================================================================
## =============================================================================
## 05_signature_GSEA_revised.R —— 【本次修订】7 个 signature 的逆转 GSEA
## 保留：HSC_LIKE_MARKERS(405) + SC hdWGCNA 5 模块 hub(M2/M4/M10/M11/M12, 各30)
## 新增：BULK_WGCNA_TURQUOISE_CORE(114)
## 删除：AML_UP / AML_DOWN / GNL3L_POS_TOP200 / GNL3L_NEG_BOTTOM200；不做 TF 网络
## 7 个均为 AML 恶性/HSC-like 干性/增殖相关程序，统一 orientation = -1
##   （负 NES = 被 TAM 抑制 = 逆转；reversal_score = NES * orientation，正=逆转）
## fgseaMultilevel，minSize=10 / maxSize=1000（沿用 03 signature 参数），排名量 = t
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))
deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))

clean_vec <- function(x) unique(x[!is.na(x) & x != ""])

## 1) HSC-like markers（cluster == HSC-like 的 gene 列）
mk <- data.table::fread(PATH$hsc_marker, data.table = FALSE)
hsc_genes <- clean_vec(mk$gene[mk$cluster == "HSC-like"])

## 2) BULK WGCNA turquoise core（每行一个基因）
bulk_core <- clean_vec(trimws(readLines(PATH$bulk_wgcna_core, warn = FALSE)))

## 3) SC hdWGCNA 5 个模块的 hub 基因（module 取值为 HSC-like-Mx，各 30）
hub <- data.table::fread(PATH$hd_hub, data.table = FALSE)
hd_modules <- c("HSC-like-M2", "HSC-like-M4", "HSC-like-M10", "HSC-like-M11", "HSC-like-M12")
hd_list <- lapply(hd_modules, function(m) clean_vec(hub$gene_name[hub$module == m]))
names(hd_list) <- paste0("SC_hdWGCNA_", sub("HSC-like-", "", hd_modules))

sig_list <- c(list(HSC_LIKE_MARKERS = hsc_genes,
                   BULK_WGCNA_TURQUOISE_CORE = bulk_core),
              hd_list)

## signature 定义/方向表（留痕）
sig_def <- data.frame(
  pathway = names(sig_list),
  n_gene = vapply(sig_list, length, integer(1)),
  orientation = -1,
  source = c("10_SC_enrich/Tab01_celltype_markers.csv (cluster=HSC-like)",
             "05_BULK_WGCNA/WGCNA_GSE30029_turquoise_core_gene_list.txt",
             rep("12_SC_hdWGCNA/Tab02_hubGene.csv (hub genes)", 5)),
  row.names = NULL
)
write.csv(sig_def, file.path(DIR$sig, "signature_definitions.csv"), row.names = FALSE)
write.csv(sig_def, file.path(DIR$tables, "signature_definitions.csv"), row.names = FALSE)
print(sig_def)

rank_symbols <- function(x) {
  y <- x |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "", is.finite(t)) |>
    dplyr::group_by(gene_symbol) |>
    dplyr::slice_max(abs(t), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  sort(setNames(y$t, y$gene_symbol), decreasing = TRUE)
}

sig_res <- lapply(names(deg), function(id) {
  fgsea::fgseaMultilevel(pathways = sig_list, stats = rank_symbols(deg[[id]]),
                         minSize = SIG_MIN, maxSize = SIG_MAX) |>
    as.data.frame() |>
    dplyr::mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";",
                                       FUN.VALUE = character(1)),
                  contrast_id = id)
}) |> dplyr::bind_rows()

sig_res <- sig_res |>
  dplyr::left_join(sig_def |> dplyr::select(pathway, orientation), by = "pathway") |>
  dplyr::mutate(reversal_score = NES * orientation,
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER))

write.csv(sig_res, file.path(DIR$sig, "Fig6F_signature_GSEA_all.csv"), row.names = FALSE)
write.csv(sig_res, file.path(DIR$tables, "Fig6F_signature_GSEA_all.csv"), row.names = FALSE)
message("05 signature GSEA done: ", length(sig_list), " signatures x ",
        dplyr::n_distinct(sig_res$contrast_id), " contrasts = ", nrow(sig_res), " rows")

# END ORIGINAL SCRIPT T06

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T07: 06_generate_figures.R
# ==============================================================================
## =============================================================================
## 06_generate_figures.R —— 出图（风格/配色/尺寸严格沿用 03）
## 沿用：Fig6A DEG计数 / Fig6B signed logFC / Fig6C1,C2 共识热图 /
##       Fig6D Hallmark NES / 补充 PCA、volcano、剂量一致性
## 修订：Fig6E 改为 7 个新 signature 逆转热图
## 新增：补充图 KEGG_LEGACY GSEA NES 热图；ORA 5 库 dotplot（术语完整）
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

V2 <- OUT; V2_MAIN <- DIR$main; V2_SUPP <- DIR$supp; V2_TAB <- DIR$tables
invisible(lapply(c(V2_MAIN, V2_SUPP, V2_TAB), dir.create, recursive = TRUE, showWarnings = FALSE))

deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))
dd <- dplyr::bind_rows(deg) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER))

label_bio_map <- lapply(c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot), function(p) {
  z <- data.table::fread(p, sep = "\t", select = c("Gene id", "Bio_type"),
                         data.table = FALSE, check.names = FALSE)
  names(z) <- c("ensembl_id", "bio_type"); z$ensembl_id <- sub("\\..*$", "", z$ensembl_id); z
}) |> dplyr::bind_rows() |>
  dplyr::arrange(ensembl_id, dplyr::desc(bio_type == "protein_coding")) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)

cell_cols <- c("MV4-11" = "#3B75AF", "THP-1" = "#4FAF88",
               "MOLM-13" = "#8C6BB1", "SKM-1" = "#D98C3F")
dose_cols <- c("TAM5" = "#E6A52C", "TAM12" = "#C9483B")

theme_v2 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "plain", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(3.3, "mm"), legend.key.width = grid::unit(3.3, "mm"),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}
save_v2 <- function(p, name, width_mm, height_mm, where = V2_MAIN, three_formats = TRUE) {
  base <- sub("\\.pdf$", "", name)
  ## PDF (矢量, 无损)
  ggplot2::ggsave(file.path(where, paste0(base, ".pdf")), p, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf, bg = "white")
  if (!three_formats) return(invisible(NULL))
  ## TIFF (300 dpi, 论文印刷标准)
  ggplot2::ggsave(file.path(where, paste0(base, ".tiff")), p, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, device = "tiff", bg = "white",
                  compression = "lzw")
  ## PNG (300 dpi, 演示/投稿备用)
  ggplot2::ggsave(file.path(where, paste0(base, ".png")), p, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, device = "png", bg = "white")
}

## ----------------------------- Fig6A -----------------------------
counts <- dd |>
  dplyr::filter(direction != "NS") |>
  dplyr::count(contrast_id, cell_line, dose, direction, name = "n") |>
  tidyr::complete(contrast_id, direction = c("Up", "Down"), fill = list(n = 0)) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                cell_line = factor(cell_line, levels = CELL_ORDER),
                dose = factor(dose, levels = c("TAM5", "TAM12")),
                direction = factor(direction, levels = c("Up", "Down")),
                signed_n = ifelse(direction == "Down", -n, n))
pA <- ggplot2::ggplot(counts, ggplot2::aes(dose, signed_n, fill = direction)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = .45, colour = "black") +
  ggplot2::geom_col(width = .66) +
  ggplot2::geom_text(ggplot2::aes(y = signed_n + ifelse(signed_n > 0, 25, -38),
                                  label = n, vjust = ifelse(signed_n > 0, 0, 1)),
                     colour = "black", size = 2.15, family = "Arial") +
  ggplot2::scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::scale_y_continuous(limits = c(-1150, 1500), breaks = seq(-1000, 1500, 500),
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::labs(x = NULL, y = "Number of DEGs") +
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(colour = "black"),
                 axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 axis.title.y = ggplot2::element_text(colour = "black"),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, face = "bold", margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.5, "mm"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(3, "mm"), legend.text = ggplot2::element_text(size = 7),
                 legend.box.just = "center", legend.justification = "center",
                 plot.margin = ggplot2::margin(3, 4, 2, 2, "mm"))
save_v2(pA, "main_A_DEG_counts.pdf", 130, 85)
write.csv(counts, file.path(V2_TAB, "Fig6A_DEG_counts.csv"), row.names = FALSE)

## ----------------------------- Fig6B -----------------------------
set.seed(20260720)
scatter <- dd |>
  dplyr::mutate(cell_line = factor(cell_line, levels = CELL_ORDER),
                dose = factor(dose, levels = c("TAM5", "TAM12")),
                x0 = as.numeric(dose), x = x0 + stats::runif(dplyr::n(), -.31, .31))
mapped_sig <- scatter |>
  dplyr::inner_join(label_bio_map, by = "ensembl_id") |>
  dplyr::filter(direction != "NS", !is.na(gene_symbol), gene_symbol != "",
                !grepl("^ENSG", gene_symbol), bio_type == "protein_coding")
lab_auto <- mapped_sig |> dplyr::group_by(contrast_id, direction) |>
  dplyr::slice_min(P.Value, n = 2, with_ties = FALSE) |> dplyr::ungroup()
lab_focus <- mapped_sig |> dplyr::filter(gene_symbol %in% FOCUS_GENES)
labs_v2 <- dplyr::bind_rows(lab_auto, lab_focus) |>
  dplyr::distinct(contrast_id, gene_symbol, .keep_all = TRUE)
ymax <- max(abs(scatter$logFC[is.finite(scatter$logFC)]), na.rm = TRUE)
## NS点下采样25%（可重复），大幅减少矢量对象数量（PDF从~4.5MB降至~1MB）；
## 显著点(Up/Down)全部保留，视觉效果几乎不变
set.seed(20260720)
ns_idx <- which(scatter$direction == "NS")
ns_keep <- sample(ns_idx, size = ceiling(length(ns_idx) * 0.25))
scatter_draw <- scatter[c(which(scatter$direction != "NS"), ns_keep), ]
pB <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = c(-DEG_LFC, DEG_LFC), linetype = 2, colour = "#666666", linewidth = .3) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = .35) +
  ggplot2::geom_point(data = scatter_draw, ggplot2::aes(x, logFC, colour = direction), size = .28, alpha = .42) +
  ggplot2::scale_colour_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#C9C9C9"),
                               breaks = c("Up", "Down", "NS")) +
  ggrepel::geom_text_repel(data = labs_v2, ggplot2::aes(x, logFC, label = gene_symbol),
                           size = 2.05, family = "Arial", max.overlaps = Inf, min.segment.length = 0,
                           seed = 20260720, box.padding = .18, point.padding = .08, segment.size = .22) +
  ggplot2::scale_x_continuous(breaks = 1:2, labels = c("5 μM", "12 μM"),
                              expand = ggplot2::expansion(mult = c(.015, .015))) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::coord_cartesian(ylim = c(-ymax, ymax), clip = "off") +
  ggplot2::labs(x = NULL, y = expression(log[2] * " fold change"), colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.5, "mm"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.text = ggplot2::element_text(size = 7),
                 legend.box.just = "center", legend.justification = "center",
                 plot.margin = ggplot2::margin(3, 4, 2, 2, "mm"))
save_v2(pB, "main_B_signed_logFC.pdf", 130, 85)

## ---------------------- main_A+B 并排拼图（已取消，不再生成） ----------------------
## pAB <- patchwork::wrap_plots(pA, pB, ncol = 2)
## save_v2(pAB, "main_AB_DEG_overview.pdf", 260, 92)

## ---------------------- Fig6C1/C2 consensus heatmaps ----------------------
annot_paths <- c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot)
bio_map <- lapply(annot_paths, function(p) {
  z <- data.table::fread(p, sep = "\t", select = c("Gene id", "Gene name", "Bio_type"),
                         data.table = FALSE, check.names = FALSE)
  names(z) <- c("ensembl_id", "gene_symbol_bio", "bio_type"); z$ensembl_id <- sub("\\..*$", "", z$ensembl_id); z
}) |> dplyr::bind_rows() |> dplyr::filter(ensembl_id != "") |>
  dplyr::arrange(ensembl_id, dplyr::desc(bio_type == "protein_coding")) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)
long <- dd |> dplyr::select(ensembl_id, gene_symbol, contrast_id, logFC, P.Value) |>
  dplyr::left_join(bio_map, by = "ensembl_id") |>
  dplyr::mutate(symbol = dplyr::coalesce(gene_symbol, gene_symbol_bio))
complete_pc <- long |>
  dplyr::filter(bio_type == "protein_coding", !is.na(symbol), symbol != "", !grepl("^ENSG", symbol)) |>
  dplyr::group_by(ensembl_id, symbol) |>
  dplyr::summarise(n_contrast = dplyr::n_distinct(contrast_id), .groups = "drop") |>
  dplyr::filter(n_contrast == 8)
rule_tbl <- read.csv(file.path(DIR$tables, "Fig6D_consensus_rule_3of4_and_2significant.csv"), check.names = FALSE)
meta_tbl <- read.csv(file.path(DIR$tables, "Fig6D_consensus_stouffer_meta.csv"), check.names = FALSE)
select_heat_ids <- function(tbl, n_per_direction = 8) {
  pool <- tbl |> dplyr::inner_join(complete_pc, by = "ensembl_id", suffix = c("", "_pc")) |>
    dplyr::mutate(response = ifelse(median_logFC >= 0, "Up", "Down")) |>
    dplyr::group_by(response, symbol) |>
    dplyr::slice_max(abs(median_logFC), n = 1, with_ties = FALSE) |> dplyr::ungroup() |>
    dplyr::arrange(response, dplyr::desc(abs(median_logFC)))
  dplyr::bind_rows(pool |> dplyr::filter(response == "Up") |> dplyr::slice_head(n = n_per_direction),
                   pool |> dplyr::filter(response == "Down") |> dplyr::slice_head(n = n_per_direction))
}
draw_consensus <- function(tbl, title_text, filename) {
  sel <- select_heat_ids(tbl, 8)
  if (!nrow(sel)) stop("No protein-coding consensus genes for ", title_text)
  hm <- long |> dplyr::filter(ensembl_id %in% sel$ensembl_id) |>
    dplyr::group_by(symbol, contrast_id) |>
    dplyr::summarise(logFC = mean(logFC), .groups = "drop") |>
    tidyr::pivot_wider(names_from = contrast_id, values_from = logFC)
  col_order <- c(paste(CELL_ORDER, "TAM5", sep = "_"), paste(CELL_ORDER, "TAM12", sep = "_"))
  mat <- as.matrix(hm[, col_order, drop = FALSE]); rownames(mat) <- make.unique(hm$symbol)
  mean_effect <- rowMeans(mat, na.rm = TRUE)
  selected_response <- sel$response[match(rownames(mat), sel$symbol)]
  row_split <- factor(ifelse(selected_response == "Up", "Upregulated", "Downregulated"),
                      levels = c("Upregulated", "Downregulated"))
  ord <- order(row_split, -abs(mean_effect)); mat <- mat[ord, , drop = FALSE]; row_split <- row_split[ord]
  dose_ann <- factor(rep(c("5 μM", "12 μM"), each = 4), levels = c("5 μM", "12 μM"))
  cell_ann <- rep(CELL_ORDER, 2)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Dose = dose_ann, Cell = cell_ann,
    col = list(Dose = c("5 μM" = dose_cols[["TAM5"]], "12 μM" = dose_cols[["TAM12"]]), Cell = cell_cols),
    annotation_name_gp = grid::gpar(fontsize = 7),
    annotation_legend_param = list(labels_gp = grid::gpar(fontsize = 7),
                                   title_gp = grid::gpar(fontsize = 7, fontface = "bold")))
  lim <- max(1, min(4, max(abs(mat), na.rm = TRUE)))
  ## 封装绘图表达式，在三种设备中分别执行
  draw_hm <- function() {
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
      mat, name = "log2FC",
      col = circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
      cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, column_split = dose_ann,
      top_annotation = ha, column_names_rot = 45, show_column_names = FALSE,
      row_names_gp = grid::gpar(fontsize = 7), row_title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
      column_title = NULL, border = TRUE,
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
                                  labels_gp = grid::gpar(fontsize = 7))),
      merge_legends = TRUE, align_heatmap_legend = "global_center",
      align_annotation_legend = "global_center", heatmap_legend_side = "right",
      annotation_legend_side = "right", padding = grid::unit(c(3, 3, 3, 3), "mm"))
  }
  base_fn <- sub("\\.pdf$", "", filename)
  ## PDF
  grDevices::cairo_pdf(file.path(V2_MAIN, paste0(base_fn, ".pdf")), width = 92 / 25.4, height = 82 / 25.4)
  draw_hm(); grDevices::dev.off()
  ## TIFF (300 dpi)
  grDevices::tiff(file.path(V2_MAIN, paste0(base_fn, ".tiff")), width = 92, height = 82,
                  units = "mm", res = 300, compression = "lzw")
  draw_hm(); grDevices::dev.off()
  ## PNG (300 dpi)
  grDevices::png(file.path(V2_MAIN, paste0(base_fn, ".png")), width = 92, height = 82,
                 units = "mm", res = 300)
  draw_hm(); grDevices::dev.off()
  write.csv(sel, file.path(V2_TAB, sub("pdf$", "csv", filename)), row.names = FALSE)
}
draw_consensus(rule_tbl, "Rule-based consensus", "main_C1_consensus_rule_heatmap.pdf")
draw_consensus(meta_tbl, "Stouffer integrated response", "main_C2_consensus_meta_heatmap.pdf")

## ----------------------------- Fig6D Hallmark (正文15条) -----------------------------
hall <- read.csv(file.path(DIR$tables, "GSEA_Hallmark_all.csv"), check.names = FALSE)
# 固定15条通路（4功能类别），全部至少在4/8比较中nominal P<0.05；非反向挑选
hall_15 <- data.frame(
  pathway = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_COMPLEMENT",
              "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  category = c(rep("Cell cycle", 4), rep("Metabolism", 4), rep("Immunity", 4), rep("Stress", 3)),
  stringsAsFactors = FALSE)
hall_category_order <- c("Cell cycle", "Metabolism", "Immunity", "Stress")
# 验证：每条至少4/8比较 nominal P<0.05
hall_check <- hall |> dplyr::filter(pathway %in% hall_15$pathway) |>
  dplyr::group_by(pathway) |> dplyr::summarise(n_p05 = sum(pval < .05, na.rm = TRUE), .groups = "drop")
if (any(hall_check$n_p05 < 4)) warning("Hallmark 15条中存在n_p05<4的通路: ",
  paste(hall_check$pathway[hall_check$n_p05 < 4], collapse = ", "))
clean_pathway <- function(x) {
  x |> stringr::str_remove("^HALLMARK_") |> stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_replace_all(c("Myc" = "MYC", "E2f" = "E2F", "G2m" = "G2M",
                               "Tnfa" = "TNFα", "Nfkb" = "NF-κB", "P53" = "p53",
                               "Mtorc1" = "mTORC1", " Via " = " via ")) |>
    stringr::str_wrap(width = 30)
}
contrast_lookup <- data.frame(contrast_id = CONTRAST_ORDER,
                              cell_line = factor(rep(CELL_ORDER, each = 2), levels = CELL_ORDER),
                              dose = factor(rep(c("TAM5", "TAM12"), 4), levels = c("TAM5", "TAM12")))
hall_15$ord <- seq_len(nrow(hall_15))  # 固定通路顺序（按类别内列表顺序）
hall_main <- hall |> dplyr::filter(pathway %in% hall_15$pathway) |>
  dplyr::left_join(hall_15, by = "pathway") |>
  dplyr::mutate(pathway_label = clean_pathway(pathway),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                category = factor(category, levels = hall_category_order),
                pathway_label = factor(pathway_label, levels = rev(clean_pathway(hall_15$pathway[order(hall_15$ord)]))),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
hall_plot <- hall_main |> dplyr::mutate(contrast_id = as.character(contrast_id)) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pD <- ggplot2::ggplot(hall_plot, ggplot2::aes(dose, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .45) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.4,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(category ~ cell_line, scales = "free", space = "free", switch = "both") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_text(size = 7), strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.background.y = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 strip.text.y.left = ggplot2::element_text(angle = 90, size = 7),
                 panel.spacing = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .35),
                 plot.margin = ggplot2::margin(2, 3, 2, 2, "mm"))
save_v2(pD, "main_D_Hallmark_NES.pdf", 135, 86)
write.csv(hall_check, file.path(V2_TAB, "Hallmark_15pathways_n_p05_check.csv"), row.names = FALSE)
write.csv(hall_main, file.path(V2_TAB, "Fig6D_Hallmark_15pathways_data.csv"), row.names = FALSE)

## Full Hallmark heatmap (supp)
hall_rank <- hall |> dplyr::group_by(pathway) |>
  dplyr::summarise(mean_abs_NES = mean(abs(NES), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_abs_NES), pathway)
hall_full <- hall |>
  dplyr::mutate(pathway_label = clean_pathway(pathway),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                pathway_label = factor(pathway_label, levels = rev(clean_pathway(hall_rank$pathway))))
pDfull <- ggplot2::ggplot(hall_full, ggplot2::aes(contrast_id, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .2) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = rep(c("5 μM", "12 μM"), 4)) +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 5.2), legend.position = "right")
save_v2(pDfull, "FigS_full_Hallmark_NES_heatmap.pdf", 155, 190, where = V2_SUPP, three_formats = FALSE)

## ------------------- Fig6D-KEGG: KEGG 正文热图（同样式） -------------------
gk <- read.csv(file.path(DIR$tables, "GSEA_KEGG_LEGACY_all.csv"), check.names = FALSE)
kegg_main_list <- data.frame(
  pathway = c("KEGG_CELL_CYCLE", "KEGG_DNA_REPLICATION", "KEGG_OOCYTE_MEIOSIS", "KEGG_SPLICEOSOME",
              "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
              "KEGG_PENTOSE_PHOSPHATE_PATHWAY", "KEGG_PURINE_METABOLISM",
              "KEGG_RIBOSOME", "KEGG_PROTEASOME",
              "KEGG_LYSOSOME", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
              "KEGG_CHEMOKINE_SIGNALING_PATHWAY"),
  label = c("Cell cycle", "DNA replication", "Oocyte meiosis", "Spliceosome",
            "Oxidative phosphorylation", "Glycolysis / gluconeogenesis",
            "Pentose phosphate pathway", "Purine metabolism",
            "Ribosome", "Proteasome",
            "Lysosome", "Cytokine-receptor interaction",
            "Chemokine signaling pathway"),
  category = c(rep("Cell cycle", 4), rep("Metabolism", 4),
               rep("Ribosome", 2), rep("Stress", 3)),
  stringsAsFactors = FALSE)
kegg_category_order <- c("Cell cycle", "Metabolism", "Ribosome", "Stress")
kegg_check <- gk |> dplyr::filter(pathway %in% kegg_main_list$pathway) |>
  dplyr::group_by(pathway) |> dplyr::summarise(n_p05 = sum(pval < .05, na.rm = TRUE), .groups = "drop")
if (any(kegg_check$n_p05 < 4)) warning("KEGG正文中存在n_p05<4的通路: ",
  paste(kegg_check$pathway[kegg_check$n_p05 < 4], collapse = ", "))
kegg_main_list$ord <- seq_len(nrow(kegg_main_list))
kegg_main <- gk |> dplyr::filter(pathway %in% kegg_main_list$pathway) |>
  dplyr::left_join(kegg_main_list, by = "pathway") |>
  dplyr::mutate(pathway_label = factor(label, levels = rev(kegg_main_list$label[order(kegg_main_list$ord)])),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                category = factor(category, levels = kegg_category_order),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
kegg_plot <- kegg_main |> dplyr::mutate(contrast_id = as.character(contrast_id)) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pDK <- ggplot2::ggplot(kegg_plot, ggplot2::aes(dose, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .45) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.4,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(category ~ cell_line, scales = "free", space = "free", switch = "both") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_text(size = 7), strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.background.y = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 strip.text.y.left = ggplot2::element_text(angle = 90, size = 6.5),
                 panel.spacing = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .35),
                 plot.margin = ggplot2::margin(2, 3, 2, 5, "mm"))
save_v2(pDK, "main_D_KEGG_NES.pdf", 140, 82)
write.csv(kegg_check, file.path(V2_TAB, "KEGG_main_pathways_n_p05_check.csv"), row.names = FALSE)
write.csv(kegg_main, file.path(V2_TAB, "Fig6D_KEGG_main_pathways_data.csv"), row.names = FALSE)

## ------------------- Fig6E: HSC-like signature 逆转（仅保留1个） -------------------
sig_res <- read.csv(file.path(DIR$tables, "Fig6F_signature_GSEA_all.csv"), check.names = FALSE)
sig_map <- data.frame(
  pathway = "HSC_LIKE_MARKERS",
  label = "HSC-like markers",
  stringsAsFactors = FALSE)
sig_main <- sig_res |> dplyr::inner_join(sig_map, by = "pathway") |>
  dplyr::mutate(reversal_score = NES * orientation,
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                label = factor(label, levels = rev(sig_map$label)),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
sig_plot <- sig_main |> dplyr::mutate(contrast_id = as.character(contrast_id),
                                      label = factor(label, levels = rev(sig_map$label))) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pE <- ggplot2::ggplot(sig_plot, ggplot2::aes(dose, label, fill = reversal_score)) +
  ggplot2::geom_tile(colour = "white", linewidth = .7) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.8,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, name = "Reversal score") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .4),
                 plot.margin = ggplot2::margin(2, 3, 2, 2, "mm"))
save_v2(pE, "main_E_HSC_reversal.pdf", 125, 28)
write.csv(sig_main, file.path(V2_TAB, "Fig6E_HSC_like_reversal_scores.csv"), row.names = FALSE)

## ----------------- Supplementary: KEGG_LEGACY GSEA heatmap -----------------
gk <- read.csv(file.path(DIR$tables, "GSEA_KEGG_LEGACY_all.csv"), check.names = FALSE)
gk_rank <- gk |> dplyr::group_by(pathway) |>
  dplyr::summarise(n_fdr = sum(padj < .05, na.rm = TRUE),
                   meanNES = mean(NES, na.rm = TRUE),
                   mean_abs = mean(abs(NES), na.rm = TRUE), .groups = "drop")
gk_sel <- gk_rank |> dplyr::filter(n_fdr >= 2) |>
  dplyr::arrange(dplyr::desc(n_fdr), dplyr::desc(mean_abs)) |>
  dplyr::slice_head(n = 60)
clean_kegg <- function(x) x |> stringr::str_remove("^KEGG_") |> stringr::str_replace_all("_", " ")
if (nrow(gk_sel)) {
  gk_df <- gk |> dplyr::filter(pathway %in% gk_sel$pathway) |>
    dplyr::mutate(plabel = clean_kegg(pathway),
                  plabel = factor(plabel, levels = rev(clean_kegg(
                    gk_sel |> dplyr::arrange(meanNES) |> dplyr::pull(pathway)))),
                  contrast_id = as.character(contrast_id)) |>
    dplyr::left_join(contrast_lookup, by = "contrast_id")
  pK <- ggplot2::ggplot(gk_df, ggplot2::aes(dose, plabel, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = .2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
    ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
    ggplot2::facet_grid(~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
    ggplot2::labs(x = NULL, y = NULL) + theme_v2(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7),
                   axis.text.y = ggplot2::element_text(size = 5.6),
                   strip.placement = "outside",
                   strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                   strip.text.x = ggplot2::element_text(size = 8),
                   panel.spacing.x = grid::unit(1, "mm"),
                   legend.position = "right")
  h_k <- max(60, nrow(gk_sel) * 4.0)
  save_v2(pK, "FigS_KEGG_LEGACY_GSEA_NES_heatmap.pdf", 170, h_k, where = V2_SUPP, three_formats = FALSE)
  write.csv(gk_sel, file.path(V2_TAB, "KEGG_LEGACY_GSEA_selected_pathways.csv"), row.names = FALSE)
}

## ----------------------- Supplementary: ORA dotplots -----------------------
ora_top <- read.csv(file.path(DIR$tables, "ORA_sig_p0.05_top10.csv"), check.names = FALSE)
ora_top$contrast_id <- factor(ora_top$contrast_id, levels = CONTRAST_ORDER)
ora_top$neglog10p <- -log10(pmax(ora_top$pvalue, 1e-300))
clean_ora_term <- function(x, coll) {
  x <- if (coll == "KEGG") sub("^KEGG_", "", x)
       else if (coll == "Hallmark") sub("^HALLMARK_", "", x) else x
  stringr::str_wrap(gsub("_", " ", x), width = 40)
}
ora_panel_levels <- paste(rep(CELL_ORDER, each = 2), rep(c("5 μM", "12 μM"), 4))
for (coll in c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Hallmark")) {
  df <- ora_top |> dplyr::filter(collection == coll)
  if (!nrow(df)) next
  # ORA 结果表本身已含 cell_line/dose；构造短面板标题（细胞系 + 剂量，固定顺序）
  df <- df |>
    dplyr::mutate(panel_lab = paste(cell_line, ifelse(dose == "TAM5", "5 μM", "12 μM")),
                  panel_lab = factor(panel_lab, levels = ora_panel_levels),
                  term_w = clean_ora_term(Description, coll))
  lv <- df |> dplyr::group_by(term_w) |>
    dplyr::summarise(mp = min(pvalue), .groups = "drop") |>
    dplyr::arrange(mp) |> dplyr::pull(term_w) |> unique()
  df$term_w <- factor(df$term_w, levels = rev(lv))
  p <- ggplot2::ggplot(df, ggplot2::aes(GeneRatioNum, term_w, colour = neglog10p, size = Count)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#F4A582", high = "#B2182B", name = expression(-log[10] * italic(P))) +
    ggplot2::scale_size_continuous(range = c(1.0, 3.2), name = "Count") +
    ggplot2::scale_x_continuous(n.breaks = 4, labels = scales::percent_format(accuracy = 1)) +
    ggplot2::facet_wrap(~ panel_lab, ncol = 2, scales = "free") +
    ggplot2::labs(x = "Gene ratio", y = NULL, title = NULL) + theme_v2(7) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.8),
                   axis.text.x = ggplot2::element_text(size = 6),
                   strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   panel.border = ggplot2::element_rect(fill = NA, linewidth = .3),
                   panel.spacing = grid::unit(3, "mm"),
                   legend.position = "right")
  hgt <- nrow(df) * 4.6 / 2 + 24
  save_v2(p, paste0("FigS_ORA_", coll, "_dotplot.pdf"), 180, max(120, hgt), where = V2_SUPP, three_formats = FALSE)
}

## --------------------------- Supplementary PCA ---------------------------
pca_list <- list()
for (id in CONTRAST_ORDER) {
  bits <- strsplit(id, "_", fixed = TRUE)[[1]]; dose_i <- tail(bits, 1)
  cell_i <- sub(paste0("_", dose_i, "$"), "", id)
  obj <- readRDS(file.path(DIR$deg, paste(cell_i, dose_i, sep = "_"), "limma_voom_object.rds"))
  e <- obj$voom$E
  vv <- matrixStats::rowVars(e)
  use <- order(vv, decreasing = TRUE)[seq_len(min(500, length(vv)))]
  pc <- stats::prcomp(t(e[use, , drop = FALSE]), center = TRUE, scale. = FALSE)
  d <- cbind(obj$metadata, as.data.frame(pc$x[, 1:2, drop = FALSE]))
  d$group <- factor(as.character(d$dose), levels = c("Control", dose_i))
  ve <- summary(pc)$importance[2, 1:2] * 100
  pca_list[[id]] <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = group)) +
    ggplot2::geom_point(size = 2.2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = replicate), size = 2.2, show.legend = FALSE, min.segment.length = 0) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#555555", dose_cols[[dose_i]]), c("Control", dose_i))) +
    ggplot2::labs(title = id, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_v2(7) + ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"), legend.position = "top")
}
pca_panel <- patchwork::wrap_plots(pca_list, ncol = 4, guides = "collect") & ggplot2::theme(legend.position = "top")
save_v2(pca_panel, "FigS_eight_matched_PCA.pdf", 180, 105, where = V2_SUPP, three_formats = FALSE)

## ------------------------- Supplementary volcano -------------------------
vol_list <- lapply(CONTRAST_ORDER, function(id) {
  z <- deg[[id]]
  la <- z |> dplyr::filter(direction != "NS", !is.na(gene_symbol), gene_symbol != "", !grepl("^ENSG", gene_symbol)) |>
    dplyr::group_by(direction) |> dplyr::slice_min(P.Value, n = 2, with_ties = FALSE) |> dplyr::ungroup()
  lf <- z |> dplyr::filter(direction != "NS", gene_symbol %in% FOCUS_GENES)
  lb <- dplyr::bind_rows(la, lf) |> dplyr::distinct(gene_symbol, .keep_all = TRUE)
  ggplot2::ggplot(z, ggplot2::aes(logFC, -log10(pmax(P.Value, 1e-300)), colour = direction)) +
    ggplot2::geom_point(size = .32, alpha = .5) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = .3) +
    ggplot2::geom_hline(yintercept = -log10(.05), linetype = 2, linewidth = .3) +
    ggrepel::geom_text_repel(data = lb, ggplot2::aes(label = gene_symbol), size = 1.8,
                             min.segment.length = 0, max.overlaps = Inf) +
    ggplot2::scale_colour_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#C7C7C7")) +
    ggplot2::labs(title = id, x = expression(log[2] * " fold change"),
                  y = expression(-log[10] * " nominal " * italic(P))) + theme_v2(7) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"), legend.position = "none")
})
vol_panel <- patchwork::wrap_plots(vol_list, ncol = 4)
save_v2(vol_panel, "FigS_eight_volcano_plots.pdf", 180, 115, where = V2_SUPP, three_formats = FALSE)

## ------------------- Supplementary dose concordance -------------------
conc <- lapply(CELL_ORDER, function(cell_i) {
  a <- deg[[paste(cell_i, "TAM5", sep = "_")]] |>
    dplyr::select(ensembl_id, gene_symbol, lfc5 = logFC, p5 = P.Value)
  b <- deg[[paste(cell_i, "TAM12", sep = "_")]] |>
    dplyr::select(ensembl_id, lfc12 = logFC, p12 = P.Value)
  dplyr::inner_join(a, b, by = "ensembl_id") |>
    dplyr::mutate(cell_line = cell_i,
                  concordant_sig = abs(lfc5) > 1 & p5 < .05 & abs(lfc12) > 1 & p12 < .05 & sign(lfc5) == sign(lfc12))
}) |> dplyr::bind_rows()
rho <- conc |> dplyr::group_by(cell_line) |>
  dplyr::summarise(rho = cor(lfc5, lfc12, method = "spearman", use = "complete.obs"), .groups = "drop")
conc <- dplyr::left_join(conc, rho, by = "cell_line")
pconc <- ggplot2::ggplot(conc, ggplot2::aes(lfc5, lfc12)) +
  ggplot2::geom_bin_2d(bins = 60) +
  ggplot2::geom_point(data = dplyr::filter(conc, concordant_sig), shape = 1, size = .45, colour = "#B2182B", alpha = .65) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .35) +
  ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = .25) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = .25) +
  ggplot2::facet_wrap(~ cell_line, scales = "free") +
  ggplot2::geom_text(data = dplyr::distinct(conc, cell_line, rho),
                     ggplot2::aes(x = -Inf, y = Inf, label = sprintf("Spearman ρ = %.2f", rho)),
                     inherit.aes = FALSE, hjust = -.08, vjust = 1.3, size = 2.4) +
  ggplot2::scale_fill_gradient(low = "#F1F5F9", high = "#2166AC", name = "Gene density") +
  ggplot2::labs(x = expression("TAM 5 μM log"[2] * "FC vs matched control"),
                y = expression("TAM 12 μM log"[2] * "FC vs matched control")) +
  theme_v2(8) + ggplot2::theme(legend.position = "right")
save_v2(pconc, "FigS_TAM5_TAM12_effect_concordance.pdf", 150, 120, where = V2_SUPP, three_formats = FALSE)
write.csv(rho, file.path(V2_TAB, "TAM5_TAM12_effect_spearman.csv"), row.names = FALSE)

writeLines(c(
  "06_analysis_revised figure generation completed (v2: Hallmark15 + KEGG main + HSC-only).",
  "01/02 rebuilt from raw counts; DEG threshold |log2FC|>1 & nominal P<0.05 (same as 03).",
  "ORA: merged up+down DEGs; universe = all tested genes; significance raw pvalue<0.05; top10 shown; GO BP/CC/MF + KEGG(MSigDB KEGG_LEGACY) + Hallmark.",
  "GSEA: Hallmark (H, 50 sets) + KEGG (C2/CP:KEGG_LEGACY, 186 sets); fgseaMultilevel minSize=10 maxSize=500; ranking statistic = moderated t.",
  "Fig6D Hallmark: 15 fixed pathways in 4 categories (Cell cycle/Metabolism/Immunity/Stress); Fig6D-KEGG: 13 main pathways in 4 categories.",
  "Fig6E: HSC_LIKE_MARKERS only (orientation=-1, red=suppressed/reversed).",
  "Asterisks in Fig6D/Fig6E indicate nominal P (pval): * <0.05, ** <0.01, *** <0.001, **** <0.0001; full FDR in supplementary tables."
), file.path(V2, "RUN_SUMMARY.txt"))
message("06 figures done. Main: ", V2_MAIN, " | Supp: ", V2_SUPP)

# END ORIGINAL SCRIPT T07

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T08: rerun_main_B_only.R
# ==============================================================================
## =============================================================================
## rerun_main_B_only.R —— 仅重绘 main_B_signed_logFC（NS点下采样25%优化版）
## 改动：NS点可重复下采样至25%，显著点全部保留，PDF体积从~4.5MB降至~1MB
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

V2_MAIN <- DIR$main

## ---- 读取DEG数据（与06_generate_figures.R一致） ----
deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))
dd <- dplyr::bind_rows(deg) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER))

label_bio_map <- lapply(c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot), function(p) {
  z <- data.table::fread(p, sep = "\t", select = c("Gene id", "Bio_type"),
                         data.table = FALSE, check.names = FALSE)
  names(z) <- c("ensembl_id", "bio_type"); z$ensembl_id <- sub("\\..*$", "", z$ensembl_id); z
}) |> dplyr::bind_rows() |>
  dplyr::arrange(ensembl_id, dplyr::desc(bio_type == "protein_coding")) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)

theme_v2 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "plain", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(3.3, "mm"), legend.key.width = grid::unit(3.3, "mm"),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}

save_v2 <- function(p, name, width_mm, height_mm, where = V2_MAIN, three_formats = TRUE) {
  base <- sub("\\.pdf$", "", name)
  ggplot2::ggsave(file.path(where, paste0(base, ".pdf")), p, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf, bg = "white")
  if (!three_formats) return(invisible(NULL))
  ggplot2::ggsave(file.path(where, paste0(base, ".tiff")), p, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, device = "tiff", bg = "white", compression = "lzw")
  ggplot2::ggsave(file.path(where, paste0(base, ".png")), p, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, device = "png", bg = "white")
}

## ----------------------------- main_B -----------------------------
set.seed(20260720)
scatter <- dd |>
  dplyr::mutate(cell_line = factor(cell_line, levels = CELL_ORDER),
                dose = factor(dose, levels = c("TAM5", "TAM12")),
                x0 = as.numeric(dose), x = x0 + stats::runif(dplyr::n(), -.31, .31))
mapped_sig <- scatter |>
  dplyr::inner_join(label_bio_map, by = "ensembl_id") |>
  dplyr::filter(direction != "NS", !is.na(gene_symbol), gene_symbol != "",
                !grepl("^ENSG", gene_symbol), bio_type == "protein_coding")
lab_auto <- mapped_sig |> dplyr::group_by(contrast_id, direction) |>
  dplyr::slice_min(P.Value, n = 2, with_ties = FALSE) |> dplyr::ungroup()
lab_focus <- mapped_sig |> dplyr::filter(gene_symbol %in% FOCUS_GENES)
labs_v2 <- dplyr::bind_rows(lab_auto, lab_focus) |>
  dplyr::distinct(contrast_id, gene_symbol, .keep_all = TRUE)
ymax <- max(abs(scatter$logFC[is.finite(scatter$logFC)]), na.rm = TRUE)

## NS点下采样25%（可重复），显著点全部保留
set.seed(20260720)
ns_idx <- which(scatter$direction == "NS")
ns_keep <- sample(ns_idx, size = ceiling(length(ns_idx) * 0.25))
scatter_draw <- scatter[c(which(scatter$direction != "NS"), ns_keep), ]

cat(sprintf("原始点数: %d | 下采样后点数: %d | 减少: %.1f%%\n",
            nrow(scatter), nrow(scatter_draw),
            (1 - nrow(scatter_draw) / nrow(scatter)) * 100))

pB <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = c(-DEG_LFC, DEG_LFC), linetype = 2, colour = "#666666", linewidth = .3) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = .35) +
  ggplot2::geom_point(data = scatter_draw, ggplot2::aes(x, logFC, colour = direction), size = .28, alpha = .42) +
  ggplot2::scale_colour_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#C9C9C9"),
                               breaks = c("Up", "Down", "NS")) +
  ggrepel::geom_text_repel(data = labs_v2, ggplot2::aes(x, logFC, label = gene_symbol),
                           size = 2.05, family = "Arial", max.overlaps = Inf, min.segment.length = 0,
                           seed = 20260720, box.padding = .18, point.padding = .08, segment.size = .22) +
  ggplot2::scale_x_continuous(breaks = 1:2, labels = c("5 μM", "12 μM"),
                              expand = ggplot2::expansion(mult = c(.015, .015))) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::coord_cartesian(ylim = c(-ymax, ymax), clip = "off") +
  ggplot2::labs(x = NULL, y = expression(log[2] * " fold change"), colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.5, "mm"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.text = ggplot2::element_text(size = 7),
                 legend.box.just = "center", legend.justification = "center",
                 plot.margin = ggplot2::margin(3, 4, 2, 2, "mm"))
save_v2(pB, "main_B_signed_logFC.pdf", 130, 85)

cat("main_B 重绘完成（NS下采样25%版）。输出目录:", V2_MAIN, "\n")

# END ORIGINAL SCRIPT T08

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T09: rerun_main_E_only.R
# ==============================================================================
## 仅重绘 main_E_HSC_reversal.pdf（修复右侧颜色条标题被截断）
## 严格照搬 06_generate_figures.R 中 main_E 代码，仅改 top margin 和图高度

Sys.setenv(TAM_ANALYSIS_ROOT = normalizePath(getwd(), winslash = "/"))
source("R/00_config.R")

V2_MAIN <- DIR$main

theme_v2 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "plain", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(3.3, "mm"), legend.key.width = grid::unit(3.3, "mm"),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}
save_v2 <- function(p, name, width_mm, height_mm, where = V2_MAIN) {
  ggplot2::ggsave(file.path(where, name), p, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf, bg = "white")
}

## contrast_lookup（与 06 脚本第 225 行一致）
contrast_lookup <- data.frame(contrast_id = CONTRAST_ORDER,
                              cell_line = factor(rep(CELL_ORDER, each = 2), levels = CELL_ORDER),
                              dose = factor(rep(c("TAM5", "TAM12"), 4), levels = c("TAM5", "TAM12")))

## 以下完全照搬 06 脚本第 333-367 行，仅 plot.margin top 2→5，height 28→32
sig_res <- read.csv(file.path(DIR$tables, "Fig6F_signature_GSEA_all.csv"), check.names = FALSE)
sig_map <- data.frame(
  pathway = "HSC_LIKE_MARKERS",
  label = "HSC-like markers",
  stringsAsFactors = FALSE)
sig_main <- sig_res |> dplyr::inner_join(sig_map, by = "pathway") |>
  dplyr::mutate(reversal_score = NES * orientation,
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                label = factor(label, levels = rev(sig_map$label)),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
sig_plot <- sig_main |> dplyr::mutate(contrast_id = as.character(contrast_id),
                                      label = factor(label, levels = rev(sig_map$label))) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")

pE <- ggplot2::ggplot(sig_plot, ggplot2::aes(dose, label, fill = reversal_score)) +
  ggplot2::geom_tile(colour = "white", linewidth = .7) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.8,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, name = "Reversal score") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 \u03bcM", "TAM12" = "12 \u03bcM")) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                 axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.title = ggplot2::element_text(colour = "black"),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.background = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 8, face = "bold", margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.0, "mm"),
                 legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8, face = "plain"),
                 legend.text = ggplot2::element_text(size = 7),
                 legend.key.height = grid::unit(3.3, "mm"),
                 legend.key.width = grid::unit(3.3, "mm"),
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .4),
                 plot.margin = ggplot2::margin(5, 3, 2, 2, "mm"))
save_v2(pE, "main_E_HSC_reversal.pdf", 140, 32)
cat("main_E_HSC_reversal.pdf re-drawn: top margin 2→5mm, height 28→32mm.\n")

# END ORIGINAL SCRIPT T09

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T10: run_all_revised.R
# ==============================================================================
## =============================================================================
## run_all_revised.R —— 02_analysis 全流程总控（从原始 count 完整重跑）
## 用法（PowerShell，在 02_analysis 目录下）：
##   & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" R\run_all_revised.R
## 顺序：00 配置 -> 01 数据准备/DEG -> 02 共识 -> 03 ORA -> 04 GSEA -> 05 signature -> 06 出图
## =============================================================================
root <- normalizePath(getwd(), winslash = "/")
Sys.setenv(TAM_ANALYSIS_ROOT = root)
t0 <- Sys.time()
scripts <- c("00_config.R", "01_prepare_and_DEG.R", "02_consensus_analysis.R",
             "03_ORA_GO_KEGG_Hallmark.R", "04_GSEA_Hallmark_KEGG.R",
             "05_signature_GSEA_revised.R", "06_generate_figures.R")
for (f in scripts) {
  message("\n==================== RUN ", f, " ====================")
  source(file.path(root, "R", f), echo = FALSE, max.deparse.length = 0)
}
writeLines(capture.output(sessionInfo()), file.path(root, "09_logs", "sessionInfo_run_all.txt"))
message("\nALL DONE in ", format(round(Sys.time() - t0, 2)))

# END ORIGINAL SCRIPT T10

# ==============================================================================
# BEGIN ORIGINAL SCRIPT T11: 06_generate_figures_backup_before_3formats.R | HISTORICAL BACKUP; not the current figure script
# ==============================================================================
## =============================================================================
## 06_generate_figures.R —— 出图（风格/配色/尺寸严格沿用 03）
## 沿用：Fig6A DEG计数 / Fig6B signed logFC / Fig6C1,C2 共识热图 /
##       Fig6D Hallmark NES / 补充 PCA、volcano、剂量一致性
## 修订：Fig6E 改为 7 个新 signature 逆转热图
## 新增：补充图 KEGG_LEGACY GSEA NES 热图；ORA 5 库 dotplot（术语完整）
## =============================================================================
source(file.path(Sys.getenv("TAM_ANALYSIS_ROOT", unset = getwd()), "R", "00_config.R"))

V2 <- OUT; V2_MAIN <- DIR$main; V2_SUPP <- DIR$supp; V2_TAB <- DIR$tables
invisible(lapply(c(V2_MAIN, V2_SUPP, V2_TAB), dir.create, recursive = TRUE, showWarnings = FALSE))

deg <- readRDS(file.path(DIR$intermediate, "DEG_all_8_contrasts.rds"))
dd <- dplyr::bind_rows(deg) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER))

label_bio_map <- lapply(c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot), function(p) {
  z <- data.table::fread(p, sep = "\t", select = c("Gene id", "Bio_type"),
                         data.table = FALSE, check.names = FALSE)
  names(z) <- c("ensembl_id", "bio_type"); z$ensembl_id <- sub("\\..*$", "", z$ensembl_id); z
}) |> dplyr::bind_rows() |>
  dplyr::arrange(ensembl_id, dplyr::desc(bio_type == "protein_coding")) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)

cell_cols <- c("MV4-11" = "#3B75AF", "THP-1" = "#4FAF88",
               "MOLM-13" = "#8C6BB1", "SKM-1" = "#D98C3F")
dose_cols <- c("TAM5" = "#E6A52C", "TAM12" = "#C9483B")

theme_v2 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "plain", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(3.3, "mm"), legend.key.width = grid::unit(3.3, "mm"),
      plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}
save_v2 <- function(p, name, width_mm, height_mm, where = V2_MAIN) {
  ggplot2::ggsave(file.path(where, name), p, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf, bg = "white")
}

## ----------------------------- Fig6A -----------------------------
counts <- dd |>
  dplyr::filter(direction != "NS") |>
  dplyr::count(contrast_id, cell_line, dose, direction, name = "n") |>
  tidyr::complete(contrast_id, direction = c("Up", "Down"), fill = list(n = 0)) |>
  dplyr::mutate(contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                cell_line = factor(cell_line, levels = CELL_ORDER),
                dose = factor(dose, levels = c("TAM5", "TAM12")),
                direction = factor(direction, levels = c("Up", "Down")),
                signed_n = ifelse(direction == "Down", -n, n))
pA <- ggplot2::ggplot(counts, ggplot2::aes(dose, signed_n, fill = direction)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = .45, colour = "black") +
  ggplot2::geom_col(width = .66) +
  ggplot2::geom_text(ggplot2::aes(y = signed_n + ifelse(signed_n > 0, 25, -38),
                                  label = n, vjust = ifelse(signed_n > 0, 0, 1)),
                     colour = "black", size = 2.15, family = "Arial") +
  ggplot2::scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::scale_y_continuous(limits = c(-1150, 1500), breaks = seq(-1000, 1500, 500),
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::labs(x = NULL, y = "Number of DEGs") +
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(colour = "black"),
                 axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 axis.title.y = ggplot2::element_text(colour = "black"),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, face = "bold", margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.5, "mm"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(3, "mm"), legend.text = ggplot2::element_text(size = 7),
                 legend.box.just = "center", legend.justification = "center",
                 plot.margin = ggplot2::margin(3, 4, 2, 2, "mm"))
save_v2(pA, "main_A_DEG_counts.pdf", 130, 85)
write.csv(counts, file.path(V2_TAB, "Fig6A_DEG_counts.csv"), row.names = FALSE)

## ----------------------------- Fig6B -----------------------------
set.seed(20260720)
scatter <- dd |>
  dplyr::mutate(cell_line = factor(cell_line, levels = CELL_ORDER),
                dose = factor(dose, levels = c("TAM5", "TAM12")),
                x0 = as.numeric(dose), x = x0 + stats::runif(dplyr::n(), -.31, .31))
mapped_sig <- scatter |>
  dplyr::inner_join(label_bio_map, by = "ensembl_id") |>
  dplyr::filter(direction != "NS", !is.na(gene_symbol), gene_symbol != "",
                !grepl("^ENSG", gene_symbol), bio_type == "protein_coding")
lab_auto <- mapped_sig |> dplyr::group_by(contrast_id, direction) |>
  dplyr::slice_min(P.Value, n = 2, with_ties = FALSE) |> dplyr::ungroup()
lab_focus <- mapped_sig |> dplyr::filter(gene_symbol %in% FOCUS_GENES)
labs_v2 <- dplyr::bind_rows(lab_auto, lab_focus) |>
  dplyr::distinct(contrast_id, gene_symbol, .keep_all = TRUE)
ymax <- max(abs(scatter$logFC[is.finite(scatter$logFC)]), na.rm = TRUE)
pB <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = c(-DEG_LFC, DEG_LFC), linetype = 2, colour = "#666666", linewidth = .3) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = .35) +
  ggplot2::geom_point(data = scatter, ggplot2::aes(x, logFC, colour = direction), size = .28, alpha = .42) +
  ggplot2::scale_colour_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#C9C9C9"),
                               breaks = c("Up", "Down", "NS")) +
  ggrepel::geom_text_repel(data = labs_v2, ggplot2::aes(x, logFC, label = gene_symbol),
                           size = 2.05, family = "Arial", max.overlaps = Inf, min.segment.length = 0,
                           seed = 20260720, box.padding = .18, point.padding = .08, segment.size = .22) +
  ggplot2::scale_x_continuous(breaks = 1:2, labels = c("5 μM", "12 μM"),
                              expand = ggplot2::expansion(mult = c(.015, .015))) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::coord_cartesian(ylim = c(-ymax, ymax), clip = "off") +
  ggplot2::labs(x = NULL, y = expression(log[2] * " fold change"), colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.5, "mm"),
                 legend.position = "top", legend.direction = "horizontal",
                 legend.text = ggplot2::element_text(size = 7),
                 legend.box.just = "center", legend.justification = "center",
                 plot.margin = ggplot2::margin(3, 4, 2, 2, "mm"))
save_v2(pB, "main_B_signed_logFC.pdf", 130, 85)

## ----------------------------- main_A+B 并排拼图 -----------------------------
pAB <- patchwork::wrap_plots(pA, pB, ncol = 2)
save_v2(pAB, "main_AB_DEG_overview.pdf", 260, 92)

## ---------------------- Fig6C1/C2 consensus heatmaps ----------------------
annot_paths <- c(PATH$b1_annot, PATH$b2_skm_annot, PATH$b2_other_annot)
bio_map <- lapply(annot_paths, function(p) {
  z <- data.table::fread(p, sep = "\t", select = c("Gene id", "Gene name", "Bio_type"),
                         data.table = FALSE, check.names = FALSE)
  names(z) <- c("ensembl_id", "gene_symbol_bio", "bio_type"); z$ensembl_id <- sub("\\..*$", "", z$ensembl_id); z
}) |> dplyr::bind_rows() |> dplyr::filter(ensembl_id != "") |>
  dplyr::arrange(ensembl_id, dplyr::desc(bio_type == "protein_coding")) |>
  dplyr::distinct(ensembl_id, .keep_all = TRUE)
long <- dd |> dplyr::select(ensembl_id, gene_symbol, contrast_id, logFC, P.Value) |>
  dplyr::left_join(bio_map, by = "ensembl_id") |>
  dplyr::mutate(symbol = dplyr::coalesce(gene_symbol, gene_symbol_bio))
complete_pc <- long |>
  dplyr::filter(bio_type == "protein_coding", !is.na(symbol), symbol != "", !grepl("^ENSG", symbol)) |>
  dplyr::group_by(ensembl_id, symbol) |>
  dplyr::summarise(n_contrast = dplyr::n_distinct(contrast_id), .groups = "drop") |>
  dplyr::filter(n_contrast == 8)
rule_tbl <- read.csv(file.path(DIR$tables, "Fig6D_consensus_rule_3of4_and_2significant.csv"), check.names = FALSE)
meta_tbl <- read.csv(file.path(DIR$tables, "Fig6D_consensus_stouffer_meta.csv"), check.names = FALSE)
select_heat_ids <- function(tbl, n_per_direction = 8) {
  pool <- tbl |> dplyr::inner_join(complete_pc, by = "ensembl_id", suffix = c("", "_pc")) |>
    dplyr::mutate(response = ifelse(median_logFC >= 0, "Up", "Down")) |>
    dplyr::group_by(response, symbol) |>
    dplyr::slice_max(abs(median_logFC), n = 1, with_ties = FALSE) |> dplyr::ungroup() |>
    dplyr::arrange(response, dplyr::desc(abs(median_logFC)))
  dplyr::bind_rows(pool |> dplyr::filter(response == "Up") |> dplyr::slice_head(n = n_per_direction),
                   pool |> dplyr::filter(response == "Down") |> dplyr::slice_head(n = n_per_direction))
}
draw_consensus <- function(tbl, title_text, filename) {
  sel <- select_heat_ids(tbl, 8)
  if (!nrow(sel)) stop("No protein-coding consensus genes for ", title_text)
  hm <- long |> dplyr::filter(ensembl_id %in% sel$ensembl_id) |>
    dplyr::group_by(symbol, contrast_id) |>
    dplyr::summarise(logFC = mean(logFC), .groups = "drop") |>
    tidyr::pivot_wider(names_from = contrast_id, values_from = logFC)
  col_order <- c(paste(CELL_ORDER, "TAM5", sep = "_"), paste(CELL_ORDER, "TAM12", sep = "_"))
  mat <- as.matrix(hm[, col_order, drop = FALSE]); rownames(mat) <- make.unique(hm$symbol)
  mean_effect <- rowMeans(mat, na.rm = TRUE)
  selected_response <- sel$response[match(rownames(mat), sel$symbol)]
  row_split <- factor(ifelse(selected_response == "Up", "Upregulated", "Downregulated"),
                      levels = c("Upregulated", "Downregulated"))
  ord <- order(row_split, -abs(mean_effect)); mat <- mat[ord, , drop = FALSE]; row_split <- row_split[ord]
  dose_ann <- factor(rep(c("5 μM", "12 μM"), each = 4), levels = c("5 μM", "12 μM"))
  cell_ann <- rep(CELL_ORDER, 2)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Dose = dose_ann, Cell = cell_ann,
    col = list(Dose = c("5 μM" = dose_cols[["TAM5"]], "12 μM" = dose_cols[["TAM12"]]), Cell = cell_cols),
    annotation_name_gp = grid::gpar(fontsize = 7),
    annotation_legend_param = list(labels_gp = grid::gpar(fontsize = 7),
                                   title_gp = grid::gpar(fontsize = 7, fontface = "bold")))
  lim <- max(1, min(4, max(abs(mat), na.rm = TRUE)))
  grDevices::cairo_pdf(file.path(V2_MAIN, filename), width = 92 / 25.4, height = 82 / 25.4)
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
    mat, name = "log2FC",
    col = circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, column_split = dose_ann,
    top_annotation = ha, column_names_rot = 45, show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 7), row_title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
    column_title = NULL, border = TRUE,
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
                                labels_gp = grid::gpar(fontsize = 7))),
    merge_legends = TRUE, align_heatmap_legend = "global_center",
    align_annotation_legend = "global_center", heatmap_legend_side = "right",
    annotation_legend_side = "right", padding = grid::unit(c(3, 3, 3, 3), "mm"))
  grDevices::dev.off()
  write.csv(sel, file.path(V2_TAB, sub("pdf$", "csv", filename)), row.names = FALSE)
}
draw_consensus(rule_tbl, "Rule-based consensus", "main_C1_consensus_rule_heatmap.pdf")
draw_consensus(meta_tbl, "Stouffer integrated response", "main_C2_consensus_meta_heatmap.pdf")

## ----------------------------- Fig6D Hallmark (正文15条) -----------------------------
hall <- read.csv(file.path(DIR$tables, "GSEA_Hallmark_all.csv"), check.names = FALSE)
# 固定15条通路（4功能类别），全部至少在4/8比较中nominal P<0.05；非反向挑选
hall_15 <- data.frame(
  pathway = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_COMPLEMENT",
              "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  category = c(rep("Cell cycle", 4), rep("Metabolism", 4), rep("Immunity", 4), rep("Stress", 3)),
  stringsAsFactors = FALSE)
hall_category_order <- c("Cell cycle", "Metabolism", "Immunity", "Stress")
# 验证：每条至少4/8比较 nominal P<0.05
hall_check <- hall |> dplyr::filter(pathway %in% hall_15$pathway) |>
  dplyr::group_by(pathway) |> dplyr::summarise(n_p05 = sum(pval < .05, na.rm = TRUE), .groups = "drop")
if (any(hall_check$n_p05 < 4)) warning("Hallmark 15条中存在n_p05<4的通路: ",
  paste(hall_check$pathway[hall_check$n_p05 < 4], collapse = ", "))
clean_pathway <- function(x) {
  x |> stringr::str_remove("^HALLMARK_") |> stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_replace_all(c("Myc" = "MYC", "E2f" = "E2F", "G2m" = "G2M",
                               "Tnfa" = "TNFα", "Nfkb" = "NF-κB", "P53" = "p53",
                               "Mtorc1" = "mTORC1", " Via " = " via ")) |>
    stringr::str_wrap(width = 30)
}
contrast_lookup <- data.frame(contrast_id = CONTRAST_ORDER,
                              cell_line = factor(rep(CELL_ORDER, each = 2), levels = CELL_ORDER),
                              dose = factor(rep(c("TAM5", "TAM12"), 4), levels = c("TAM5", "TAM12")))
hall_15$ord <- seq_len(nrow(hall_15))  # 固定通路顺序（按类别内列表顺序）
hall_main <- hall |> dplyr::filter(pathway %in% hall_15$pathway) |>
  dplyr::left_join(hall_15, by = "pathway") |>
  dplyr::mutate(pathway_label = clean_pathway(pathway),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                category = factor(category, levels = hall_category_order),
                pathway_label = factor(pathway_label, levels = rev(clean_pathway(hall_15$pathway[order(hall_15$ord)]))),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
hall_plot <- hall_main |> dplyr::mutate(contrast_id = as.character(contrast_id)) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pD <- ggplot2::ggplot(hall_plot, ggplot2::aes(dose, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .45) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.4,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(category ~ cell_line, scales = "free", space = "free", switch = "both") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_text(size = 7), strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.background.y = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 strip.text.y.left = ggplot2::element_text(angle = 90, size = 7),
                 panel.spacing = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .35),
                 plot.margin = ggplot2::margin(2, 3, 2, 2, "mm"))
save_v2(pD, "main_D_Hallmark_NES.pdf", 135, 86)
write.csv(hall_check, file.path(V2_TAB, "Hallmark_15pathways_n_p05_check.csv"), row.names = FALSE)
write.csv(hall_main, file.path(V2_TAB, "Fig6D_Hallmark_15pathways_data.csv"), row.names = FALSE)

## Full Hallmark heatmap (supp)
hall_rank <- hall |> dplyr::group_by(pathway) |>
  dplyr::summarise(mean_abs_NES = mean(abs(NES), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_abs_NES), pathway)
hall_full <- hall |>
  dplyr::mutate(pathway_label = clean_pathway(pathway),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                pathway_label = factor(pathway_label, levels = rev(clean_pathway(hall_rank$pathway))))
pDfull <- ggplot2::ggplot(hall_full, ggplot2::aes(contrast_id, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .2) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = rep(c("5 μM", "12 μM"), 4)) +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 5.2), legend.position = "right")
save_v2(pDfull, "FigS_full_Hallmark_NES_heatmap.pdf", 155, 190, where = V2_SUPP)

## ------------------- Fig6D-KEGG: KEGG 正文热图（同样式） -------------------
gk <- read.csv(file.path(DIR$tables, "GSEA_KEGG_LEGACY_all.csv"), check.names = FALSE)
kegg_main_list <- data.frame(
  pathway = c("KEGG_CELL_CYCLE", "KEGG_DNA_REPLICATION", "KEGG_OOCYTE_MEIOSIS", "KEGG_SPLICEOSOME",
              "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
              "KEGG_PENTOSE_PHOSPHATE_PATHWAY", "KEGG_PURINE_METABOLISM",
              "KEGG_RIBOSOME", "KEGG_PROTEASOME",
              "KEGG_LYSOSOME", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
              "KEGG_CHEMOKINE_SIGNALING_PATHWAY"),
  label = c("Cell cycle", "DNA replication", "Oocyte meiosis", "Spliceosome",
            "Oxidative phosphorylation", "Glycolysis / gluconeogenesis",
            "Pentose phosphate pathway", "Purine metabolism",
            "Ribosome", "Proteasome",
            "Lysosome", "Cytokine-receptor interaction",
            "Chemokine signaling pathway"),
  category = c(rep("Cell cycle", 4), rep("Metabolism", 4),
               rep("Ribosome", 2), rep("Stress", 3)),
  stringsAsFactors = FALSE)
kegg_category_order <- c("Cell cycle", "Metabolism", "Ribosome", "Stress")
kegg_check <- gk |> dplyr::filter(pathway %in% kegg_main_list$pathway) |>
  dplyr::group_by(pathway) |> dplyr::summarise(n_p05 = sum(pval < .05, na.rm = TRUE), .groups = "drop")
if (any(kegg_check$n_p05 < 4)) warning("KEGG正文中存在n_p05<4的通路: ",
  paste(kegg_check$pathway[kegg_check$n_p05 < 4], collapse = ", "))
kegg_main_list$ord <- seq_len(nrow(kegg_main_list))
kegg_main <- gk |> dplyr::filter(pathway %in% kegg_main_list$pathway) |>
  dplyr::left_join(kegg_main_list, by = "pathway") |>
  dplyr::mutate(pathway_label = factor(label, levels = rev(kegg_main_list$label[order(kegg_main_list$ord)])),
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                category = factor(category, levels = kegg_category_order),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
kegg_plot <- kegg_main |> dplyr::mutate(contrast_id = as.character(contrast_id)) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pDK <- ggplot2::ggplot(kegg_plot, ggplot2::aes(dose, pathway_label, fill = NES)) +
  ggplot2::geom_tile(colour = "white", linewidth = .45) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.4,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(category ~ cell_line, scales = "free", space = "free", switch = "both") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_text(size = 7), strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.background.y = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 strip.text.y.left = ggplot2::element_text(angle = 90, size = 6.5),
                 panel.spacing = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .35),
                 plot.margin = ggplot2::margin(2, 3, 2, 5, "mm"))
save_v2(pDK, "main_D_KEGG_NES.pdf", 140, 82)
write.csv(kegg_check, file.path(V2_TAB, "KEGG_main_pathways_n_p05_check.csv"), row.names = FALSE)
write.csv(kegg_main, file.path(V2_TAB, "Fig6D_KEGG_main_pathways_data.csv"), row.names = FALSE)

## ------------------- Fig6E: HSC-like signature 逆转（仅保留1个） -------------------
sig_res <- read.csv(file.path(DIR$tables, "Fig6F_signature_GSEA_all.csv"), check.names = FALSE)
sig_map <- data.frame(
  pathway = "HSC_LIKE_MARKERS",
  label = "HSC-like markers",
  stringsAsFactors = FALSE)
sig_main <- sig_res |> dplyr::inner_join(sig_map, by = "pathway") |>
  dplyr::mutate(reversal_score = NES * orientation,
                contrast_id = factor(contrast_id, levels = CONTRAST_ORDER),
                label = factor(label, levels = rev(sig_map$label)),
                star = dplyr::case_when(pval < 1e-4 ~ "****", pval < 1e-3 ~ "***",
                                        pval < 1e-2 ~ "**", pval < .05 ~ "*", TRUE ~ ""))
sig_plot <- sig_main |> dplyr::mutate(contrast_id = as.character(contrast_id),
                                      label = factor(label, levels = rev(sig_map$label))) |>
  dplyr::left_join(contrast_lookup, by = "contrast_id")
pE <- ggplot2::ggplot(sig_plot, ggplot2::aes(dose, label, fill = reversal_score)) +
  ggplot2::geom_tile(colour = "white", linewidth = .7) +
  ggplot2::geom_text(ggplot2::aes(label = star), colour = "black", size = 2.8,
                     family = "Arial", fontface = "bold", vjust = .7) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, name = "Reversal score") +
  ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
  ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
  ggplot2::labs(x = NULL, y = NULL) + theme_v2(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, margin = ggplot2::margin(t = 2)),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 strip.placement = "outside",
                 strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                 strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2, b = 2)),
                 panel.spacing.x = grid::unit(1.0, "mm"), legend.position = "right",
                 legend.title = ggplot2::element_text(size = 8), legend.box.just = "center",
                 panel.border = ggplot2::element_rect(fill = NA, linewidth = .4),
                 plot.margin = ggplot2::margin(2, 3, 2, 2, "mm"))
save_v2(pE, "main_E_HSC_reversal.pdf", 125, 28)
write.csv(sig_main, file.path(V2_TAB, "Fig6E_HSC_like_reversal_scores.csv"), row.names = FALSE)

## ----------------- Supplementary: KEGG_LEGACY GSEA heatmap -----------------
gk <- read.csv(file.path(DIR$tables, "GSEA_KEGG_LEGACY_all.csv"), check.names = FALSE)
gk_rank <- gk |> dplyr::group_by(pathway) |>
  dplyr::summarise(n_fdr = sum(padj < .05, na.rm = TRUE),
                   meanNES = mean(NES, na.rm = TRUE),
                   mean_abs = mean(abs(NES), na.rm = TRUE), .groups = "drop")
gk_sel <- gk_rank |> dplyr::filter(n_fdr >= 2) |>
  dplyr::arrange(dplyr::desc(n_fdr), dplyr::desc(mean_abs)) |>
  dplyr::slice_head(n = 60)
clean_kegg <- function(x) x |> stringr::str_remove("^KEGG_") |> stringr::str_replace_all("_", " ")
if (nrow(gk_sel)) {
  gk_df <- gk |> dplyr::filter(pathway %in% gk_sel$pathway) |>
    dplyr::mutate(plabel = clean_kegg(pathway),
                  plabel = factor(plabel, levels = rev(clean_kegg(
                    gk_sel |> dplyr::arrange(meanNES) |> dplyr::pull(pathway)))),
                  contrast_id = as.character(contrast_id)) |>
    dplyr::left_join(contrast_lookup, by = "contrast_id")
  pK <- ggplot2::ggplot(gk_df, ggplot2::aes(dose, plabel, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = .2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "NES") +
    ggplot2::scale_x_discrete(labels = c("TAM5" = "5 μM", "TAM12" = "12 μM")) +
    ggplot2::facet_grid(~ cell_line, scales = "free_x", space = "free_x", switch = "x") +
    ggplot2::labs(x = NULL, y = NULL) + theme_v2(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 7),
                   axis.text.y = ggplot2::element_text(size = 5.6),
                   strip.placement = "outside",
                   strip.background.x = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .3),
                   strip.text.x = ggplot2::element_text(size = 8),
                   panel.spacing.x = grid::unit(1, "mm"),
                   legend.position = "right")
  h_k <- max(60, nrow(gk_sel) * 4.0)
  save_v2(pK, "FigS_KEGG_LEGACY_GSEA_NES_heatmap.pdf", 170, h_k, where = V2_SUPP)
  write.csv(gk_sel, file.path(V2_TAB, "KEGG_LEGACY_GSEA_selected_pathways.csv"), row.names = FALSE)
}

## ----------------------- Supplementary: ORA dotplots -----------------------
ora_top <- read.csv(file.path(DIR$tables, "ORA_sig_p0.05_top10.csv"), check.names = FALSE)
ora_top$contrast_id <- factor(ora_top$contrast_id, levels = CONTRAST_ORDER)
ora_top$neglog10p <- -log10(pmax(ora_top$pvalue, 1e-300))
clean_ora_term <- function(x, coll) {
  x <- if (coll == "KEGG") sub("^KEGG_", "", x)
       else if (coll == "Hallmark") sub("^HALLMARK_", "", x) else x
  stringr::str_wrap(gsub("_", " ", x), width = 40)
}
ora_panel_levels <- paste(rep(CELL_ORDER, each = 2), rep(c("5 μM", "12 μM"), 4))
for (coll in c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Hallmark")) {
  df <- ora_top |> dplyr::filter(collection == coll)
  if (!nrow(df)) next
  # ORA 结果表本身已含 cell_line/dose；构造短面板标题（细胞系 + 剂量，固定顺序）
  df <- df |>
    dplyr::mutate(panel_lab = paste(cell_line, ifelse(dose == "TAM5", "5 μM", "12 μM")),
                  panel_lab = factor(panel_lab, levels = ora_panel_levels),
                  term_w = clean_ora_term(Description, coll))
  lv <- df |> dplyr::group_by(term_w) |>
    dplyr::summarise(mp = min(pvalue), .groups = "drop") |>
    dplyr::arrange(mp) |> dplyr::pull(term_w) |> unique()
  df$term_w <- factor(df$term_w, levels = rev(lv))
  p <- ggplot2::ggplot(df, ggplot2::aes(GeneRatioNum, term_w, colour = neglog10p, size = Count)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#F4A582", high = "#B2182B", name = expression(-log[10] * italic(P))) +
    ggplot2::scale_size_continuous(range = c(1.0, 3.2), name = "Count") +
    ggplot2::scale_x_continuous(n.breaks = 4, labels = scales::percent_format(accuracy = 1)) +
    ggplot2::facet_wrap(~ panel_lab, ncol = 2, scales = "free") +
    ggplot2::labs(x = "Gene ratio", y = NULL, title = NULL) + theme_v2(7) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.8),
                   axis.text.x = ggplot2::element_text(size = 6),
                   strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   panel.border = ggplot2::element_rect(fill = NA, linewidth = .3),
                   panel.spacing = grid::unit(3, "mm"),
                   legend.position = "right")
  hgt <- nrow(df) * 4.6 / 2 + 24
  save_v2(p, paste0("FigS_ORA_", coll, "_dotplot.pdf"), 180, max(120, hgt), where = V2_SUPP)
}

## --------------------------- Supplementary PCA ---------------------------
pca_list <- list()
for (id in CONTRAST_ORDER) {
  bits <- strsplit(id, "_", fixed = TRUE)[[1]]; dose_i <- tail(bits, 1)
  cell_i <- sub(paste0("_", dose_i, "$"), "", id)
  obj <- readRDS(file.path(DIR$deg, paste(cell_i, dose_i, sep = "_"), "limma_voom_object.rds"))
  e <- obj$voom$E
  vv <- matrixStats::rowVars(e)
  use <- order(vv, decreasing = TRUE)[seq_len(min(500, length(vv)))]
  pc <- stats::prcomp(t(e[use, , drop = FALSE]), center = TRUE, scale. = FALSE)
  d <- cbind(obj$metadata, as.data.frame(pc$x[, 1:2, drop = FALSE]))
  d$group <- factor(as.character(d$dose), levels = c("Control", dose_i))
  ve <- summary(pc)$importance[2, 1:2] * 100
  pca_list[[id]] <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = group)) +
    ggplot2::geom_point(size = 2.2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = replicate), size = 2.2, show.legend = FALSE, min.segment.length = 0) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#555555", dose_cols[[dose_i]]), c("Control", dose_i))) +
    ggplot2::labs(title = id, x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_v2(7) + ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"), legend.position = "top")
}
pca_panel <- patchwork::wrap_plots(pca_list, ncol = 4, guides = "collect") & ggplot2::theme(legend.position = "top")
save_v2(pca_panel, "FigS_eight_matched_PCA.pdf", 180, 105, where = V2_SUPP)

## ------------------------- Supplementary volcano -------------------------
vol_list <- lapply(CONTRAST_ORDER, function(id) {
  z <- deg[[id]]
  la <- z |> dplyr::filter(direction != "NS", !is.na(gene_symbol), gene_symbol != "", !grepl("^ENSG", gene_symbol)) |>
    dplyr::group_by(direction) |> dplyr::slice_min(P.Value, n = 2, with_ties = FALSE) |> dplyr::ungroup()
  lf <- z |> dplyr::filter(direction != "NS", gene_symbol %in% FOCUS_GENES)
  lb <- dplyr::bind_rows(la, lf) |> dplyr::distinct(gene_symbol, .keep_all = TRUE)
  ggplot2::ggplot(z, ggplot2::aes(logFC, -log10(pmax(P.Value, 1e-300)), colour = direction)) +
    ggplot2::geom_point(size = .32, alpha = .5) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = .3) +
    ggplot2::geom_hline(yintercept = -log10(.05), linetype = 2, linewidth = .3) +
    ggrepel::geom_text_repel(data = lb, ggplot2::aes(label = gene_symbol), size = 1.8,
                             min.segment.length = 0, max.overlaps = Inf) +
    ggplot2::scale_colour_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#C7C7C7")) +
    ggplot2::labs(title = id, x = expression(log[2] * " fold change"),
                  y = expression(-log[10] * " nominal " * italic(P))) + theme_v2(7) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"), legend.position = "none")
})
vol_panel <- patchwork::wrap_plots(vol_list, ncol = 4)
save_v2(vol_panel, "FigS_eight_volcano_plots.pdf", 180, 115, where = V2_SUPP)

## ------------------- Supplementary dose concordance -------------------
conc <- lapply(CELL_ORDER, function(cell_i) {
  a <- deg[[paste(cell_i, "TAM5", sep = "_")]] |>
    dplyr::select(ensembl_id, gene_symbol, lfc5 = logFC, p5 = P.Value)
  b <- deg[[paste(cell_i, "TAM12", sep = "_")]] |>
    dplyr::select(ensembl_id, lfc12 = logFC, p12 = P.Value)
  dplyr::inner_join(a, b, by = "ensembl_id") |>
    dplyr::mutate(cell_line = cell_i,
                  concordant_sig = abs(lfc5) > 1 & p5 < .05 & abs(lfc12) > 1 & p12 < .05 & sign(lfc5) == sign(lfc12))
}) |> dplyr::bind_rows()
rho <- conc |> dplyr::group_by(cell_line) |>
  dplyr::summarise(rho = cor(lfc5, lfc12, method = "spearman", use = "complete.obs"), .groups = "drop")
conc <- dplyr::left_join(conc, rho, by = "cell_line")
pconc <- ggplot2::ggplot(conc, ggplot2::aes(lfc5, lfc12)) +
  ggplot2::geom_bin_2d(bins = 60) +
  ggplot2::geom_point(data = dplyr::filter(conc, concordant_sig), shape = 1, size = .45, colour = "#B2182B", alpha = .65) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .35) +
  ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = .25) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = .25) +
  ggplot2::facet_wrap(~ cell_line, scales = "free") +
  ggplot2::geom_text(data = dplyr::distinct(conc, cell_line, rho),
                     ggplot2::aes(x = -Inf, y = Inf, label = sprintf("Spearman ρ = %.2f", rho)),
                     inherit.aes = FALSE, hjust = -.08, vjust = 1.3, size = 2.4) +
  ggplot2::scale_fill_gradient(low = "#F1F5F9", high = "#2166AC", name = "Gene density") +
  ggplot2::labs(x = expression("TAM 5 μM log"[2] * "FC vs matched control"),
                y = expression("TAM 12 μM log"[2] * "FC vs matched control")) +
  theme_v2(8) + ggplot2::theme(legend.position = "right")
save_v2(pconc, "FigS_TAM5_TAM12_effect_concordance.pdf", 150, 120, where = V2_SUPP)
write.csv(rho, file.path(V2_TAB, "TAM5_TAM12_effect_spearman.csv"), row.names = FALSE)

writeLines(c(
  "06_analysis_revised figure generation completed (v2: Hallmark15 + KEGG main + HSC-only).",
  "01/02 rebuilt from raw counts; DEG threshold |log2FC|>1 & nominal P<0.05 (same as 03).",
  "ORA: merged up+down DEGs; universe = all tested genes; significance raw pvalue<0.05; top10 shown; GO BP/CC/MF + KEGG(MSigDB KEGG_LEGACY) + Hallmark.",
  "GSEA: Hallmark (H, 50 sets) + KEGG (C2/CP:KEGG_LEGACY, 186 sets); fgseaMultilevel minSize=10 maxSize=500; ranking statistic = moderated t.",
  "Fig6D Hallmark: 15 fixed pathways in 4 categories (Cell cycle/Metabolism/Immunity/Stress); Fig6D-KEGG: 13 main pathways in 4 categories.",
  "Fig6E: HSC_LIKE_MARKERS only (orientation=-1, red=suppressed/reversed).",
  "Asterisks in Fig6D/Fig6E indicate nominal P (pval): * <0.05, ** <0.01, *** <0.001, **** <0.0001; full FDR in supplementary tables."
), file.path(V2, "RUN_SUMMARY.txt"))
message("06 figures done. Main: ", V2_MAIN, " | Supp: ", V2_SUPP)

# END ORIGINAL SCRIPT T11
