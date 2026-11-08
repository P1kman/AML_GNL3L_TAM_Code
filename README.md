# GNL3L–TAM analysis code in AML

Code accompanying *Integrative Bulk and Single-Cell Transcriptomics for Primitive Cell State–Informed Drug Repurposing in Acute Myeloid Leukemia*.

## Files

- `Public_Data_Analysis.R`: compilation of 22 original public-data analysis and plotting scripts.
- `TAM_RNAseq_Analysis.R`: compilation of 11 original TAM RNA-seq scripts, including focused reruns, the original driver, and a clearly identified historical figure-script backup.
- `README.md`: this guide and the source-section index.

## Preservation of the original scripts

Every original script is included once. Original script contents are preserved byte for byte, including comments, paths, thresholds, input readers, source() calls, and line endings. Only comment separators and a short compilation header were added outside those contents. Both files use UTF-8. Each section was checked against its original source file after writing.

This compilation uses the original local scripts, not the path-adjusted scripts in an earlier repository package. Previously added wrappers, recreated integration scripts, example inputs, figures, tables, and other supporting files are not included. No full analysis rerun or new result generation was performed for this compilation.

## Scope and inputs

The public-data scripts cover GSE30029 bulk expression analysis and GSE116256 single-cell analysis, including differential expression, enrichment, WGCNA/hdWGCNA, cell annotation, composition, CytoTRACE2, CellChat, gene-set intersection visualization, single-gene GSEA, and candidate-regulator analysis. The included scripts do not constitute an automated implementation of every external drug-resource query used in the study.

The TAM scripts begin with processed gene-count/annotation matrices and cover differential expression, consensus analysis, enrichment, signature reversal, and figure generation. FASTQ processing is outside this collection. Required matrices, curated objects, annotation files, database exports, and third-party helper scripts must be obtained separately according to the original read/load/source statements. They are not bundled here.

## How to use

**These are source-code compilations, not standalone one-command pipelines. Do not source either complete file as a single workflow.** Numbered sections follow analysis stages for navigation; this does not replace prerequisites, interactive annotation, or the original driver logic.

1. Open the relevant compilation in RStudio with UTF-8 encoding. Use the section IDs and line ranges below to find the original script.
2. Read its input paths, package imports, and prerequisites before running the relevant stage. Original machine-specific paths have deliberately not been changed. Adapt paths in your own working copy or reproduce the expected directory structure.
3. Original `source()` calls still refer to separate files. When needed, copy the referenced original-script section into a local file with its original name and expected relative path. External scripts not present in the index must be obtained separately. Merely running a helper section earlier does not make a later file-based `source()` call unnecessary.
4. TAM sections repeatedly load `R/00_config.R`. To use the original driver, restore the original TAM sections under the expected `R/` directory, supply the input data, and follow `run_all_revised.R`. This is a local execution option; the distributed repository itself contains only the three files listed above.
5. Focused reruns are alternatives for selected figures. The historical figure backup is retained for completeness, not as an additional stage of the current analysis. Public scripts may also contain repeated or historical interactive blocks; review their comments rather than executing every block indiscriminately.

Install packages used by the selected original scripts. Package imports and installation instructions, where present, are preserved in those sections; this package does not introduce a new dependency lockfile or certify an execution environment.

## Source-section index

Paths below are relative to the original public-analysis directory or TAM `R` directory. Line ranges include only original script contents, excluding added separator comments.

### Public_Data_Analysis.R

| Section | Original source | Lines in compilation |
|---|---|---|
| P01 | `00_pub_theme_colors.R` | 10–409 |
| P02 | `02_BULK_DIFF.R` | 416–934 |
| P03 | `03_BULK_GO_KEGG.R` | 941–1117 |
| P04 | `04_BULK_GSEA.R` | 1124–2429 |
| P05 | `05_BULK_WGCNA.R` | 2436–3443 |
| P06 | `07_01_SC_annotation_数据清洗.R` | 3450–4097 |
| P07 | `07_02_SC_annotation_聚类与注释.R` | 4104–4951 |
| P08 | `07_03_SC_annotation_名称修改与marker气泡图绘制.R` | 4958–5233 |
| P09 | `07_SC_annotation/1_marker_cluster_output.R` | 5240–5263 |
| P10 | `07_SC_annotation/1_new_cluster_ids_output.R` | 5270–5293 |
| P11 | `07_SC_annotation/plot_annotation_pub.R` | 5300–5651 |
| P12 | `08_SC_cell_ratio.R` | 5658–6901 |
| P13 | `08_SC_cell_ratio/plot_cell_ratio_pub.R` | 6908–7379 |
| P14 | `09_SC_CYTOTRACE2.R` | 7386–8752 |
| P15 | `09_SC_CYTOTRACE2/plot_cytotrace2_pub.R` | 8759–9518 |
| P16 | `10_SC_enrich.R` | 9525–9756 |
| P17 | `11_SC_cellchat.R` | 9763–10040 |
| P18 | `11_SC_cellchat/plot_cellchat_pub.R` | 10047–11060 |
| P19 | `12_SC_hdWGCNA.R` | 11067–11479 |
| P20 | `13_single_gene_GSEA.R` | 11486–12734 |
| P21 | `13_BULK_SC_gene_intersection_AND_single_gene_GSEA/BULK_SC_gene_upset_inputs/upset_evidence_plots.R` | 12741–13716 |
| P22 | `14_TF_CistromeDB_decoupleR.R` | 13722–14447 |

### TAM_RNAseq_Analysis.R

| Section | Original source | Lines in compilation |
|---|---|---|
| T01 | `00_config.R` | 10–139 |
| T02 | `01_prepare_and_DEG.R` | 146–241 |
| T03 | `02_consensus_analysis.R` | 248–301 |
| T04 | `03_ORA_GO_KEGG_Hallmark.R` | 308–440 |
| T05 | `04_GSEA_Hallmark_KEGG.R` | 447–498 |
| T06 | `05_signature_GSEA_revised.R` | 505–576 |
| T07 | `06_generate_figures.R` | 583–1138 |
| T08 | `rerun_main_B_only.R` | 1145–1243 |
| T09 | `rerun_main_E_only.R` | 1250–1323 |
| T10 | `run_all_revised.R` | 1330–1347 |
| T11 | `06_generate_figures_backup_before_3formats.R` | 1354–1880 |

## Upload to GitHub / 上传方法

Unzip the archive and upload the three files directly to the repository root using **Add file → Upload files**, then commit the changes. Upload the extracted files rather than only the ZIP, so reviewers can browse the code. If replacing an earlier package, remove obsolete repository files that are no longer intended for distribution.

解压后，将两个 R 文件和 README.md 直接上传到 GitHub 仓库根目录；不要只把 ZIP 上传成一个附件。原代码未改写，合并文件并不等于可从头一键运行的脚本。原有本地路径、source() 调用、环境清理和历史分支均保留，使用时应按对应原始分析步骤准备数据并运行。TAM 最后一节是历史绘图备份，不应与当前绘图脚本连续重复执行。

Public dataset identifiers are GSE30029 and GSE116256. No unassigned GEO accession is invented here for newly generated sequencing data; cite the final data availability statement when available. No software license has been newly assigned by this packaging step; third-party code and packages retain their applicable terms.
