### top 8 and NES and padj on the Vis plot 
# ---------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

has_gseavis <- requireNamespace("GseaVis", quietly = TRUE)

# Enforce R to register extra system font families if available
# This keeps system text mapping stable during device rendering
options(bitmapType = 'cairo')  

# ---------------------------------------------------------
# 2. Read input files
# ---------------------------------------------------------
file1_path <- "/mnt/localstorage/ramm/data/Hari_volcano/Jeehye/New_for_paper/ATF4_targets_CBlist/ATF4Targets.txt"
file2_path <- "/mnt/localstorage/ramm/data/Hari_volcano/Jeehye/Motorneurons/ATF4_targets_MNlist/DEG_TR_WT_Week4_6_100_Batch1_exclude_PCG.DEG.txt"

df1 <- read.delim(file1_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
df2 <- read.delim(file2_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

week_label <- regmatches(basename(file2_path), regexpr("Week\\d+", basename(file2_path)))
if (length(week_label) == 0 || is.na(week_label) || week_label == "") {
  week_label <- "Week4-6"
}

# ---------------------------------------------------------
# 3. Identify ATF4 target column
# ---------------------------------------------------------
candidate_cols <- c("GeneSymbol", "gene", "Gene", "SYMBOL", colnames(df1)[1])
target_col <- candidate_cols[candidate_cols %in% colnames(df1)][1]

if (is.na(target_col) || target_col == "") {
  stop("Could not identify gene symbol column in ATF4 target file.")
}

atf4_targets <- unique(trimws(as.character(df1[[target_col]])))
atf4_targets <- atf4_targets[!is.na(atf4_targets) & atf4_targets != ""]

cat("ATF4 target column used:", target_col, "\n")
cat("Number of ATF4 target genes:", length(atf4_targets), "\n")

# ---------------------------------------------------------
# 4. Prepare ranked gene list
# ---------------------------------------------------------
required_cols <- c("GeneSymbol", "log2FoldChange")
missing_cols <- setdiff(required_cols, colnames(df2))
if (length(missing_cols) > 0) {
  stop("Missing required columns in DEG file: ", paste(missing_cols, collapse = ", "))
}

gsea_df <- df2[!is.na(df2$GeneSymbol) & df2$GeneSymbol != "" & !is.na(df2$log2FoldChange), ]
gsea_df$absFC <- abs(gsea_df$log2FoldChange)
gsea_df <- gsea_df[order(gsea_df$absFC, decreasing = TRUE), ]
gsea_df <- gsea_df[!duplicated(gsea_df$GeneSymbol), ]

gene_list <- gsea_df$log2FoldChange
names(gene_list) <- gsea_df$GeneSymbol
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Ranked gene list size:", length(gene_list), "\n")

# ---------------------------------------------------------
# 5. Run GO BP GSEA
# ---------------------------------------------------------
cat("Running GO BP GSEA...\n")

gsea_results <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

gsea_res_df <- as.data.frame(gsea_results)

write.csv(
  gsea_res_df,
  paste0("GSEA_Results_Table_", week_label, ".csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 6. Custom ATF4 target GSEA
# ---------------------------------------------------------
term2gene <- data.frame(
  term = "ATF4_Targets",
  gene = atf4_targets,
  stringsAsFactors = FALSE
)
set.seed(42)
custom_gsea <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = term2gene,
  minGSSize     = 5,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

custom_res_df <- as.data.frame(custom_gsea)

write.csv(
  custom_res_df,
  paste0("GSEA_ATF4_Targets_", week_label, ".csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 7. Helper function
# ---------------------------------------------------------
save_gsea_plot <- function(gsea_obj, top_id, title_text, color, out_file,
                           genes_to_mark = NULL, width = 7.5, height = 6.2) {

  p <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID = top_id,
    title = title_text,
    color = color,
    subplots = 1:2,
    pvalue_table = FALSE,
    ES_geom = "line"
  )

  # Explicitly intercept and trace text elements to override with Arial family
  p <- p + theme(text = element_text(family = "Arial"), 
  axis.text.x  = element_text(size = 12),  # <--- Size of X-axis tick numbers
    axis.text.y  = element_text(size = 12),  # <--- Size of Y-axis tick numbers
    axis.title.x = element_text(size = 14),  # <--- Size of X-axis title/label
    axis.title.y = element_text(size = 14)   # <--- Size of Y-axis title/label
  )

  if (!is.null(genes_to_mark) && length(genes_to_mark) > 0) {
    genes_to_mark <- unique(genes_to_mark[!is.na(genes_to_mark) & genes_to_mark != ""])
    if (length(genes_to_mark) > 0) {
      p <- p + enrichplot::geom_gsea_gene(
        genes = genes_to_mark,
        geneSet = top_id,
        size = 8.5              
      )
    }
  }

  ggsave(
    filename = out_file,
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  invisible(p)
}

# ---------------------------------------------------------
# 8. Save GO BP plot with highlighted genes
# ---------------------------------------------------------
if (nrow(gsea_res_df) > 0) {
  top_id   <- gsea_res_df$ID[1]
  top_desc <- gsea_res_df$Description[1]
  top_nes  <- round(gsea_res_df$NES[1], 2)
  top_padj <- format(gsea_res_df$p.adjust[1], scientific = TRUE, digits = 3)

  # Dynamic core enrichment extraction for standard GO pathways
  go_core_str <- gsea_res_df$core_enrichment[1]
  if (!is.na(go_core_str) && go_core_str != "") {
    genes_to_mark <- head(unlist(strsplit(go_core_str, "/")), 8)
  } else {
    genes_to_mark <- NULL
  }

  cat("Generating GO GSEA plot for:", top_desc, "\n")

  go_title <- paste0(
    "ATF4 Targets\n",
    top_desc,
    "\nNES = ", top_nes,
    " | padj = ", top_padj
  )

  save_gsea_plot(
    gsea_obj      = gsea_results,
    top_id        = top_id,
    title_text    = go_title,
    color         = "#E6B800",
    out_file      = paste0("GSEA_PublicationPlot_GO_", week_label, ".png"),
    genes_to_mark = genes_to_mark,
    width         = 7.5,
    height        = 6.2
  )

} else {
  cat("No significant GO BP pathways found.\n")
}

# ---------------------------------------------------------
# 9. Save ATF4 target custom GSEA plot with highlighted genes
# ---------------------------------------------------------
if (nrow(custom_res_df) > 0) {
  top_custom_id   <- custom_res_df$ID[1]
  top_custom_desc <- if ("Description" %in% colnames(custom_res_df)) custom_res_df$Description[1] else top_custom_id
  top_custom_nes  <- round(custom_res_df$NES[1], 2)
  top_custom_padj <- format(custom_res_df$p.adjust[1], scientific = TRUE, digits = 3)

  cat("Generating custom ATF4 target GSEA plot...\n")

  custom_title <- paste0(
    "ATF4 targets\n",
    top_custom_desc,
    "\nNES = ", top_custom_nes,
    " | padj = ", top_custom_padj
  )

  # DYNAMIC CORE ENRICHMENT EXTRACTION FROM YOUR REAL RESULTS
  custom_core_str <- custom_res_df$core_enrichment[1]
  
  if (!is.na(custom_core_str) && custom_core_str != "") {
    all_core_genes <- unlist(strsplit(custom_core_str, "/"))
    genes_to_mark_custom <- head(all_core_genes, 8) 
  } else {
    genes_to_mark_custom <- NULL
  }

  save_gsea_plot(
    gsea_obj      = custom_gsea,
    top_id        = top_custom_id,
    title_text    = custom_title,
    color         = "purple4",
    out_file      = paste0("GSEA_PublicationPlot_ATF4Targets_", week_label, ".png"),
    genes_to_mark = genes_to_mark_custom,
    width         = 7.5,
    height        = 6.2
  )

} else {
  cat("No significant enrichment detected for custom ATF4 target set.\n")
}


# ---------------------------------------------------------
# 10. Optional GseaVis plot if gseaNb exists
# ---------------------------------------------------------
if (has_gseavis) {
  gseavis_exports = getNamespaceExports("GseaVis")

  if ("gseaNb" %in% gseavis_exports && nrow(custom_res_df) > 0) {
    cat("Saving optional GseaVis::gseaNb ATF4 plot...\n")

    png(
      filename = paste0("GSEA_PublicationPlot_ATF4Targets_GseaVis_NES", week_label, ".png"),
      width = 2200,
      height = 1800,
      res = 300
    )

    # 1. Format metrics for clean display text
    stat_label <- paste0("NES = ", top_custom_nes, "\np.adj = ", top_custom_padj)

    # 2. Build the base plot
    annotated_vis_plot <- GseaVis::gseaNb(
      object      = custom_gsea,
      geneSetID   = top_custom_id,
      subPlot     = 2,
      addGene     = TRUE,
      markTopgene = TRUE,
      topGeneN    = min(length(genes_to_mark_custom), 10),
      geneSize    = 8
    ) 
    
    # Inject the annotation text into the top curve layer (subplot 1)
    annotated_vis_plot[[1]] <- annotated_vis_plot[[1]] + 
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = stat_label, 
        hjust = 1.1,    
        vjust = 1.5,    
        family = "Arial", 
        fontface = "plain", 
        size = 6,        
        color = "black"
      )

    # 3. SET AXIS LABELS DIRECTLY ON SUBPLOTS AND STRIP MULTI-PANEL INHERITED TITLES
    # Wiping out title/subtitle layouts here clears out the nested subplot layers completely
    annotated_vis_plot[[1]] <- annotated_vis_plot[[1]] + 
      labs(title = NULL, subtitle = NULL, x = NULL, y = "Enrichment Score") +
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
      
    annotated_vis_plot[[2]] <- annotated_vis_plot[[2]] + 
      labs(title = NULL, subtitle = NULL, x = "Rank in Ordered Dataset", y = NULL) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank())

    # 4. REMOVE THE CACHED MAIN CONTAINER TITLE VIA PLOT_LAYOUT()
    # This prevents the pre-packaged "ATF4 Targets" layout title from sticking to the top
    final_output <- annotated_vis_plot + patchwork::plot_layout(guides = "collect") & 
      theme(
        text          = element_text(family = "Arial"),
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),
        axis.title.x  = element_text(size = 20), 
        axis.title.y  = element_text(size = 20)
      )

    # 5. GENERATE THE SINGLE Master TITLE FOR THE PUBLICATION
    custom_plot_title <- paste0("ATF4 Targets\n", top_custom_desc)
    
    final_output <- final_output + patchwork::plot_annotation(
      title = custom_plot_title,
      theme = theme(
        plot.title = element_text(family = "Arial", size = 22, face = "plain", hjust = 0.5)
      )
    )

    print(final_output)
    dev.off()
    
  } else {
    cat("GseaVis installed, but gseaNb not exported or no custom GSEA results.\n")
  }
} else {
  cat("GseaVis not installed; skipped optional GseaVis plot.\n")
}