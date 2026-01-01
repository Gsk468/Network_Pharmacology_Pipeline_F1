
# --- R Cell ---
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(repr)

# --- R Cell ---
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.xlsx( "05_ppi/selected_proteins.xlsx", rowNames = F,colNames = T)
print(info)
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
print(gene)

# --- R Cell ---
#设置导出文件夹路径
output_dir <- "06_enrichment"
#如果文件夹不存在，则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
#进行GO富集分析
GO <- enrichGO(gene$ENTREZID,
               OrgDb = GO_database,
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
print(GO)
#将结果保存为CSV文件
output_file <- file.path(output_dir, "GO_enrichment_results.csv")
write.csv(GO, file = output_file)

# --- R Cell ---
# 设置导出文件夹路径
output_dir <- "06_enrichment"


KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
# Ensure readable gene symbols
KEGG <- setReadable(KEGG, OrgDb = GO_database, keyType = "ENTREZID")
print(KEGG)
# 将结果保存为CSV文件
output_file <- file.path(output_dir, "KEGG_enrichment_results.csv")
write.csv(KEGG, file = output_file)

# --- R Cell ---
#label_format=150表示标签长度为150，即每行显示150个字符，避免自动换行
p1 <- barplot(GO, split="ONTOLOGY", label_format=150)+facet_grid(ONTOLOGY~., scale="free")#柱状图
p2 <- barplot(KEGG,showCategory = 30,title = 'KEGG Pathway', label_format=150)
p3 <- dotplot(GO, split="ONTOLOGY", label_format=150)+facet_grid(ONTOLOGY~., scale="free")#点状图
p4 <- dotplot(KEGG, label_format=150)
p1
p2
p3
p4


# --- R Cell ---
p5 <- enrichplot::cnetplot(GO, circular = FALSE, colorEdge = FALSE,
         shadowtext = "all",
         color.params = list(foldChange = NULL, edge = FALSE, category = "#e13a38", gene = "#2fa7f3"))
p6 <- enrichplot::cnetplot(KEGG, circular = FALSE, colorEdge = FALSE,
         shadowtext = "all",
         color.params = list(foldChange = NULL, edge = FALSE, category = "#e13a38", gene = "#2fa7f3"))
p5
p6

# --- R Cell ---
p7 <- enrichplot::heatplot(GO,showCategory = 50, label_format=150)#基因-通路关联热图
p8 <- enrichplot::heatplot(KEGG,showCategory = 50, label_format=150)
p7
p8

# --- R Cell ---
GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
p9 <- enrichplot::emapplot(GO2,showCategory = 15, color = "p.adjust", layout = "kk", cex_label_category = 1, cex_line = 0.3)#通路间关联网络图
p10 <- enrichplot::emapplot(KEGG2,showCategory = 15, color = "p.adjust", layout = "kk", cex_label_category = 1, cex_line = 0.3)
p9
p10

# --- R Cell ---
GO_BP<-enrichGO( gene$ENTREZID,#GO富集分析BP模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
p11 <- plotGOgraph(GO_BP)#GO-BP功能网络图
GO_CC<-enrichGO( gene$ENTREZID,#GO富集分析CC模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
p12 <- plotGOgraph(GO_CC)#GO-CC功能网络图
GO_MF<-enrichGO( gene$ENTREZID,#GO富集分析MF模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
p13 <- plotGOgraph(GO_MF)#GO-MF功能网络图
p11
p12
p13

# --- R Cell ---
output_dir <- "07_figure"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 设置画布宽度和高度（以英寸为单位）
width <- 10
height <- 10

# 保存GO柱状图
ggsave(file.path(output_dir, "GO_barplot.svg"), p1, device = "svg", width = width, height = height)

# 保存KEGG柱状图
ggsave(file.path(output_dir, "KEGG_barplot.svg"), p2, device = "svg", width = width, height = height)

# 保存GO点状图
ggsave(file.path(output_dir, "GO_dotplot.svg"), p3, device = "svg", width = width, height = height)

# 保存KEGG点状图
ggsave(file.path(output_dir, "KEGG_dotplot.svg"), p4, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "GO_cnetplot.svg"), p5, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "KEGG_cnetplot.svg"), p6, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "GO_heatplot.svg"), p7, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "KEGG_heatplot.svg"), p8, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "GO_emapplot.svg"), p9, device = "svg", width = width, height = height)

ggsave(file.path(output_dir, "KEGG_emapplot.svg"), p10, device = "svg", width = width, height = height)

# --- Sankey Diagram (Drug-Target-Pathway) ---
library(ggalluvial)
library(dplyr)

# Load Drug-Target Interactions
dt_file <- "05_ppi/drug_target_interactions.csv"
if (file.exists(dt_file)) {
  dt_data <- read.csv(dt_file)

  # Load Pathway-Target (KEGG results)
  kegg_res_file <- file.path(output_dir, "KEGG_enrichment_results.csv")
  if (file.exists(kegg_res_file)) {
    kegg_data <- read.csv(kegg_res_file)

    # Process KEGG data to get Target-Pathway pairs
    tp_edges <- data.frame(Target = character(), Pathway = character(), stringsAsFactors = FALSE)

    # Only top 10 pathways to avoid clutter
    top_kegg <- head(kegg_data, 10)

    for (i in 1:nrow(top_kegg)) {
      pathway <- top_kegg$Description[i]
      genes <- strsplit(as.character(top_kegg$geneID[i]), "/")[[1]]
      for (gene in genes) {
        tp_edges <- rbind(tp_edges, data.frame(Target = gene, Pathway = pathway))
      }
    }

    # Merge Drug-Target and Target-Pathway
    # Join by Target
    sankey_data <- merge(dt_data, tp_edges, by = "Target")

    if (nrow(sankey_data) > 0) {
      # Plot Sankey
      p_sankey <- ggplot(sankey_data,
                         aes(axis1 = Drug, axis2 = Target, axis3 = Pathway)) +
        geom_alluvium(aes(fill = Drug), width = 1/12) +
        geom_stratum(width = 1/12, fill = "black", color = "grey") +
        geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("Drug", "Target", "Pathway"), expand = c(.05, .05)) +
        scale_fill_viridis_d() +
        theme_minimal() +
        ggtitle("Drug-Target-Pathway Sankey Diagram")

      ggsave(file.path(output_dir, "Sankey_Diagram.svg"), p_sankey, device = "svg", width = 12, height = 8)
    }
  }
}

# --- Pathway Chord Plot (pathway_ccplot) ---
library(GOplot)

# Prepare Gene Data (ID, logFC)
# We use centrality 'sum' as proxy for logFC for coloring
genedata <- data.frame(ID = info$SYMBOL, logFC = info$sum)

# Prepare Term Data (category, ID, term, adj_pval, genes)
# KEGG result has: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
if (exists("KEGG") && !is.null(KEGG)) {
  kegg_df <- as.data.frame(KEGG)
  if (nrow(kegg_df) > 0) {
    # GOplot needs 'genes' column to be comma separated? No, circle_dat merges.
    # circle_dat(terms, genes)
    # terms: category, ID, term, adj_pval
    # genes: ID, logFC
    # But how does it know which genes belong to which term?
    # Actually, circle_dat merges if terms has 'genes' column?
    # Standard usage: circle_dat(terms, genes)
    # terms should have 'genes' column?
    # Let's check documentation or common usage.
    # Usually circle_dat builds the object from enrichment result + gene data.
    # If we use clusterProfiler result, we might need conversion.

    # Simple approach: Construct the data frame expected by GOChord directly if possible, or use circle_dat.
    # circle_dat expects 'genes' column in 'terms' data frame to be a list of gene symbols?
    # Or 'geneID' column (slash separated).
    # GOplot::circle_dat tries to match.

    # Format for circle_dat:
    # terms: category, ID, term, adj_pval, genes (data frame)
    #   genes column: Uppercase gene symbols
    # Let's try to adapt kegg_df

    terms_df <- data.frame(
      category = "KEGG",
      ID = kegg_df$ID,
      term = kegg_df$Description,
      adj_pval = kegg_df$p.adjust,
      genes = gsub("/", ",", kegg_df$geneID), # GOplot expects comma separated
      stringsAsFactors = FALSE
    )

    # GOplot::circle_dat
    tryCatch({
      circ <- circle_dat(terms_df, genedata)

      # Plot GOChord
      # limit: number of terms
      # gene.order: 'logFC'
      p_chord <- GOChord(circ, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

      ggsave(file.path(output_dir, "pathway_ccplot.svg"), p_chord, device = "svg", width = 12, height = 12)
    }, error = function(e) {
      message("Error generating Chord plot: ", e$message)
    })
  }
}
