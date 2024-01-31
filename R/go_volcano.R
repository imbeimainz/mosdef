# library(AnnotationDbi)
# library(ggplot2)
# library("org.Hs.eg.db")
# library(rlang)
# 
# # res_de$symbol <- mapIds(org.Hs.eg.db,
# #                         keys = row.names(res_de),
# #                         column = "SYMBOL",
# #                         keytype = "ENSEMBL",
# #                         multiVals = "first"
# # )
# 
# df <- deseqresult2df(res_airway)
# # has to be removed when turning this into a functioj
# L2FC_cutoff = 1
# res_enrich <- topgoDE_airway
# 
# # finding the highest value in the log2FoldChange column and rounding it up to get a nice symetric plot
# x_limit <- ceiling(max(abs(range(df$log2FoldChange, na.rm = TRUE))))
# 
# 
# 
# df$diffexpressed <- "NO"
# # if log2Foldchange > Cutoff and pvalue < 0.05, set as "UP"
# df$diffexpressed[df$log2FoldChange > L2FC_cutoff & df$pvalue < 0.05] <- "UP"
# # if log2Foldchange < -Cutoff and pvalue < 0.05, set as "DOWN"
# df$diffexpressed[df$log2FoldChange < -L2FC_cutoff & df$pvalue < 0.05] <- "DOWN"
# 
# # calculate top 30 degenes based on pvalue
# #df$delabel <- ifelse(df$symbol %in% head(df[order(df$pvalue), "symbol"], labeled_genes), df$symbol, NA)
# 
# 
# # Either get indexes of the term name
# grep("cell adhesion", res_enrich[["Term"]], fixed =TRUE)
# 
# 
# #or have the user tell us at which index the term of interest is
# index_test <- 201
# test_vec <- res_enrich$genes[index_test]
# test_vec <- strsplit(test_vec, ",")
# test_vec <- as.vector(test_vec)
# test_for_plot <- test_vec[[1]][1:30]
# 
# for (i in 1:length(df$id)) {
# 
# 
#   if (df$symbol[i] %in% test_vec[[1]]){
#     df$de_label[i] <- df$symbol[i]
# 
# 
#   } else{
#     df$de_label[i] <- NA
#   }
# 
# 
# }
# #df$delabel <- ifelse(df$symbol %in% head(df[order(df$de_label), "symbol"], labeled_genes), df$symbol, NA)
# 
# p <-ggplot(data = df, aes(
#   x = .data$log2FoldChange, y = -log10(.data$pvalue),
#   colour = .data$diffexpressed, label = .data$de_label
# )) +
#   geom_vline(xintercept = c(-L2FC_cutoff, L2FC_cutoff), col = "gray", linetype = "dashed") +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(
#     values = c("skyblue", "grey", "tomato"), # to set the colours of our variable
#     labels = c("Downregulated", "Not significant", "Upregulated")
#   ) +
#   coord_cartesian(ylim = c(0, max(-log10(df$pvalue))), xlim = c(-x_limit, x_limit)) +
#   scale_x_continuous(breaks = seq(-x_limit, x_limit, 2)) +
#   geom_label_repel(max.overlaps = Inf) # To show all labels
# # or   geom_text_repel(max.overlaps = Inf) not quite sure which is better
# p
# 
# # Define the genes to highlight
# genes_to_highlight <- c(test_vec[[1]])
# # Define a custom color scale
# custom_color_scale <- c("DOWN" = "skyblue", "NO" = "grey", "UP" = "tomato", "Highlighted" = "green")
# 
# # Add another geom_point layer to highlight specific genes
# q <- p +
#   geom_point(data = df[df$symbol %in% genes_to_highlight, ],
#              aes(color = "Highlighted"), size = 3, shape = 16) +
#   scale_color_manual(values = custom_color_scale,
#                      labels = c("Downregulated", "GO term", "Not significant", "Upregulated")
#   )
# 
# # Show the plot
# print(q)
# 
# 
# #Ideas:
# # dont split into up and down (by color) but by goterm and not goterm in maybe red and black
# # add an alpha = 0.5 to the original alyer and alpha 1 to the second layer
# # somehow filter terms with more gnees for the most significant by padj
#   # maybe a loop that extraxts all genes from test vec out of the res_de together with their padj
#   # then filter by that padj and take top 30 + add a message maybe even a parameter
