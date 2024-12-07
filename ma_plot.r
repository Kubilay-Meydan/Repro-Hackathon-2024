library(DESeq2)
library(ggplot2)
packageVersion("DESeq2")


# On récupère la table et on nettoie pour garder juste les colonnes intéressantes
count_data <- read.table("merged_feature_counts.txt", header = TRUE, skip = 1, row.names = 1)

clean_count_data <- count_data[,c(6:11)]
clean_count_data <- as.matrix(clean_count_data)

#Partie pour créer le metadata
metadata <- data.frame(
  Sample = c("SRR10379726.fastq_trimmed.fastq.bam", "SRR10379723.fastq_trimmed.fastq.bam", "SRR10379722.fastq_trimmed.fastq.bam",
             "SRR10379725.fastq_trimmed.fastq.bam", "SRR10379724.fastq_trimmed.fastq.bam", "SRR10379721.fastq_trimmed.fastq.bam"),
  Condition = c("control", "persister", "persister", "control", "control", "persister")
)
#write.csv(metadata, "metadata.csv", row.names = FALSE)

#metadata <- read.csv("metadata.csv", row.names = 1)
rownames(metadata) <- colnames(clean_count_data)
metadata$Condition <- factor(metadata$Condition, levels = c("control", "persister"))

#class(clean_count_data)
#str(clean_count_data)
#class(metadata)
#str(metadata)

# Créer l'objet DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = clean_count_data,
                              colData = metadata,
                              design = ~ Condition)

# Exécuter DESeq
dds <- DESeq(dds)

# Obtenir les résultats
res <- results(dds)
summary(res)

# Vérification de la distribution des pvalues ajustées
png("Hist_pvalAdj.png", width = 800, height = 600, res = 300)
hist(res$padj)
dev.off()


##### MA-PLOT complete dataset

res_df <- as.data.frame(res)
# Identification des gènes différentiellements exprimés
res_df$Significance <- ifelse(res_df$padj < 0.05 & !is.na(res_df$padj), "Significant", "Not Significant")

# Représentation sous forme de MA-plot
plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
  geom_point(size = 1.2, alpha = 0.5) +  # Points (gènes)
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black"),
                     guide = "none") +  # Couleurs
  scale_x_log10(breaks = c(10^0, 10^2, 10^4),  # Axe des abscisses
                labels = c(expression(10^0), expression(10^2), expression(10^4))) +  # Labels
  coord_cartesian(ylim = c(-4, 4)) +  # Axe des ordonnées
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Ligne horizontale
  theme_gray(base_size = 14) +  # Grille de fond et taille du texte
  labs(title = "MA-plot of complete RNA-seq dataset",
      x = "Mean of normalized counts",
      y = "Log2 Fold Change") +
  theme(panel.grid.major = element_line(color = "white"),  # Grille majeure
      panel.grid.minor = element_line(color = "white"),  # Grille mineure
      legend.position = "top")

ggsave("plot_complete_dataset.png", plot = plot)

##### Volcano-plot complete dataset

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Catégoriser les gènes en fonction de leur régulation
res_df$Regulation <- "No Change"
res_df$Regulation[res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "Up-regulated"
res_df$Regulation[res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "Down-regulated"

# Significativité des pvalues ajustées
res_df$Significance <- "Not significant"
res_df$Significance[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant"

# Volcanon plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("red","gray","blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot of complete dataset", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("volcano_plot.png", plot = volcano_plot)