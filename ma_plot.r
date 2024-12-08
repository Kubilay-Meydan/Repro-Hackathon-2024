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
hist_plot <- ggplot(as.data.frame(res$padj), aes(x = res$padj)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Adjusted p-values",
       x = "Adjusted p-value",
       y = "Frequency")

ggsave("Hist_pvalAdj.pdf", plot = hist_plot, width = 7, height = 8)


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
  theme_gray(base_size = 14) +  # Couleur du fond
  labs(title = "MA-plot of complete RNA-seq dataset",
      x = "Mean of normalized counts",
      y = "Log2 Fold Change") +
  theme(panel.grid.major = element_line(color = "white"),  # Grille
      panel.grid.minor = element_line(color = "white"),  # Grille
      legend.position = "top")

ggsave("MA_plot_complete_dataset.png", plot = plot)

##### Volcano-plot complete dataset
res_df <- as.data.frame(res)

# Identification des gènes significatifs différentiellement exprimés
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

#Identification des gènes up-regulated et down-regulated
res_df$regulation <- ifelse(res_df$log2FoldChange > 1 & res$padj < 0.05, "Up-regulated", 
                            ifelse(res_df$log2FoldChange < -1 & res$padj < 0.05, "Down-regulated", "Not Significant"))

# Construction du graphique
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
                        geom_point(alpha = 0.7, size = 2) +  # Modifier la taille et la transparence des points
                        scale_color_manual(values = c("Not Significant" = "gray", "Up-regulated" = "blue", "Down-regulated" = "red"), na.translate = FALSE) +
                        theme_minimal() +
                        labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)") +
                        theme(legend.position = "right", plot.title = element_text(size = 16),
                              axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
                        xlim(c(-5, 5)) +  # Axe des abscisses
                        ylim(c(0, 50)) +  # Axe des ordonnées
                        # Délimitation verticale
                        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
                        # Délimitation horizontale
                        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

ggsave("Volcano_plot.pdf", plot = volcano_plot)


##### ACP 

# Transformation des données
rld <- rlog(dds)

# Extraction des données de plotPCA
pca_data <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)

# Calculer le pourcentage de variance expliquée
percent_variance <- round(100 * attr(pca_data, "percentVar"))

# Modifier les noms des échantillons (couper avant le ".")
pca_data$name <- sapply(pca_data$name, function(x) strsplit(x, "\\.")[[1]][1])

ACP <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = name)) +
                geom_point(size = 4, alpha = 0.8) +                     
                geom_text(vjust = -1, hjust = 0.5, size = 3) +         
                labs(title = "Principal Component Analysis of samples (PCA)", x = paste0("PC1: ", percent_variance[1], "% variance"),  
                      y = paste0("PC2: ", percent_variance[2], "% variance"), color = "Condition") +                            
                scale_x_continuous(limits = c(-25, 21)) +  
                scale_y_continuous(limits = c(-15, 10)) +  
                theme(legend.position = "bottom", plot.title = element_text(size = 20),    
                      axis.title = element_text(size = 16), axis.text = element_text(size = 14),                   
                      legend.text = element_text(size = 14), legend.title = element_text(size = 16)) 

ggsave("ACP_samples.png", plot = ACP)

