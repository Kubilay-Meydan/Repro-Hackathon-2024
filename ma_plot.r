library(DESeq2)
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

pdf("MA_plot.pdf")
plotMA(res, main = "MA plot",ylab = "Log2 Fold Change",xlab ="Mean of normalized counts", ylim = c(-4.5, 4.5), colSig = "red", colNonSig = "black", colLine = NA, xlim = c(1e-01, 1e+06), log = "x")
# colLine en NA pour l'enlever et mettre celle en pointillée avec abline

abline(h = 0, col = "black", lty = 2)
dev.off()
