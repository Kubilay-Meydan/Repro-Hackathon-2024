############ GENE NAMES EXTRACTION #######################################################################
# Définition des chemins des fichiers
merged_feature_counts_path <- "merged_feature_counts.txt"
gene_specific_info_path <- "GeneSpecificInformation_NCTC8325.tsv"
output_path <- "Merged_counts_gene_names.txt"

# Lecture des fichiers
merged_feature_counts <- read.delim(merged_feature_counts_path, header = TRUE, sep = "\t", comment.char = "#")
gene_specific_info <- read.delim(gene_specific_info_path, header = TRUE, sep = "\t")

# Sélection des colonnes nécessaires dans gene_specific_info
gene_specific_info_reduced <- gene_specific_info[, c("locus.tag", "pan.gene.symbol")]

# Renommer les colonnes pour éviter des confusions
colnames(gene_specific_info_reduced) <- c("Geneid", "Gene_name")

# Initialiser la nouvelle colonne Gene_name avec NA
merged_feature_counts$Gene_name <- NA

# Remplir la colonne Gene_name par une boucle (pour compatibilité avec les anciennes versions de R)
for (i in 1:nrow(merged_feature_counts)) {
  match_index <- which(gene_specific_info_reduced$Geneid == merged_feature_counts$Geneid[i])
  if (length(match_index) > 0) {
    merged_feature_counts$Gene_name[i] <- gene_specific_info_reduced$Gene_name[match_index]
  }
}

# Sauvegarde du fichier fusionné
write.table(
  merged_feature_counts,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(paste("Fichier fusionné enregistré dans :", output_path, "\n"))


###########. TRANSLATION GENE KEGG API EXTRACTION  #########################################################################################
get_translation_genes <- function() {
  # URL pour accéder à la hiérarchie BRITE de KEGG pour "sao"
  url_brite <- "https://rest.kegg.jp/get/br:sao00001"
  
  # Télécharger les données BRITE
  raw_data <- readLines(url_brite, warn = FALSE)
  
  # Initialiser la liste des gènes
  translation_genes <- character()
  
  # Identifier les index de début et de fin pour la section "09122 Translation"
  translation_start <- grep("^B\\s+09122 Translation", raw_data)
  next_b <- grep("^B\\s+", raw_data)
  
  # Déterminer la fin de la section "09122 Translation"
  translation_end <- ifelse(length(next_b[next_b > translation_start[1]]) > 0,
                            min(next_b[next_b > translation_start[1]]) - 1,
                            length(raw_data))
  
  # Extraire les gènes pour "09122 Translation"
  if (length(translation_start) > 0) {
    for (i in seq(translation_start + 1, translation_end)) {
      line <- raw_data[i]
      if (grepl("^D\\s+SAOUHSC_", line)) {
        gene_id <- sub("^D\\s+(SAOUHSC_\\S+).*", "\\1", line)
        translation_genes <- c(translation_genes, gene_id)
      }
    }
  }
  
  # Vérifications intermédiaires
  cat("Nombre de gènes extraits de '09122 Translation':", length(translation_genes), "\n")
  
  return(translation_genes)
}

# Charger le fichier merged_feature_counts
merged_feature_counts_path <- "Merged_counts_gene_names.txt"
merged_feature_counts <- read.delim(merged_feature_counts_path, header = TRUE, sep = "\t")

# Vérifier les colonnes du fichier
cat("Colonnes dans Merged_counts_gene_names.txt:\n")
print(colnames(merged_feature_counts))

# Récupérer la liste des gènes associés à "09122 Translation"
translation_genes <- get_translation_genes()

# Vérifier combien de gènes extraits correspondent dans merged_feature_counts
matching_genes <- translation_genes[translation_genes %in% merged_feature_counts$Geneid]
cat("Nombre de gènes correspondant dans le fichier source:", length(matching_genes), "\n")

# Filtrer les lignes de merged_feature_counts pour ces gènes
filtered_counts <- merged_feature_counts[
  merged_feature_counts$Geneid %in% matching_genes, 
]

# Afficher les correspondances trouvées
cat("Gènes correspondants trouvés dans le fichier (premiers 5):\n")
print(head(filtered_counts))

# Sauvegarder la matrice filtrée dans un fichier
output_path <- "Translation_feature_counts.tsv"
write.table(
  filtered_counts,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(paste("Matrice filtrée enregistrée dans :", output_path, "\n"))





# Fichier d'entrée et de sortie
input_file <- "Translation_feature_counts.tsv"  # Nom du fichier source
output_file <- "tRNA_lines.tsv"                # Nom du fichier de sortie

# Lire le fichier source
data <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Vérifier les colonnes disponibles
cat("Colonnes disponibles dans le fichier:\n")
print(colnames(data))

# Filtrer les lignes où Gene_name se termine par "S" majuscule
filtered_data <- data[grepl("S$", data$Gene_name), ]

# Vérifier les lignes extraites
cat("Nombre de lignes extraites:", nrow(filtered_data), "\n")
cat("Aperçu des lignes extraites (premiers 5):\n")
print(head(filtered_data))

# Sauvegarder la matrice filtrée dans un fichier
write.table(
  filtered_data,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(paste("Matrice filtrée enregistrée dans :", output_file, "\n"))






################# DIFFERENTIAL EXPRESSION PLOT     ###################################################################################################



# Charger les bibliothèques nécessaires
library(DESeq2)
library(ggplot2)

# Chargement des métadonnées
metadata <- data.frame(
  Sample = c("SRR10379726.fastq_trimmed.fastq.bam", "SRR10379723.fastq_trimmed.fastq.bam", 
             "SRR10379722.fastq_trimmed.fastq.bam", "SRR10379725.fastq_trimmed.fastq.bam", 
             "SRR10379724.fastq_trimmed.fastq.bam", "SRR10379721.fastq_trimmed.fastq.bam"),
  Condition = c("control", "persister", "persister", "control", "control", "persister")
)
#write.csv(metadata, "metadata.csv", row.names = FALSE)

# Chargement des données de comptage
counts_file <- "Translation_feature_counts.tsv"
counts_data <- read.delim(counts_file, header = TRUE, row.names = 1)

# Lecture des métadonnées
#metadata <- read.csv("metadata.csv")

# Nettoyage des données de comptage pour ne garder que les colonnes pertinentes
clean_count_data <- counts_data[, c(6:11)]

# Création de l'objet DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = clean_count_data,
  colData = metadata,
  design = ~ Condition
)

# Exécution de l'analyse différentielle
dds <- DESeq(dds)
res <- results(dds)

# Conversion des résultats en dataframe pour ggplot2
res_df <- as.data.frame(res)

# Ajouter les noms des gènes depuis le fichier de comptage
res_df$gene_name <- counts_data$Gene_name[match(rownames(res_df), rownames(counts_data))]

# Ajouter une colonne pour identifier les gènes significatifs
res_df$point_type <- ifelse(res_df$padj < 0.05 & !is.na(res_df$padj), "Significant", "Non-significant")

# Liste des gènes à annoter
genes_to_annotate <- c("infA", "infB", "infC", "pth", "frr", "tsf")

# Filtrer les gènes pour l'annotation
annot_data <- res_df[res_df$gene_name %in% genes_to_annotate, ]

# Liste des gènes correspondants aux AA-tRNA-synthétases
trna_genes <- res_df[grep("S$", res_df$gene_name), ]

# Ajouter une colonne 'Annotation' pour identifier les AA-tRNA-synthétases
res_df$Annotation <- ifelse(res_df$gene_name %in% trna_genes$gene_name, "AA-tRNA-synthétases", "Autres")

# Création du ggplot
plot <- ggplot(res_df, aes(x = log2(baseMean + 1), y = log2FoldChange, color = point_type)) +
  geom_point(alpha = 1, size = 2) +  
  geom_point(data = trna_genes, aes(x = log2(baseMean + 1), y = log2FoldChange, color = "AA-tRNA-synthétases"), 
             size = 2, shape = 1, stroke = 1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth =1) +  
  geom_text(
    data = res_df[res_df$gene_name %in% genes_to_annotate, ],  
    aes(label = gene_name),
    hjust = -0.5, vjust = -1.5, size = 3, color = "black"  
  ) + 
  geom_segment(
    data = annot_data,
    aes(
      x = log2(baseMean), xend = log2(baseMean) + 0.5, 
      y = log2FoldChange, yend = log2FoldChange + 0.3 
    ),
    color = "black"
  ) +
  labs(
    title = "MA-plot of genes related to translation",
    x = "Log2 base Mean",
    y = "Log2 Fold Change",
    color = NULL  
  ) +
  scale_color_manual(values = c("Non-significant" = "darkgrey", "Significant" = "red", 
                                "AA-tRNA-synthétases" = "black")) + 
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 2)) +  # Limites et ticks de l'axe X
  scale_y_continuous(limits = c(-6, 5), breaks = seq(-6, 5, by = 1)) +  # Limites et ticks de l'axe Y
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth  = 1),  # Ajouter un encadré
    legend.title = element_blank() 
  )
ggsave("plot_translation.png", plot = plot)