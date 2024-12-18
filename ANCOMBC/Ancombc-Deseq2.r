library(vegan)
library(ggplot2)
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(reshape2)
library(cowplot)
library(DESeq2)
library(tibble)


## Differential abundance analysis of Pinto Saltillo vs Pinto Saltillo Soil

archivo_biom_12 <- ("Agricola_NonAgricola_12Sam.biom")
data12 <- import_biom(archivo_biom_12, parseFunction=parse_taxonomy_default)
metdata12 <- read.table("MetaData_Agricola_NonAgricola_12Sam.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)

merged12 <- phyloseq(otu_table(data12), tax_table(data12), sample_data(metdata12)) 
merged12@tax_table@.Data <- substring(merged12@tax_table@.Data, 4)
colnames(merged12@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

merged12_Bacteria <- subset_taxa(merged12, Kingdom == "Bacteria") 
Genus12 <-  tax_glom(merged12_Bacteria, "Genus")
Genus12
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1474 taxa and 12 samples ]
#sample_data() Sample Data:       [ 12 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1474 taxa by 7 taxonomic ranks ]

#View(otu_table(Genus12))
#View(tax_table(Genus12))

#reads_num <- cbind(data.frame(otu_table(Genus12)), data.frame(tax_table(Genus12)))
#write.table(reads_num, file = "Pinto_saltillo_12_Abundancia_genus.txt", row.names=TRUE, sep="\t")

Deq_B = phyloseq_to_deseq2(Genus12, ~ Sowing)
Deq_B
Deq_B = DESeq(Deq_B)
resB = results(Deq_B, contrast=c("Sowing", "After-sowing", "Before-sowing")  , alpha= 0.05) # cooksCutoff=FALSE)
resB = resB[order(resB$padj),]
tabla = cbind(as(resB, "data.frame"), as(tax_table(Genus12)[rownames(resB), ], "matrix"))
#write.table(tabla, file ="Pinto_saltillo_12_LFC.tsv" , row.names=TRUE, sep="\t")

# To rename duplicate taxonomic genera
taxdata12 <- as.data.frame(tax_table(Genus12))  # objeto phyloseq de Gris con 12 muestras Pinto Saltillo
duplicatedgen12 <- duplicated(taxdata12$Genus)
#taxdata12[duplicatedgen12, ]
taxdata12$Genus <- make.unique(taxdata12$Genus, sep = "_")
tax_table(Genus12) <- tax_table(as.matrix(taxdata12))

# To filter (remove) taxa with more than 90% zeros
otudata12 <- as.data.frame(otu_table(Genus12))
sampledata12 <- as.data.frame(sample_data(Genus12))
taxazeroperc12 <- rowSums(otudata12 == 0) / ncol(otudata12)
#zeropercent12 <- cbind(taxazeroperc12, taxdata12)
#zeropercent12 <- tibble::rownames_to_column(zeropercent12, var = "IDs")
#write.table(zeropercent12, file = "Pinto_saltillo_12_otu_zero_percent.txt", row.names=FALSE, sep="\t")
otu_filtered12 <- otudata12[taxazeroperc12 < 0.90, ]
otu_table(Genus12) <- otu_table(as.matrix(otu_filtered12), taxa_are_rows = TRUE)
Genus12
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1469 taxa and 12 samples ]
#sample_data() Sample Data:       [ 12 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1469 taxa by 7 taxonomic ranks ]

Deq_Bf = phyloseq_to_deseq2(Genus12, ~ Sowing)
Deq_Bf = DESeq(Deq_Bf)
resBf = results(Deq_Bf, contrast=c("Sowing", "After-sowing", "Before-sowing")  , alpha= 0.05) # cooksCutoff=FALSE)
resBf = resBf[order(resBf$padj),]
tablaf = cbind(as(resBf, "data.frame"), as(tax_table(Genus12)[rownames(resBf), ], "matrix"))
#write.table(tablaf, file ="Pinto_saltillo_12_LFC_filtrado.tsv" , row.names=TRUE, sep="\t")

filtered_deseq_pinto12 <- tablaf %>%
  filter(
    (padj < 0.05 & (log2FoldChange >= 2 | log2FoldChange <= -2)) )
#filtered_deseq_pinto12 <- tibble::rownames_to_column(filtered_deseq_pinto12, var = "IDs")
#View(filtered_deseq_pinto12)
#write.table(filtered_deseq_pinto12, file = "Pinto_saltillo_12_LFC_filtrado-2.tsv", row.names=FALSE, sep="\t", quote=FALSE)


# ANCOM-BC with filtered phyloseq object

# Change the reference level to "Before-sowing"
Genus12@sam_data$Sowing <- relevel(factor(Genus12@sam_data$Sowing), ref = "Before-sowing")

# Run ANCOM-BC with "Before-sowing" as reference
ancombc_res_before <- ancombc( data = Genus12, formula = "Sowing", p_adj_method = "holm", lib_cut = 1000, group = "Sowing")

# Extract results
ancombc_pinto12 <- ancombc_res_before$res
# Extracting taxonomy data from `phyloseq`
taxonomy12 <- as.data.frame(tax_table(Genus12))
View(taxonomy12)
# Make sure taxon IDs are in a column named "TaxonID"
taxonomy12$TaxonID <- rownames(taxonomy12)
# Convert `ancombc_pinto12` to a data frame
ancombc_df_pinto12 <- as.data.frame(ancombc_pinto12)
View(ancombc_df_pinto12)
# Rename taxon identifier column in `ancombc_df_pinto12`
colnames(ancombc_df_pinto12)[colnames(ancombc_df_pinto12) == "lfc.taxon"] <- "TaxonID"
View(ancombc_df_pinto12)
# Join `ancombc_df_pinto12` with `taxonomy12` by column "TaxonID"
ancombc_merged_pinto12 <- dplyr::left_join(ancombc_df_pinto12, taxonomy12[, c("TaxonID", "Genus")], by = "TaxonID")
View(ancombc_merged_pinto12)
#write.table(ancombc_merged_pinto12, file="Pinto_saltillo_12_ancombc.txt", row.names=FALSE, sep="\t", quote=FALSE)

filtered_ancombc_pinto12 <- ancombc_merged_pinto12 %>%
  filter(
    (q_val.SowingAfter.sowing < 0.05 & (lfc.SowingAfter.sowing >= 2 | lfc.SowingAfter.sowing <= -2)) )
#View(filtered_ancombc_pinto12)
#write.table(filtered_ancombc_pinto12, file = "Pinto_saltillo_12_ancombc_filtered.txt", row.names=FALSE, sep="\t", quote=FALSE)

# Filter only the columns of interest for comparisons
lfc_data12 <- ancombc_merged_pinto12 %>%
  select(Genus, lfc.SowingAfter.sowing, q_val.SowingAfter.sowing, se.SowingAfter.sowing)
# Rename columns
colnames(lfc_data12)[colnames(lfc_data12) == "lfc.SowingAfter.sowing"] <- "LFC"
colnames(lfc_data12)[colnames(lfc_data12) == "q_val.SowingAfter.sowing"] <- "q_val"
colnames(lfc_data12)[colnames(lfc_data12) == "se.SowingAfter.sowing"] <- "SE"
View(lfc_data12)

# Filter the data to include only genera with LFC >= 2 or LFC <= -2 and q_val < 0.05
lfc_filtered12 <- lfc_data12 %>%
  filter((LFC >= 2 | LFC <= -2), q_val < 0.05) %>%
  arrange(desc(LFC)) %>% mutate(Genus=factor(Genus,levels=unique(Genus)))
View(lfc_filtered12)

# Create the plot
if (nrow(lfc_filtered12) > 0) {
  plot <- ggplot(lfc_filtered12, aes(x = Genus, y = LFC, fill = LFC > 0)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = LFC - SE, ymax = LFC + SE), width = 0.1) +
    scale_fill_manual(values = c("red", "blue"), labels = c("Negative LFC", "Positive LFC")) +
    theme_bw() +
    labs(
      x = "Genus",
      y = "Log Fold Change (LFC)",
      fill = "",
      title = "ANCOM-BC: Pinto Saltillo vs Soil Pinto Saltillo (reference)"
    ) +
    scale_y_continuous(breaks=seq(-10,10,2)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      text = element_text(size = 10),
      legend.position = "top"
    ) + coord_flip() +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave("Pinto_saltillo_12_ancombc_filtrado.pdf", plot=plot, width=6, height=11, family="Times")
} else {
  message("There is no data that meets the filter criteria.")
}


## Differential abundance analysis of Black bean vs. bayo_and_black_soil

archivo_biomNegro <- ("Negro_Rhizo_Soil.biom")
dataNegro <- import_biom(archivo_biomNegro, parseFunction=parse_taxonomy_default)
metdataNegro <- read.table("MetaData_Negro_Rhizo_Soil.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
merged_metagenomesNegro <- phyloseq(otu_table(dataNegro), tax_table(dataNegro), sample_data(metdataNegro)) 
merged_metagenomesNegro@tax_table@.Data <- substring(merged_metagenomesNegro@tax_table@.Data, 4)
colnames(merged_metagenomesNegro@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
merged_metagenomes_BacNegro <- subset_taxa(merged_metagenomesNegro, Kingdom == "Bacteria") 
GenusNegro <-  tax_glom(merged_metagenomes_BacNegro, "Genus")
GenusNegro
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1464 taxa and 8 samples ]
#sample_data() Sample Data:       [ 8 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1464 taxa by 7 taxonomic ranks ]

#readsNegro <- cbind(data.frame(otu_table(GenusNegro)), data.frame(tax_table(GenusNegro)))
#write.table(readsNegro, file = "Negro_Abundancia_genus.txt", row.names=TRUE, sep="\t")

DeqB_Negro = phyloseq_to_deseq2(GenusNegro, ~ Origin_sample)
DeqB_Negro = DESeq(DeqB_Negro)
resB_Negro = results(DeqB_Negro, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05)       # cooksCutoff=FALSE)
resB_Negro = resB_Negro[order(resB_Negro$padj),]
tablaNegro = cbind(as(resB_Negro, "data.frame"), as(tax_table(GenusNegro)[rownames(resB_Negro), ], "matrix"))
tablaNegro <- tibble::rownames_to_column(tablaNegro, var = "IDs")
#write.table(tablaNegro, file ="Negro_LFC.tsv" , row.names=FALSE, sep="\t")

taxdataNegro <- as.data.frame(tax_table(GenusNegro))
duplicatedgenNegro <- duplicated(taxdataNegro$Genus)
#taxdataNegro[duplicatedgenNegro, ]
taxdataNegro$Genus <- make.unique(taxdataNegro$Genus, sep = "_")
tax_table(GenusNegro) <- tax_table(as.matrix(taxdataNegro))

otudataNegro <- as.data.frame(otu_table(GenusNegro))
sampledataNegro <- as.data.frame(sample_data(GenusNegro))
taxazeropercNegro <- rowSums(otudataNegro == 0) / ncol(otudataNegro)
#zeropercentNegro <- cbind(taxazeropercNegro, taxdataNegro)
#zeropercentNegro <- tibble::rownames_to_column(zeropercentNegro, var = "IDs")
#write.table(zeropercentNegro, file = "Negro_otu_zero_percent.txt", row.names=FALSE, sep="\t")
otu_filteredNegro <- otudataNegro[taxazeropercNegro < 0.90, ]
otu_table(GenusNegro) <- otu_table(as.matrix(otu_filteredNegro), taxa_are_rows = TRUE)
GenusNegro
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1464 taxa and 8 samples ]
#sample_data() Sample Data:       [ 8 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1464 taxa by 7 taxonomic ranks ]

DeqBNegro = phyloseq_to_deseq2(GenusNegro, ~ Origin_sample)
DeqBNegro = DESeq(DeqBNegro)
resBNegro = results(DeqBNegro, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05)       # cooksCutoff=FALSE)
resBNegro = resB_Negro[order(resBNegro$padj),]
tablNegro = cbind(as(resBNegro, "data.frame"), as(tax_table(GenusNegro)[rownames(resBNegro), ], "matrix"))
tablNegro <- tibble::rownames_to_column(tablNegro, var = "IDs")
#write.table(tablNegro, file ="Negro_LFC_filtrado.tsv" , row.names=FALSE, sep="\t")

filtered_deseq_Negro <- tablNegro %>%
  filter(
    (padj < 0.05 & (log2FoldChange >= 2 | log2FoldChange <= -2)) )
#filtered_deseq_Negro <- tibble::rownames_to_column(filtered_deseq_Negro, var = "IDs")
#View(filtered_deseq_Negro)
#write.table(filtered_deseq_Negro, file = "Negro_LFC_filtrado-2.tsv", row.names=FALSE, sep="\t", quote=FALSE)

# ANCOM-BC with filtered phyloseq object

# Change the reference level to "Soil"
GenusNegro@sam_data$Origin_sample <- relevel(factor(GenusNegro@sam_data$Origin_sample), ref = "Soil")

# Run ANCOM-BC with "Soil" as reference
ancombcNegro <- ancombc( data = GenusNegro, formula = "Origin_sample", p_adj_method = "holm", lib_cut = 1000, group = "Origin_sample")

ancombc_Negro <- ancombcNegro$res
taxonomyNegro <- as.data.frame(tax_table(GenusNegro))
taxonomyNegro$TaxonID <- rownames(taxonomyNegro)
ancombc_df_Negro <- as.data.frame(ancombc_Negro)
colnames(ancombc_df_Negro)[colnames(ancombc_df_Negro) == "lfc.taxon"] <- "TaxonID"
ancombc_merged_Negro <- dplyr::left_join(ancombc_df_Negro, taxonomyNegro[, c("TaxonID", "Genus")], by = "TaxonID")
#View(ancombc_merged_Negro)
#write.table(ancombc_merged_Negro, file="Negro_ancombc.txt", row.names=FALSE, sep="\t", quote=FALSE)

filtered_ancombc_Negro <- ancombc_merged_Negro %>%
  filter(
    (q_val.Origin_sampleRhizosphere < 0.05 & (lfc.Origin_sampleRhizosphere >= 2 | lfc.Origin_sampleRhizosphere <= -2)) )
#View(filtered_ancombc_Negro)
#write.table(filtered_ancombc_Negro, file = "Negro_ancombc_filtered.txt", row.names=FALSE, sep="\t", quote=FALSE)

lfc_dataNegro <- ancombc_merged_Negro %>%
  select(Genus, lfc.Origin_sampleRhizosphere, q_val.Origin_sampleRhizosphere, se.Origin_sampleRhizosphere)
colnames(lfc_dataNegro)[colnames(lfc_dataNegro) == "lfc.Origin_sampleRhizosphere"] <- "LFC"
colnames(lfc_dataNegro)[colnames(lfc_dataNegro) == "q_val.Origin_sampleRhizosphere"] <- "q_val"
colnames(lfc_dataNegro)[colnames(lfc_dataNegro) == "se.Origin_sampleRhizosphere"] <- "SE"
#View(lfc_dataNegro)
lfc_filteredNegro <- lfc_dataNegro %>%
  filter((LFC >= 2 | LFC <= -2), q_val < 0.05) %>%
  arrange(desc(LFC)) %>% mutate(Genus=factor(Genus,levels=unique(Genus)))
#View(lfc_filteredNegro)

plotNegro <- ggplot(lfc_filteredNegro, aes(x = Genus, y = LFC, fill = LFC > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = LFC - SE, ymax = LFC + SE), width = 0.1) +
  scale_fill_manual(values = c("blue", "red"), labels = c("Positive LFC", "Negative LFC")) +
  theme_bw() +
  labs(
    x = "Genus",
    y = "Log Fold Change (LFC)",
    fill = "",
    title = "ANCOM-BC: Negro vs Soil other varieties (reference)"
  ) +
  scale_y_continuous(breaks=seq(-10,10,2)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size=10),
    text = element_text(size = 10),
    legend.position = "top"
  ) + coord_flip() +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave("Negro_ancombc_filtrado.pdf", plot=plotNegro, width=6, height=11, family="Times")


## Differential abundance analysis of Bayo bean vs. bayo_and_black_soil

archivo_biomBayo <- ("Bayo_Rhizo_Soil.biom")
dataBayo <- import_biom(archivo_biomBayo, parseFunction=parse_taxonomy_default)
metdataBayo <- read.table("MetaData_Bayo_Rhizo_Soil.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
merged_metagenomesBayo <- phyloseq(otu_table(dataBayo), tax_table(dataBayo), sample_data(metdataBayo)) 
merged_metagenomesBayo@tax_table@.Data <- substring(merged_metagenomesBayo@tax_table@.Data, 4)
colnames(merged_metagenomesBayo@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
merged_metagenomes_BacBayo <- subset_taxa(merged_metagenomesBayo, Kingdom == "Bacteria") 
GenusBayo <-  tax_glom(merged_metagenomes_BacBayo, "Genus")
GenusBayo
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1465 taxa and 7 samples ]
#sample_data() Sample Data:       [ 7 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1465 taxa by 7 taxonomic ranks ]

#readsBayo <- cbind(data.frame(otu_table(GenusBayo)), data.frame(tax_table(GenusBayo)))
#readsBayo <- tibble::rownames_to_column(readsBayo, var = "IDs")
#write.table(readsBayo, file = "Bayo_Abundancia_genus.txt", row.names=FALSE, sep="\t")

DeqB_Bayo = phyloseq_to_deseq2(GenusBayo, ~ Origin_sample)
DeqB_Bayo = DESeq(DeqB_Bayo)
resB_Bayo = results(DeqB_Bayo, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05)       # cooksCutoff=FALSE)
resB_Bayo = resB_Bayo[order(resB_Bayo$padj),]
tablaBayo = cbind(as(resB_Bayo, "data.frame"), as(tax_table(GenusBayo)[rownames(resB_Bayo), ], "matrix"))
tablaBayo <- tibble::rownames_to_column(tablaBayo, var = "IDs")
#write.table(tablaBayo, file ="Bayo_LFC.tsv" , row.names=FALSE, sep="\t")

taxdataBayo <- as.data.frame(tax_table(GenusBayo))  # objeto phyloseq de Gris con 7 muestras Bayo
duplicatedgenBayo <- duplicated(taxdataBayo$Genus)
#taxdataBayo[duplicatedgenBayo, ]
taxdataBayo$Genus <- make.unique(taxdataBayo$Genus, sep = "_")
tax_table(GenusBayo) <- tax_table(as.matrix(taxdataBayo))

otudataBayo <- as.data.frame(otu_table(GenusBayo))
sampledataBayo <- as.data.frame(sample_data(GenusBayo))
taxazeropercBayo <- rowSums(otudataBayo == 0) / ncol(otudataBayo)
#zeropercentBayo <- cbind(taxazeropercBayo, taxdataBayo)
#zeropercentBayo <- tibble::rownames_to_column(zeropercentBayo, var = "IDs")
#write.table(zeropercentBayo, file = "Bayo_otu_zero_percent.txt", row.names=FALSE, sep="\t")
otu_filteredBayo <- otudataBayo[taxazeropercBayo < 0.90, ]
otu_table(GenusBayo) <- otu_table(as.matrix(otu_filteredBayo), taxa_are_rows = TRUE)
GenusBayo
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1465 taxa and 7 samples ]
#sample_data() Sample Data:       [ 7 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1465 taxa by 7 taxonomic ranks ]

DeqBBayo = phyloseq_to_deseq2(GenusBayo, ~ Origin_sample)
DeqBBayo = DESeq(DeqBBayo)
resBBayo = results(DeqBBayo, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05)       # cooksCutoff=FALSE)
resBBayo = resB_Bayo[order(resBBayo$padj),]
tablBayo = cbind(as(resBBayo, "data.frame"), as(tax_table(GenusBayo)[rownames(resBBayo), ], "matrix"))
tablBayo <- tibble::rownames_to_column(tablBayo, var = "IDs")
#write.table(tablBayo, file ="Bayo_LFC_filtrado.tsv" , row.names=FALSE, sep="\t")

filtered_deseq_Bayo <- tablBayo %>%
  filter(
    (padj < 0.05 & (log2FoldChange >= 2 | log2FoldChange <= -2)) )
#View(filtered_deseq_Bayo)
#write.table(filtered_deseq_Bayo, file = "Bayo_LFC_filtrado-2.tsv", row.names=FALSE, sep="\t", quote=FALSE)

# ANCOM-BC with filtered phyloseq object

# Change the reference level to "Soil"
GenusBayo@sam_data$Origin_sample <- relevel(factor(GenusBayo@sam_data$Origin_sample), ref = "Soil")

# Run ANCOM-BC with "Soil" as reference
ancombcBayo <- ancombc( data = GenusBayo, formula = "Origin_sample", p_adj_method = "holm", lib_cut = 1000, group = "Origin_sample")

ancombc_Bayo <- ancombcBayo$res
taxonomyBayo <- as.data.frame(tax_table(GenusBayo))
taxonomyBayo$TaxonID <- rownames(taxonomyBayo)
ancombc_df_Bayo <- as.data.frame(ancombc_Bayo)
colnames(ancombc_df_Bayo)[colnames(ancombc_df_Bayo) == "lfc.taxon"] <- "TaxonID"
ancombc_merged_Bayo <- dplyr::left_join(ancombc_df_Bayo, taxonomyBayo[, c("TaxonID", "Genus")], by = "TaxonID")
#View(ancombc_merged_Bayo)
#write.table(ancombc_merged_Bayo, file="Bayo_ancombc.txt", row.names=FALSE, sep="\t", quote=FALSE)

filtered_ancombc_Bayo <- ancombc_merged_Bayo %>%
  filter(
    (q_val.Origin_sampleRhizosphere < 0.05 & (lfc.Origin_sampleRhizosphere >= 2 | lfc.Origin_sampleRhizosphere <= -2)) )
#View(filtered_ancombc_Bayo)
#write.table(filtered_ancombc_Bayo, file = "Bayo_ancombc_filtered.txt", row.names=FALSE, sep="\t", quote=FALSE)

lfc_dataBayo <- ancombc_merged_Bayo %>%
  select(Genus, lfc.Origin_sampleRhizosphere, q_val.Origin_sampleRhizosphere, se.Origin_sampleRhizosphere)
colnames(lfc_dataBayo)[colnames(lfc_dataBayo) == "lfc.Origin_sampleRhizosphere"] <- "LFC"
colnames(lfc_dataBayo)[colnames(lfc_dataBayo) == "q_val.Origin_sampleRhizosphere"] <- "q_val"
colnames(lfc_dataBayo)[colnames(lfc_dataBayo) == "se.Origin_sampleRhizosphere"] <- "SE"
lfc_filteredBayo <- lfc_dataBayo %>%
  filter((LFC >= 2 | LFC <= -2), q_val < 0.05) %>%
  arrange(desc(LFC)) %>% mutate(Genus=factor(Genus,levels=unique(Genus)))
#View(lfc_filteredBayo)

plotBayo <- ggplot(lfc_filteredBayo, aes(x = Genus, y = LFC, fill = LFC > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = LFC - SE, ymax = LFC + SE), width = 0.1) +
  scale_fill_manual(values = c("blue", "red"), labels = c("Positive LFC", "Negative LFC")) +
  theme_bw() +
  labs(
    x = "Genus",
    y = "Log Fold Change (LFC)",
    fill = "",
    title = "ANCOM-BC: Bayo vs Soil other varieties (reference)"
  ) +
  scale_y_continuous(breaks=seq(-10,10,2)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size=10),
    text = element_text(size = 10),
    legend.position = "top"
  ) + coord_flip() +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave("Bayo_ancombc_filtrado.pdf", plot=plotBayo, width=6, height=11, family="Times")


