# Create a Venn diagram focusing on the genera enriched with a log2FoldChange greater than or equal to 2 and a 
# p-value less than 0.05. The objective is to analyze the effect of the three cultivars on the soil microbiota.

library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("ggrepel") # For improved label visualization
library("EnhancedVolcano") # For creating volcano plots
library("DESeq2") # For differential abundance analysis
library("ggvenn") # For creating Venn diagrams
library("gplots") # For Venn diagrams visualizations

# --- Analysis of Pinto Saltillo samples ---

# Import BIOM data
PS_file_biom <- ("Agricola_NonAgricola_12Sam.biom")
PS_data <- import_biom(PS_file_biom, parseFunction=parse_taxonomy_default)
print(PS_data) # Display imported data

# Read PS_metdata file
PS_metdata <- read.table("MetaData_Agricola_NonAgricola_12Sam.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(PS_metdata) # Preview PS_metdata structure

# Combine OTU table, taxonomic table, and PS_metdata into a single phyloseq object
PS_merged_metagenomes <- phyloseq(otu_table(PS_data), tax_table(PS_data), sample_data(PS_metdata)) 
print(PS_merged_metagenomes) # Display summary of the combined object

# Rename columns in the taxonomic table
PS_merged_metagenomes@tax_table@.Data <- substring(PS_merged_metagenomes@tax_table@.Data, 4)
colnames(PS_merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
PS_merged_metagenomes_Bacteria <- subset_taxa(PS_merged_metagenomes, Kingdom == "Bacteria") 
print(PS_merged_metagenomes_Bacteria) 
 
# Select data at the Genus level
PS_Genus <-  tax_glom(PS_merged_metagenomes_Bacteria, "Genus")
print(PS_Genus) # Display aggregated data

# Convert phyloseq data to DESeq2 format for differential abundance analysis
PS_Deq_B = phyloseq_to_deseq2(PS_Genus, ~ Sowing) # Set design formula for DESeq2
print(PS_Deq_B) # Display summary of the DESeq2 format

# Perform differential abundance analysis
PS_Deq_B = DESeq(PS_Deq_B)

# Extract results with specific contrast and significance cutoff
PS_resB = results(PS_Deq_B, contrast=c("Sowing", "After-sowing", "Before-sowing")  , alpha= 0.05) 
             
# Adjust significance level
PS_resB = PS_resB[order(PS_resB$padj),]

PS_tabla = cbind(as(PS_resB, "data.frame"), as(tax_table(PS_Genus)[rownames(PS_resB), ], "matrix"))

PS_table_LOG <- subset(PS_tabla, log2FoldChange >= 2)

PS_table_pvalue <- subset(PS_table_LOG, pvalue <=0.05)

PS_GENUS <- PS_table_pvalue[,12]


# --- Repeat the analysis for Bayo and Black bean varieties ---
# Similar steps are performed for "Bayo_Rhizo_Soil.biom" and "Negro_Rhizo_Soil.biom", 
# yielding BAYO_GENUS and BLACK_GENUS for the Bayo and Black bean varieties, respectively.

# --- Analysis of Bayo_bean samples ---


bayo_file_biom <- ("Bayo_Rhizo_Soil.biom")
bayo_data <- import_biom(bayo_file_biom, parseFunction=parse_taxonomy_default)
bayo_data



bayo_metdata <- read.table("MetaData_Bayo_Rhizo_Soil.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(bayo_metdata)

bayo_metdata

bayo_merged_metagenomes <- phyloseq(otu_table(bayo_data), tax_table(bayo_data), sample_data(bayo_metdata)) 

bayo_merged_metagenomes


bayo_merged_metagenomes@tax_table@.Data <- substring(bayo_merged_metagenomes@tax_table@.Data, 4)
colnames(bayo_merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

bayo_merged_metagenomes_Bacteria <- subset_taxa(bayo_merged_metagenomes, Kingdom == "Bacteria") 

bayo_merged_metagenomes_Bacteria 



bayo_Genus <-  tax_glom(bayo_merged_metagenomes_Bacteria, "Genus")

bayo_Genus


bayo_Deq_B_Genus = phyloseq_to_deseq2(bayo_Genus, ~ Origin_sample)

bayo_Deq_B_Genus = DESeq(bayo_Deq_B_Genus)




bayo_resB_Genus = results(bayo_Deq_B_Genus, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05) 
bayo_resB_Genus = bayo_resB_Genus[order(bayo_resB_Genus$padj),]


bayo_tabla = cbind(as(bayo_resB_Genus, "data.frame"), as(tax_table(bayo_Genus)[rownames(bayo_resB_Genus), ], "matrix"))

bayo_table_LOG <- subset(bayo_tabla, log2FoldChange >= 2)

bayo_table_pvalue <- subset(bayo_table_LOG, pvalue <=0.05)


BAYO_GENUS <- bayo_table_pvalue[,12]


# --- Analysis of Black_bean samples ---



black_file_biom <- ("Negro_Rhizo_Soil.biom")
black_data <- import_biom(black_file_biom, parseFunction=parse_taxonomy_default)
black_data


black_metdata <- read.table("MetaData_Negro_Rhizo_Soil.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(black_metdata)
black_metdata

black_merged_metagenomes <- phyloseq(otu_table(black_data), tax_table(black_data), sample_data(black_metdata)) 
black_merged_metagenomes



black_merged_metagenomes@tax_table@.Data <- substring(black_merged_metagenomes@tax_table@.Data, 4)
colnames(black_merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

black_merged_metagenomes_Bacteria <- subset_taxa(black_merged_metagenomes, Kingdom == "Bacteria") 
black_merged_metagenomes_Bacteria 

black_Genus <-  tax_glom(black_merged_metagenomes_Bacteria, "Genus")
black_Genus


black_Deq_B_Genus = phyloseq_to_deseq2(black_Genus, ~ Origin_sample)
black_Deq_B_Genus = DESeq(black_Deq_B_Genus)



black_resB_Genus = results(black_Deq_B_Genus, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05) 
black_resB_Genus = black_resB_Genus[order(black_resB_Genus$padj),]

black_tabla = cbind(as(black_resB_Genus, "data.frame"), as(tax_table(black_merged_metagenomes_Bacteria)[rownames(black_resB_Genus), ], "matrix"))

black_table_LOG <- subset(black_tabla, log2FoldChange >= 2)

black_table_pvalue <- subset(black_table_LOG, pvalue <=0.05)

BLACK_GENUS <- black_table_pvalue[,12]


# --- Create a Venn diagram of shared genera ---

Venn_Up <- list(Pinto_Saltillo =sample(PS_GENUS), Bayo_bean = sample(BAYO_GENUS), Black_bean=sample(BLACK_GENUS))

ggvenn(Venn_Up)
# Save the Venn diagram as a PDF
pdf("D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Volcano_Plot/Volcano_Plot_modificado/Venn_Diagram/R_Result/Venn_Diagram_Genus_LFC2_pvalue005.pdf",         # Nombre del archivo
    width = 8, height = 7, # Ancho y alto en pulgadas
    bg = "white",          # Color de fondo
    paper = "A4")   
ggvenn(Venn_Up)

dev.off()

# Extract intersections and their lengths
ItemsList <- venn(Venn_Up, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)

Intersection_content <- attributes(ItemsList)$intersections