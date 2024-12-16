setwd("D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Bray_Curtis_Var_PS/R_Data/")


# This script calculates and plot the betha diversity of microbial communities
# in different common bean accessions. 
 #Beta diversity mesures the variaton in community composition between samlpes. 


library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("vegan") # For ecological analysis

# Import BIOM data
file_biom <- ("All_PintoSalt_Accecion.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays importes data

# Read metadata file
metdata <- read.table("Matadata_all_PintoSAlt_Acc_Edit.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Adjust the order of levels in "Type_sample" factor
levels(metdata$Type_sample) # View factor levels
metdata$Type_sample <- factor(metdata$Type_sample, levels = c("MetaG-PAA","MetaG-PAD","MetaG-PAI","MetaG-PNA","MetaG-PND","MetaG-PNI","VP1-2","VP3-4","VP5-6","MetaG-AA","MetaG-AB","MetaG-AD","MetaG-AE",
"MetaG-AH","MetaG-AI","MetaG-NA","MetaG-ND","MetaG-NE","MetaG-NH","MetaG-NI","Bayo_V22","Bayo_V47","Flor_de_Mayo_V43","Negro_Criollo_V10","Negro_Criollo_V2",
"Negro_Criollo_V4","Negro_San_Luis_V26","Negro_V45","Peruano_V16"))
levels(metdata$Type_sample) # Verify factor levels

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
print(merged_metagenomes) # Display summary of the combined object

# Rename columns in the taxonomic table
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
merged_metagenomes_Bacteria <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria") 
merged_metagenomes_Bacteria

# Claculate beta diversity using Bray-Curtis distance 
BC <- distance(merged_metagenomes_Bacteria, method = "bray")

# Perform Principal Coordinates Analysis (PCoA) based on Bray-Curtis distances
ordBC = ordinate(merged_metagenomes_Bacteria, method = "PCoA", distance = BC)

# Export beta diversity  results to a CSV file
braycurtis.pcoa.export <- as.data.frame(ordBC$vectors, row.names=NULL, optional=FALSE, cut.names=FALSE, col.names=names(ordBC$vectors), fix.empty.names=TRUE, stringsAsFactors=default.stringsAsFactors())
write.csv(braycurtis.pcoa.export, file="D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Bray_Curtis_Var_PS/R_Results/braycurtis_29Samples.csv")

# Set default theme for plots
theme_set(theme_bw())

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Plot betha diversity (ordination results)
Plot <- plot_ordination(merged_metagenomes_Bacteria, ordBC, color="Bean_Type",shape="Bean_Type") +
ggtitle("PCoA: Bray-Curtis") +
theme(aspect.ratio=1) +
geom_point(alpha = 0.7 ,size=3) +
scale_colour_manual( values = c("orange","green","blue","purple","red"))

# Save plot to PDF
pdf("D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Bray_Curtis_Var_PS/R_Results/betha_diversity_palette1.pdf",         # Nombre del archivo
    width = 10, height = 19, # Ancho y alto en pulgadas
    bg = "white",          # Color de fondo
    paper = "letter",family= "Times")
print(Plot)
dev.off()


# Statistical analysis using PERMANOVA (adonis2)
adonis_results <- adonis2(BC ~sample_data(merged_metagenomes_Bacteria)$Origin_sample, 
    permutations = 999, method = "bray", p.adjst.m="holm")

print(adonis_results) # Display results