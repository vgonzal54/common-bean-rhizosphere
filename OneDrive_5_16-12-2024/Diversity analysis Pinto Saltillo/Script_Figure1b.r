# This script calculates and plot the betha diversity of microbial communities
# in envviromental samples using taxonomic data. 
#Beta diversity mesures the variaton in community composition between samlpes. 

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("vegan") # For ecological analysis

# Import BIOM data
file_biom <- ("17_bean_metasamples.biom" )
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays importes data

# Read metadata file
metdata <- read.table("MetaData_Mod.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Adjust the order of levels in "Sowing" factor
levels(metdata$Sowing)  # View factor levels
metdata$Sowing <- factor(metdata$Sowing, levels = c('Before-sowing','After-sowing')) # Con esta línea se modifica el orden
levels(metdata$Sowing) # Verify factor levels

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
write.csv(braycurtis.pcoa.export, file="D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Pinto_Saltillo_Graph/R_Result/Bray_Curtis_default.txt")


# Set default theme for plots
theme_set(theme_bw())

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Plot betha diversity (ordination results)
beta_plot <- plot_ordination(merged_metagenomes_Bacteria, ordBC, shape= "Soil_Stage_Seeding", color="Soil_Stage_Seeding") +
scale_shape_manual(values=c(4, 15,3,18 )) +
scale_color_manual(values=c("#CD950C","#009ACD", "#CD2626", "#8B5A00")) +
theme(aspect.ratio=1) + geom_point(size = 3) + 
theme_bw() + theme(text=element_text(colour="black", size=12), 
axis.text.x=element_text(colour="black", size = 10), 
axis.text.y=element_text(colour="black", size = 10))

# Save plot to PDF
pdf(file= "D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Pinto_Saltillo_Graph/R_Result/Beta_diversidad.pdf",
family= "Times", width = 8 , height = 6 , paper= "letter")
print(beta_plot)
dev.off()

# Statistical analysis using PERMANOVA (adonis2)
adonis_results <- adonis2(BC ~sample_data(merged_metagenomes_Bacteria)$Sowing, 
    permutations = 999, method = "bray", p.adjst.m="holm")

print(adonis_results) # Display results