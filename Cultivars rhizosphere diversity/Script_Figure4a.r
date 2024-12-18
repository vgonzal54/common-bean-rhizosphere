# This script processes metagenomic data to analyze and visualize the relative abundance of bacterial genera
# in different sample types. It integrates filters for bacterial taxa, aggregates data
# at the genus level, and generates bar plots of relative abundance percentages.

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting

# Import BIOM data
file_biom <- ("acecciones.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays imported data

# Read metadata file
metdata <- read.table("Metadata_aceccions.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata)# Preview metadata structure

# Adjust the order of levels in "Sample_Type" factor
levels(metdata$Coloquial_names) # View factor levels
metdata$Coloquial_names <- factor(metdata$Coloquial_names, levels = c('VP1-2', 'VP3-4', 'VP5-6', 'Black_V10', 'Black_V26', 'Black_V2', 
                        'Black_V45', 'Black_V4','Bayo_V16', 'Bayo_V22', 'Bayo_V43', 'Bayo_V47'))
levels(metdata$Coloquial_names) # Verify factor levels

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
merged_metagenomes # Display summary of the combined object

# Rename columns in the taxonomic table
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
merged_metagenomes_Bacteria  <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
merged_metagenomes_Bacteria

# Aggregate data at the Genus level
Bacteria_Genus <- tax_glom(merged_metagenomes_Bacteria, taxrank = 'Genus')

# Transform abundance values to percentages
percentages_Genus <-  transform_sample_counts(Bacteria_Genus, function(x) x*100 / sum(x) )

# Convert phyloseq data to a data frame
Tabla_Genus <- psmelt(percentages_Genus) # Melt data to long format

# Replace genera with low abundance (<1%) with a common label
Tabla_Genus$Genus <- as.character(Tabla_Genus$Genus)  # Ensure Genus column is a character vector
Tabla_Genus$Genus[Tabla_Genus$Abundance < 1] <- "Abundance ( < 1%) " # Replace low abundance genera
unique(Tabla_Genus$Genus) # Verify updated Genus names

# Define custom color palette
palette <- c("#E6E6E6","#9EDAE5","#CC9B7A","#497E00","#BCFFBC","#8F7700","#003C67",
"#868686","#D092A7","#8971E1","#EFC000","#492900","#17BECF","#FFBB78",
"#BDBDBD","#D3BA68","#FFF1FF","#E7C9C6","#D33C32","#F100F1","#FFFF2E",
"#C6DBEF","#f46716","#4A6990","#FF0080","#6f04f9","#F4B7C0","#FFFFBF",
"#CCEBC5","#C3BC3F","#835B82","#027B8E","#8DBFA8","#27c657")

# Set default theme for plots
theme_set(theme_bw())

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Generate bar plot of relative abundances
Varieties_plot <- ggplot(data=Tabla_Genus, aes(x=Coloquial_names , y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  labs(title = "",
        x = "",
        y = "Relative Abundance (%)") + theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
  scale_fill_manual(name = "Genus", values = c(palette))

# Save plot as a PDF
pdf(file="Var_Relative_abundance_genus.pdf", 
    family= "Times", width = 10 , height = 7 ,pointsize = 10 , paper = "USr") 
print(Varieties_plot)
dev.off()


## Further modifications to the PDF file were done using Adobe Illustrator


