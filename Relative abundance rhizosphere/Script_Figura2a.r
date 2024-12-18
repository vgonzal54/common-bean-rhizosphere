# This script processes metagenomic data to analyze and visualize the relative abundance of bacterial genera
# in different sample types. It integrates filters for bacterial taxa, aggregates data
# at the genus level, and generates bar plots of relative abundance percentages.

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting

# Import BIOM data
file_biom <- ("17_bean_metasamples.biom" )
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) 

# Read metadata file
metdata <- read.table("MetaData_17.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Adjust the order of levels in "Sample_Type" factor
levels(metdata$Sample_Type)  # View factor levels
metdata$Sample_Type <- factor(metdata$Sample_Type, levels = c('MetaG-PAA','MetaG-PAD','MetaG-PAI','MetaG-PNA','MetaG-PND','MetaG-PNI',
'MetaG-AA','MetaG-AB','MetaG-AD','MetaG-AE','MetaG-AH','MetaG-AI','MetaG-NA','MetaG-ND',
'MetaG-NE','MetaG-NH','MetaG-NI')) 
levels(metdata$Sample_Type)  # Verify factor levels

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
merged_metagenomes # Display summary of the combined object

# Rename columns in the taxonomic table
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
merged_metagenomes_Bacteria <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria") 
merged_metagenomes_Bacteria

# Aggregate data at the Genus level
Bacteria_Genus <- tax_glom(merged_metagenomes_Bacteria, taxrank = 'Genus')

# Transform abundance values to percentages
percentages_Genus <-  transform_sample_counts(Bacteria_Genus, function(x) x*100 / sum(x) )

# Convert phyloseq data to a data frame
Tabla_Genus <- psmelt(percentages_Genus) # Melt data to long format
unique(Tabla_Genus$Genus) # View unique Genus names

# Replace genera with low abundance (<1%) with a common label
Tabla_Genus$Genus <- as.character(Tabla_Genus$Genus)  # Ensure Genus column is a character vector
Tabla_Genus$Genus[Tabla_Genus$Abundance < 1] <- "Abundance ( < 1%) " # Replace low abundance genera
unique(Tabla_Genus$Genus) # Verify updated Genus names

# Define custom color palette
palette <- c("#E6E6E6","#CC9B7A","#497E00","#BCFFBC","#8F7700","#003C67","#7AA6DC",
"#CD534C","#868686","#EFC000","#0073C2","#D62728","#492900","#17BECF",
"#FFBB78","#FF5A5A","#BDBDBD","#A73030","#8F7700","#FFF1FF","#E7C9C6",
"#ABCD72","#FA9107","#D33C32","#F100F1","#FFFF2E","#C6DBEF","#14FFB1",
"#f46716","#4A6990","#6f04f9","#D5D5FF","#F4B7C0","#C46487","#C3BC3F",
"#835B82","#027B8E","#8DBFA8","#27c657")

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Generate bar plot of relative abundances
genus_abundance <- ggplot(data=Tabla_Genus, aes(x=Sample_Type , y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  labs(title = "",
        x = "",
        y = "Relative Abundance (%)") + theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_manual(name = "Genus", values = c(palette))

# Save plot as a PDF
pdf(file="PintoSaltillo_Relative_abundance_genus.pdf", 
    family= "Times", width = 10 , height = 7 ,pointsize = 10 , paper = "USr") #pointsize = 10)
print(genus_abundance)
dev.off()


## Further modifications to the PDF file were done using Adobe Illustrator


