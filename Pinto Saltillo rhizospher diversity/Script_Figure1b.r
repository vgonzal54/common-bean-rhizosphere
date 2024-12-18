# This script calculates and plot the betha diversity of microbial communities in envviromental samples using taxonomic data. 
# Beta diversity mesures the variaton in community composition between samlpes. 

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("vegan") # For ecological analysis

# Import BIOM data
file_biom <- ("17_bean_metasamples.biom" )
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)

# Read metadata file
metdata <- read.table("MetaData_17.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Adjust the order of levels in "Soil_Stage_Seeding" factor
levels(metdata$Soil_Stage_Seeding)
metdata$Soil_Stage_Seeding <- factor(metdata$Soil_Stage_Seeding, levels = c('A_Soil','A_Rhizo', 'N_Soil', 'N_Rhizo'))
levels(metdata$Soil_Stage_Seeding)

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
merged_metagenomes

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
write.csv(braycurtis.pcoa.export, file="Bray_Curtis_data.txt")

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Plot beta diversity (ordination results)
beta_plot <- plot_ordination(merged_metagenomes_Bacteria, ordBC, shape= "Soil_Stage_Seeding", color="Soil_Stage_Seeding") +
   scale_shape_manual(values=c(4,15,3,18)) +
   scale_color_manual(values=c('#D39200','#619CFF','#FF0000','#9A6A00')) +
   theme(aspect.ratio=1) + geom_point(size = 3) + 
   theme_bw() + theme(text=element_text(colour="black", size=12), 
   axis.text.x=element_text(colour="black", size = 10), 
   axis.text.y=element_text(colour="black", size = 10))

# Save plot to PDF
pdf(file= "Beta_diversity.pdf",
family= "Times", width = 8 , height = 6 , paper= "letter")
print(beta_plot)
dev.off()

# Statistical analysis using PERMANOVA (adonis2)

levels(metdata$Sowing)
#[1] "Before-sowing" "After-sowing" 
adonis2(BC~sample_data(merged_metagenomes_Bacteria)$Sowing, permutations = 999, method = "bray", p.adjst.m="holm")
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = BC ~ sample_data(merged_metagenomes_Bacteria)$Sowing, permutations = 999, method = "bray", p.adjst.m = "holm")
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     1   3.0567 0.81476 65.975  0.001 ***
#Residual 15   0.6950 0.18524                  
#Total    16   3.7517 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Further modifications to the PDF file were done using Adobe Illustrator


