# This script calculates and plot the betha diversity of microbial communities in different common bean accessions. 
# Beta diversity mesures the variaton in community composition between samlpes. 

library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("vegan") # For ecological analysis
library("pairwiseAdonis") # For permanova pairwise analysis

# Import BIOM data
file_biom <- ("All_PintoSalt_Accecion.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays imported data

# Read metadata file
metdata <- read.table("Matadata_all_PintoSAlt_Acc.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Adjust the order of levels in "Type" factor
levels(metdata$Type)
metdata$Type <- factor(metdata$Type, levels = c('Bayo_bean', 'Black_bean', 'Pinto_Saltillo', 'Bayo_and_Black_soil', 'Pinto_Saltillo_soil'))
levels(metdata$Type)

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
write.csv(braycurtis.pcoa.export, file="bray-curtis_29Samples.csv")

# Set default theme for plots
theme_set(theme_bw())

# Define custom discrete fill scale
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

# Plot beta diversity (ordination results)
Plot <- plot_ordination(merged_metagenomes_Bacteria, ordBC, color="Type", shape="Type") + 
   ggtitle("PCoA: Bray-Curtis") + 
   theme(aspect.ratio=1) + geom_point(size = 3) + 
   scale_color_manual(values=c('#F8766D','#00BA38','#619CFF','#DB72FB','#D39200', '#FF61C3')) + 
   scale_shape_manual(values=c(16,17,15,18,4,25)) + theme_bw() + 
   theme(text=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size = 10), 
    axis.text.y=element_text(colour="black", size = 10))+ theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

# Save plot to PDF
pdf("beta_diversity_29samples.pdf",
    width = 10, height = 19, bg = "white", paper = "letter", family= "Times")
print(Plot)
dev.off()


# Statistical analysis using PERMANOVA (adonis2)

adonis2(BC~sample_data(merged_metagenomes_Bacteria)$Type, permutations = 999, method = "bray", p.adjst.m="holm")
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = BC ~ sample_data(merged_metagenomes_Bacteria)$Type, permutations = 999, method = "bray", p.adjst.m = "holm")
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     4   4.7945 0.77542 20.716  0.001 ***
#Residual 24   1.3886 0.22458                  
#Total    28   6.1832 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.adonis(BC,factors=sample_data(merged_metagenomes_Bacteria)$Type)

dispersion <- betadisper(d=BC, group=sample_data(merged_metagenomes_Bacteria)$Type, type="centroid")
dispersion
permutest(dispersion, permutations=999, pairwise=TRUE)


## Further modifications to the PDF file were done using Adobe Illustrator


