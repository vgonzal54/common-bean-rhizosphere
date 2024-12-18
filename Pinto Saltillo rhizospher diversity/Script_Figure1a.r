# Alpha Diversity Analysis Script
# Calculate and visualize alpha diversity (Chao1, Shannon, InvSimpson) for metagenomic samples.


# Charge the libraries 
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("ggpubr")  # For statistical annotations on plots
library("gridExtra") # For combining multiple plots

# Import BIOM data  
file_biom <- ("17_bean_metasamples.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays imported data

# Read metadata file
metdata <- read.table("MetaData_17.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata)

# Adjust the order of factor levels for "Soil_stage_Seeding"
levels(metdata$Soil_Stage_Seeding) # View factor levels
metdata$Soil_Stage_Seeding <- factor(metdata$Soil_Stage_Seeding, levels = c('A_Soil','A_Rhizo', 'N_Soil', 'N_Rhizo'))
levels(metdata$Soil_Stage_Seeding) # Verify factor levels

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
print(merged_metagenomes) # Display summary of the combined object

# Rename columns in the taxonomic table 
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
merged_metagenomes_Bacteria <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria") 
merged_metagenomes_Bacteria

# Estimate alpha diversity 
richness <- estimate_richness(merged_metagenomes_Bacteria)

# Save alpha diversity results to a CSV file 
write.csv(richness, file="Alfa_diversity_Bacteria.csv")

# Define pairwise comparisons for statistical testing  
custom_ord <- c('A_Soil','A_Rhizo', 'N_Soil', 'N_Rhizo')
comparition <- list(c('A_Soil','A_Rhizo','N_Soil','N_Rhizo'))
compare <- list(c('N_Soil', 'N_Rhizo'), c('A_Soil', 'A_Rhizo'), c('A_Soil', 'N_Soil'), c('A_Rhizo', 'N_Rhizo'))

# Plot alpha diversity with statistical test 
palette <- c("#533D43", "#DABA22", "#AF1B97", "#258B08")
a <- plot_richness(merged_metagenomes_Bacteria, x='Soil_Stage_Seeding', measures=c("Chao1","Shannon","InvSimpson", nrow=1)) + 
   geom_boxplot(alpha=0.4, aes(color=Soil_Stage_Seeding)) + geom_point(alpha=0.6, aes(color=Soil_Stage_Seeding)) + 
   scale_fill_manual(values=palette) + scale_color_manual(values=palette, guide= "none") + theme_bw() + 
   scale_x_discrete(limits=c(custom_ord)) + theme(panel.background=element_blank(), strip.background=element_blank()) + 
   theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4), 
     plot.margin = margin(0.25, 0.25, 1, 0.25, "cm")) + theme(text=element_text(colour="black", size=10), 
     axis.text.x=element_text(colour="black", size=8), axis.text.y=element_text(colour="black", size=8)) + 
   theme(plot.margin = margin(1.5, 1.5, 1.5, 1, "cm")) + 
   stat_compare_means(method="wilcox.test", comparisons=compare, size=2)

pdf("Alfa_diversity_Bacteria_boxplot.pdf", family="Times", width=7, height=5)
a
dev.off()


## Further modifications to the PDF file were done using Adobe Illustrator


