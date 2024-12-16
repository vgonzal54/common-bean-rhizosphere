# Alpha Diversity Analysis Script
# Calculate and visualize alpha diversity (Chao1, Shannon, InvSimpson) for metagenomic samples.


# Charge the libraries 
library("phyloseq") 
library("ggplot2") 
library("ggpubr")  # For statistical annotations on plots
library("gridExtra") # For combining multiple plots

# Import BIOM data  
file_biom <- ("17_bean_metasamples.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Displays importes data

# Read metadata file
metdata <- read.table("MetaData_Mod.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata)

# Adjust the order of factor levels for "Soil_stage_Seeding"
levels(metdata$Soil_Stage_Seeding) # View factor levels
metdata$Soil_Stage_Seeding <- factor(metdata$Soil_Stage_Seeding, levels = c('A_Soil','A_Rhizo', 'N_Soil', 'N_Rhizo')) # Con esta línea se modifica el orden
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
write.csv(richness, file="D:/Documentos/Doctorado/Doctorado_Octavo_Semestre/Manuscrito/Manuscrito_Pinto_Saltillo/Graficas_R/R_Result/Alfa_diversity_Dominio_Bacteria.csv")

# Define pairwise comparisons for statistical testing  
comparition <- list(c('A_Soil','A_Rhizo'))
compare <- list(c('N_Soil', 'N_Rhizo'))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# Plot Chao1 alpha diversity with statistical test 
Observed <- plot_richness(merged_metagenomes_Bacteria, x= "Soil_Stage_Seeding", measures= "Chao1") + #,"Shannon", "Simpson","InvSimpson")) +
            geom_boxplot(alpha=0.4, aes(color= Soil_Stage_Seeding))+
            geom_point(alpha=0.7, size =2.5, aes( color= Soil_Stage_Seeding)) +
            theme(legend.position="none",
                text=element_text(colour="black", size=12),
                axis.text.x=element_text(colour="black", size = 10), 
                axis.text.y=element_text(colour="black", size = 10),
                #axis.ticks.x = element_blank(),
                #axis.text.x =  element_blank(), 
                plot.margin = margin(1, 1, 1, 1, "cm")) +
            theme_bw() +
            stat_compare_means(method = "wilcox.test", comparisons = comparition, label = "p.signif", 
                symnum.args = symnum.args) +
            stat_compare_means(method = "wilcox.test", comparisons = compare, label = "p.signif", 
                symnum.args = symnum.args)

Graph_Observed <- Observed + 
     labs(x = "" , y = "") +
     #scale_fill_manual(values = c("#997700", "#004488"), guide="none") +
     #scale_fill_manual(values = c("#997700", "#004488"))
     scale_color_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"), guide= "none") +
     scale_fill_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"))

# Plot Shannon Index alpha diversity
Shannon <- plot_richness(merged_metagenomes_Bacteria, x= "Soil_Stage_Seeding", measures= "Shannon") +
            geom_boxplot(alpha=0.4, aes(color= Soil_Stage_Seeding))+
            geom_point(alpha=0.6, aes( color= Soil_Stage_Seeding)) +
            theme(legend.position="none",
                text=element_text(colour="black", size=12),
                axis.text.x=element_text(colour="black", size = 10), 
                axis.text.y=element_text(colour="black", size = 10),
                #axis.text.x = "none" ,#element_text (size= 20),
                #axis.title.y =  "none",  #element_text(size = 12), 
                plot.margin = margin(1, 1, 1, 1, "cm")) +
            theme_bw()+
            stat_compare_means(method = "wilcox.test", comparisons = comparition, label = "p.signif", 
                symnum.args = symnum.args)+
            stat_compare_means(method = "wilcox.test", comparisons = compare, label = "p.signif", 
                symnum.args = symnum.args)

# Plot InvSimpson alpha diversity
Graph_Shannon <- Shannon + 
     labs(x = "" , y = "") +
     #scale_color_manual(values = c("#997700", "#004488"), guide= "none")+
     #scale_fill_manual(values = c("#997700", "#004488"))
     scale_color_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"), guide= "none") +
     scale_fill_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"))


InvSimpson <- plot_richness(merged_metagenomes_Bacteria, x= "Soil_Stage_Seeding", measures= "InvSimpson") +
            geom_boxplot(alpha=0.4, aes(color= Soil_Stage_Seeding))+
            geom_point(alpha=0.6, aes( color= Soil_Stage_Seeding)) +
            theme(legend.position="none",
                text=element_text(colour="black", size=12),
                axis.text.x=element_text(colour="black", size = 10), 
                axis.text.y=element_text(colour="black", size = 10),
                #axis.text.x = element_text (size= 20),
                #axis.title.y = element_text(size = 12), 
                plot.margin = margin(1, 1, 1, 1, "cm")) +
            theme_bw()+
            stat_compare_means(method = "wilcox.test", comparisons = comparition, label = "p.signif", 
                symnum.args = symnum.args) +
            stat_compare_means(method = "wilcox.test", comparisons = compare, label = "p.signif", 
                symnum.args = symnum.args)

Graph_InvSimpson <- InvSimpson + 
     labs(x = "" , y = "") +
     #scale_color_manual(values = c("#997700", "#004488"), guide= "none")+
     #scale_fill_manual(values = c("#997700", "#004488"))
     scale_color_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"), guide= "none") +
     scale_fill_manual(values = c("#DABA22", "#533D43", "#258B08", "#AF1B97"))

# Generate combined plots into a single figure
combined_plot <- grid.arrange(Graph_Observed, Graph_Shannon, Graph_InvSimpson, nrow =1)

# Save combined plot to a PDF file
pdf(file="D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Pinto_Saltillo_Graph/R_Result/Diversidad_Alfa.pdf",
family= "Times", width =  8, height = 5 , paper= "USr")
print(combined_plot)
 dev.off()
