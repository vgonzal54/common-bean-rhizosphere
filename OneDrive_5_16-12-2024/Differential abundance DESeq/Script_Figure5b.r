# This script is to investigate whether the rhizosphere effect influences bacterial community 
# composition in the common bean (Phaseolus vulgaris) cultivar 'Bayo', 
# using differential abundance analysis.
             
# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("ggrepel") # For improved label visualization
library("EnhancedVolcano")  # For creating volcano plots
library("DESeq2") # For differential abundance analysis

# Import BIOM data
file_biom <- ("Bayo_Rhizo_Soil.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Display imported data

# Read metadata file
metdata <- read.table("MetaData_Bayo_Rhizo_Soil.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
head(metdata) # Preview metadata structure

# Combine OTU table, taxonomic table, and metadata into a single phyloseq object
merged_metagenomes <- phyloseq(otu_table(data), tax_table(data), sample_data(metdata)) 
print(merged_metagenomes) # Display summary of the combined object

# Rename columns in the taxonomic table
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter data to include only the "Bacteria" Kingdom
merged_metagenomes_Bacteria <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria") 
print(merged_metagenomes_Bacteria) 
 
# Select data at the Genus level
Genus <-  tax_glom(merged_metagenomes_Bacteria, "Genus")
print(Genus) # Display aggregated data

# Convert phyloseq data to DESeq2 format for differential abundance analysis
Deq_B_Genus = phyloseq_to_deseq2(Genus, ~ Origin_sample) # Set design formula for DESeq2
print(Deq_B_Genus) # Display summary of the DESeq2 format

# Perform differential abundance analysis
Deq_B_Genus = DESeq(Deq_B_Genus)

# Extract results with specific contrast and significance cutoff
resB_Genus = results(Deq_B_Genus, contrast=c("Origin_sample", "Rhizosphere", "Soil")  , alpha= 0.05) 

# Adjust significance level              
resB_Genus = resB_Genus[order(resB_Genus$padj),] # Sort results by adjusted p-value
dim(resB_Genus) # Preview top results

# Map Genus names for significant taxa
specie_table <- read.table("IDGenus_Bayo", header = TRUE, sep = "\t")
names(specie_table)[1] <- "ID" 
specie_map <- specie_table[match(rownames(resB_Genus), specie_table$ID),]


# Generate a volcano plot
Bayo_volcano_Genus <- EnhancedVolcano(resB_Genus,
    lab = specie_map$Genus,
    x = 'log2FoldChange',
    y = 'pvalue',
    #xlim = c(-2, 10),
    selectLab = c( 'Flavobacterium','Rhizobium', 'Stenotrophomonas',
                'Erwinia', 'Neorhizobium','Xanthomonas', 'Pseudomonas', 'Lysobacter', 'Ensifer', 'Variovorax', 'Pantoea'
                     ),
    pCutoff = 0.05,
    drawConnectors = TRUE,
    legendPosition = 'top', 
    FCcutoff=log2(2),
    legendLabSize = 10,
    title = "Bayo",
    subtitle = ""
    )
    

# Save the volcano plot to a PDF
pdf(file="D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Volcano_Plot/Volcano_Plot_modificado/R_Results/Volcano_Plot_Bayo_GenNames.pdf",
family= "Times", width = 15 , height = 25 , paper = "USr")
print(Bayo_volcano_Genus)
dev.off()
