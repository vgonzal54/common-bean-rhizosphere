# This script processes metagenomic data to analyze the bacterial communities in agricultural and non-agricultural soils.
# It utilizes differential abundance analysis and generates a volcano plot to compare bacterial genera based on their relative abundance.

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("ggrepel") # For improved label visualization
library("EnhancedVolcano") # For creating volcano plots 
library("DESeq2") # For differential abundance analysis

# Import BIOM data
file_biom <- ("Agricola_NonAgricola_12Sam.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Display imported dat

# Read metadata file
metdata <- read.table("MetaData_Agricola_NonAgricola_12Sam.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)
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

# Convert phyloseq data to DESeq2 format for differential abundance analysis
Deq_B = phyloseq_to_deseq2(merged_metagenomes_Bacteria, ~ Sowing) # Set design formula for DESeq2
print(Deq_B) # Display summary of the DESeq2 format

# Perform differential abundance analysis
Deq_B = DESeq(Deq_B)

# Extract results with specific contrast and significance cutoff
resB = results(Deq_B, contrast=c("Sowing", "After-sowing", "Before-sowing")
	, alpha= 0.05) 

# Adjust significance level
resB = resB[order(resB$padj),] # Sort results by adjusted p-value
dim(resB)

# Combine results with taxonomy information
tabla = cbind(as(resB, "data.frame"), as(tax_table(merged_metagenomes_Bacteria)[rownames(resB), ], "matrix")) # Add taxonomy to results

# Save results and export results to TSV format
write.table(tabla, file ="D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Volcano_Plot/Volcano_Plot_modificado/R_Results/Pinto_saltillo_LFC.tsv" , row.names=TRUE, sep="\t")

# Preview top results
head(resB)

# Generate a volcano plot
Volcano_plot <- EnhancedVolcano(resB,
     			lab = NA,  # Exclude labels for data points
     			legendPosition = 'top',
    			x = 'log2FoldChange',
     			y = 'pvalue',
     			pCutoff=0.05, # Significance cutoff
     			FCcutoff=log2(2), # Fold-change cutoff
     			legendLabSize = 10,
    			title = "Bulk soil Vs Pinto Saltilo rhizosphere",
     			subtitle = ""
    			)

# Save the volcano plot to a PDF
pdf(file= "D:/Documentos/Doctorado/Doctorado_Decimo_semestre/Gráficas_Variedades/Volcano_Plot/Volcano_Plot_modificado/R_Results/Volcano_Plot_Pinto_Saltillo.pdf", 
 	family= "Times", width = 20 , height = 15 , paper= "letter")
print(Volcano_plot) # Print the volcano plot
dev.off()