# We analyzed shifts in bacterial genera from the bulk soil to the rhizosphere using DESeq2 to determine 
# if the rhizosphere effect on bacterial community
# composition in the common bean (Phaseolus vulgaris) cultivar 'Pinto Saltillo'

# Load required libaries
library("phyloseq") # For microbiome data analysis
library("ggplot2") # For plotting
library("ggrepel") # For improved label visualization
library("EnhancedVolcano") # For creating volcano plots
library("DESeq2") # For differential abundance analysis

# Import BIOM data
file_biom <- ("Agricola_NonAgricola_12Sam.biom")
data <- import_biom(file_biom, parseFunction=parse_taxonomy_default)
print(data) # Display imported data

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
 
# Select data at the Genus level
Genus <-  tax_glom(merged_metagenomes_Bacteria, "Genus")
print(Genus) # Display aggregated data

# Convert phyloseq data to DESeq2 format for differential abundance analysis
Deq_B = phyloseq_to_deseq2(Genus, ~ Sowing) # Set design formula for DESeq2
print(Deq_B) # Display summary of the DESeq2 format

# Perform differential abundance analysis
Deq_B = DESeq(Deq_B)

# Extract results with specific contrast and significance cutoff
resB = results(Deq_B, contrast=c("Sowing", "After-sowing", "Before-sowing")  , alpha= 0.05) 
             
# Adjust significance level
resB = resB[order(resB$padj),] # Sort results by adjusted p-value
dim(resB) # Preview top results

# Map Genus names for significant taxa
specie_table <- read.table("IDGenus_Pinto_Saltillo", header = TRUE, sep = "\t")
names(specie_table)[1] <- "ID" 
specie_map <- specie_table[match(rownames(resB), specie_table$ID),]


# Generate a volcano plot
PS_Volcano_Genus <- EnhancedVolcano(resB,
    lab = specie_map$Genus,
    x = 'log2FoldChange',
    y = 'pvalue',
   selectLab = c( 'Flavobacterium','Rhizobium', 'Stenotrophomonas',
                'Erwinia', 'Neorhizobium','Xanthomonas', 'Pseudomonas', 'Lysobacter', 'Ensifer', 'Variovorax','Pantoea'
                     ),
    pCutoff = 0.05,
    drawConnectors = TRUE,
    legendPosition = 'top', 
    FCcutoff=log2(2),
    legendLabSize = 10,
    title = "Pinto Saltilo",
    subtitle = ""
    )
    
# Save the volcano plot to a PDF
pdf(file="Volcano_Plot_PintoSaltillo_GenNames.pdf",
family= "Times", width = 15, height = 25, paper = "USr")
print(PS_Volcano_Genus)
dev.off()


## Further modifications to the PDF file were done using Adobe Illustrator


