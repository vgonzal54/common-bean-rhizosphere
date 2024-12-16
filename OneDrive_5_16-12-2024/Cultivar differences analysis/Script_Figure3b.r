# This script performs a Principal Component Analysis (PCA) to explore morphological similarities among common bean (Phaseolus vulgaris)
# accessions. 

# Load required libraries

library("dplyr") # For data manipuation
library("FactoMineR") # For performing Principal Component Analysis (PCA)
library("factoextra") # For vizualizing PCA results
library("ggsci") # For color palettes in visualizations

# Load the accessions data table 
Aseccion_Table <- read.table(file= "asecciones_replicas.txt", sep="\t", dec=".", header=TRUE, stringsAsFactors = TRUE) 

# Assign the first column as row names
row.names(Aseccion_Table) <- Aseccion_Table[,1]

#  Remove the first column to retain only the variables
PS_Table <- Aseccion_Table[,-1]


# Filter the variables of interest
Ps_Table_Filter <- select(PS_Table, PDH,CDH,CCE,IPH,HDC,CDE,CDA,CDV,CDS,CSS,DCS,NCS,PTS,LDS,ADS,EDS,FTS,FLS,PDS)

# Create a vector with the common names of the bean types
Bean_Type <- c("Negro","Negro","Negro","Negro","Negro","Negro","Negro","Negro","Negro",
						"Bayo","Bayo","Bayo","Bayo","Bayo","Bayo","Negro","Negro","Negro",
						"Bayo","Bayo","Bayo","Negro","Negro","Negro",
						"Bayo","Bayo","Bayo","Pinto_Saltillo")

# Add the common names vector to the filtered table
New_table_rep <- cbind(Ps_Table_Filter,Bean_Type)

# Perform Principal Component Analysis (PCA)
table_rep.pca <- PCA(New_table_rep[,-20], graph = FALSE)

# Create the PCA plot
PCA_plot <- fviz_pca_ind(table_rep.pca,
             geom.ind = "point" , # show points only (but not "text")
             palette = c("red", "green", "Blue"), # Color palette by groups
             fill.ind= New_table_rep$ Bean_Type, # Color point borders basen on beans type
             col.ind= New_table_rep$ Bean_Type, # color by groups
             legend.title = "Bean common names", # Legend title
             mean.point = FALSE, # No display the centroid
             pointshape = 10, # Shape of the points
             pointsize = 3, # Size of the points 
             repel = FALSE, # Deactivate label overlap prevention
             title = "Bean aseccions PCA") # plot title

# Save the plot to a PDF fil
pdf(file= "D:/Documentos/Doctorado/Doctorado_Noveno_Semestre/Asecciones/R_Result/PCA_nombres_comunes.pdf", 
    family= "Times", width = 8 , height = 8 , paper= "letter")
print(PCA_plot)
dev.off()