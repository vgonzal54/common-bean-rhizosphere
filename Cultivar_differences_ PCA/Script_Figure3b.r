# This script performs a Principal Component Analysis (PCA) to explore morphological similarities among
# common bean (Phaseolus vulgaris) accessions. 

# Load required libraries

library("dplyr") # For data manipuation
library("FactoMineR") # For performing Principal Component Analysis (PCA)
library("factoextra") # For vizualizing PCA results
library("ggsci") # For color palettes in visualizations
library("vegan") # For ecological analysis (permanova analysis)

# Load the accessions data table 
Aseccion_Table <- read.table(file= "asecciones_replicas.txt", sep="\t", dec=".", header=TRUE, stringsAsFactors = TRUE) 

# Assign the first column as row names
row.names(Aseccion_Table) <- Aseccion_Table[,1]

#  Remove the first column to retain only the variables
PS_Table <- Aseccion_Table[,-1]

# Filter the variables of interest
Ps_Table_Filter <- select(PS_Table, PDH,CDH,CCE,IPH,HDC,CDE,CDA,CDV,CDS,CSS,DCS,NCS,PTS,LDS,ADS,EDS,FTS,FLS,PDS)

# Create a vector with the common names of the bean types
Cultivar <- c("Black_bean","Black_bean","Black_bean","Black_bean","Black_bean","Black_bean","Black_bean","Black_bean","Black_bean",
						"Bayo_bean","Bayo_bean","Bayo_bean","Bayo_bean","Bayo_bean","Bayo_bean","Black_bean","Black_bean","Black_bean",
						"Bayo_bean","Bayo_bean","Bayo_bean","Black_bean","Black_bean","Black_bean",
						"Bayo_bean","Bayo_bean","Bayo_bean","Pinto_Saltillo")

# Add the common names vector to the filtered table
New_table_rep <- cbind(Ps_Table_Filter,Cultivar)

# Perform Principal Component Analysis (PCA)
table_rep.pca <- PCA(New_table_rep[,-20], graph = FALSE)

# Create the PCA plot
PCA_plot <- fviz_pca_ind(table_rep.pca,
             geom.ind = "point" , # show points only (but not "text")
             fill.ind= New_table_rep$Cultivar, # Color point borders basen on beans type
             col.ind= New_table_rep$Cultivar, # color by groups
             legend.title = "Cultivar", # Legend title
             mean.point = FALSE, # No display the centroid
             pointsize = 3, # Size of the points 
             repel = FALSE, # Deactivate label overlap prevention
             title = "") + # plot title
             theme_bw() + scale_shape_manual(values = c(16, 17, 15)) # Shape of the points 

# Save the plot to a PDF fil
pdf(file= "PCA_morphological_beans.pdf", 
    family= "Times", width = 8 , height = 8 , paper= "letter")
print(PCA_plot)
dev.off()


## Further modifications to the PDF file were done using Adobe Illustrator


# Statistical analysis using PERMANOVA (adonis2)

variables <- New_table_rep[,-20]
group <- as.factor(New_table_rep$Cultivar)

adonis2(variables ~ group, data = New_table_rep, method = "euclidean", permutations = 999, p.adjust.m='holm')
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = variables ~ group, data = New_table_rep, permutations = 999, method = "euclidean", p.adjust.m = "holm")
#         Df SumOfSqs      R2     F Pr(>F)    
#Model     2   2627.6 0.60026 18.77  0.001 ***
#Residual 25   1749.9 0.39974                 
#Total    27   4377.5 1.00000                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

