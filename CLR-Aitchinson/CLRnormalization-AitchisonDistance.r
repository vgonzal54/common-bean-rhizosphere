library(vegan)
library(ape)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(compositions)
library(pairwiseAdonis)
library("QsRutils")
library(tibble)

archivo_biom <- ("All_PintoSalt_Accecion.biom")
alldata <- import_biom(archivo_biom, parseFunction=parse_taxonomy_default)
metdata <- read.table("Matadata_all_PintoSAlt_Acc.tsv", row.names=1, header=TRUE, sep='\t', check.names=TRUE, stringsAsFactors=TRUE)

levels(metdata$Bean_name) 
metdata$Bean_name <- factor(metdata$Bean_name, levels = c('MetaG-PAA','MetaG-PAD','MetaG-PAI','MetaG-PNA','MetaG-PND','MetaG-PNI',
'MetaG-AA','MetaG-AB','MetaG-AD','MetaG-AE','MetaG-AH','MetaG-AI','MetaG-NA','MetaG-ND',
'MetaG-NE','MetaG-NH','MetaG-NI', 'VP1-2', 'VP3-4', 'VP5-6', 'Negro_V10', 'Negro_V26', 'Negro_V2', 'Negro_V45', 'Negro_V4',
'Bayo_V16', 'Bayo_V22', 'Bayo_V43', 'Bayo_V47')) # This line changes the order
levels(metdata$Bean_name)
levels(metdata$Type) 
metdata$Type <- factor(metdata$Type, levels = c('Bayo_bean','Black_bean','Pinto_Saltillo','Bayo_and_Black_soil','Pinto_Saltillo_soil'))
levels(metdata$Type) 

merged_metagenomes <- phyloseq(otu_table(alldata), tax_table(alldata), sample_data(metdata)) 
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
merged_metagenomes
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10580 taxa and 29 samples ]
## sample_data() Sample Data:       [ 29 samples by 3 sample variables ]
## tax_table()   Taxonomy Table:    [ 10580 taxa by 7 taxonomic ranks ]

merged_metagenomes_Bacteria <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria") 
merged_metagenomes_Bacteria
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 7083 taxa and 29 samples ]
## sample_data() Sample Data:       [ 29 samples by 3 sample variables ]
## tax_table()   Taxonomy Table:    [ 7083 taxa by 7 taxonomic ranks ]


## Distance analysis with Aitchison's method

gloglom_Genus <- tax_glom(merged_metagenomes_Bacteria, taxrank = 'Genus')
gloglom_Genus
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1471 taxa and 29 samples ]
#sample_data() Sample Data:       [ 29 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1471 taxa by 7 taxonomic ranks ]

#View(otu_table(gloglom_Genus))
#View(tax_table(gloglom_Genus))

# To save the table of OTUs and taxa of the genus level
reads_num <- cbind(data.frame(otu_table(gloglom_Genus)), data.frame(tax_table(gloglom_Genus)))
View(reads_num)
reads_num <- tibble::rownames_to_column(reads_num, var = "TaxonID")
View(reads_num)
write.table(reads_num, file = "29_Abundancia_genus.txt", row.names=FALSE, sep="\t")

# Flip the data table (because if I do it as is, an error appears when running subsequent analyses)
mat <- veganotu(gloglom_Genus)
otu_table(gloglom_Genus) <- otu_table(mat, taxa_are_rows=F)
View(otu_table(gloglom_Genus))

# Convert otus to a data frame before doing the transformation to CLR
clr_data <- as.data.frame(otu_table(gloglom_Genus))
View(clr_data)
# If there are zeros in the data, add a small value to the abundances. Replace zeros with a small value
clr_data[clr_data == 0] <- 1e-6
View(clr_data)
# Transform data to CLR
clr_data <- microbiome::transform(clr_data, "clr")
View(clr_data)

# Convert to an object of type `acomp` from the compositions package
clr_comp <- acomp(clr_data)
# Calculate the Aitchison distance
distance_matrix <- dist(clr_comp)
View(as.matrix(distance_matrix))


## PCoA analysis on the Aitchison distance matrix

pcoa_result <- pcoa(distance_matrix)

# Extracting the coordinates of the samples
pcoa_scores <- as.data.frame(pcoa_result$vectors)
View(pcoa_scores)
pcoa_scores$Sample <- rownames(pcoa_scores)
View(pcoa_scores)
pcoa_scores <- cbind(pcoa_scores, metdata)
View(pcoa_scores)

# Obtain the percentage of variance explained by the axes
variance_explained <- pcoa_result$values$Relative_eig * 100
var_exp_axis1 <- round(variance_explained[1], 2)
var_exp_axis2 <- round(variance_explained[2], 2)

# Plot the two principal coordinates
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, color=Type, shape=Type)) + geom_point(size=3) + 
scale_color_manual(values=c('#F8766D','#00BA38','#619CFF','#DB72FB','#D39200','#FF61C3')) + 
scale_shape_manual(values=c(16,17,15,18,4,25)) + theme_bw() + 
labs(title="PCoA of Aitchison Distance Matrix\n", x=paste0("Axis1 (", var_exp_axis1, "%)"), y=paste0("Axis2 (", var_exp_axis2, "%)")) + 
theme(text=element_text(colour="black", size=12), axis.text=element_text(colour="black", size=10)) + 
theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
pdf("29_CLR_transformado_PCoA_Aitchison.pdf", width=8, height=7, family="Times")
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, color=Type, shape=Type)) + geom_point(size=3) + 
scale_color_manual(values=c('#F8766D','#00BA38','#619CFF','#DB72FB','#D39200','#FF61C3')) + 
scale_shape_manual(values=c(16,17,15,18,4,25)) + theme_bw() + 
labs(title="PCoA of Aitchison Distance Matrix\n", x=paste0("Axis1 (", var_exp_axis1, "%)"), y=paste0("Axis2 (", var_exp_axis2, "%)")) + 
theme(text=element_text(colour="black", size=12), axis.text=element_text(colour="black", size=10)) + 
theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
dev.off()


# PERMANOVA statistical analysis

adonis2(distance_matrix ~ metdata$Type, perm=999, p.adjust.m='holm')
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = distance_matrix ~ metdata$Type, permutations = 999, p.adjust.m = "holm")
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     4   6379.7 0.35991 3.3736  0.001 ***
#Residual 24  11346.4 0.64009                  
#Total    28  17726.1 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.adonis(distance_matrix,factors=metdata$Type, p.adjust.m='holm')
#Set of permutations < 'minperm'. Generating entire set.
#                                         pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1        Pinto_Saltillo vs Soil_Pinto_Saltillo  1 2228.6645 7.0658890 0.32021773   0.001      0.010   *
#2                      Pinto_Saltillo vs Negro  1 1116.7454 2.4580041 0.14935007   0.003      0.021   .
#3                       Pinto_Saltillo vs Bayo  1 1183.8507 2.5302810 0.16292564   0.005      0.025   .
#4       Pinto_Saltillo vs Soil_other_varieties  1 1463.3413 4.3856283 0.26765091   0.002      0.018   .
#5                 Soil_Pinto_Saltillo vs Negro  1 2571.4371 5.1006253 0.36173043   0.002      0.018   .
#6                  Soil_Pinto_Saltillo vs Bayo  1 2562.2792 4.8129183 0.37563014   0.003      0.021   .
#7  Soil_Pinto_Saltillo vs Soil_other_varieties  1  456.6415 1.4658358 0.17314720   0.037      0.111    
#8                                Negro vs Bayo  1  631.1983 0.7503484 0.09681479   1.000      1.000    
#9                Negro vs Soil_other_varieties  1 1792.2685 2.8223910 0.31991225   0.046      0.111    
#10                Bayo vs Soil_other_varieties  1 1853.5296 2.6240314 0.34417899   0.024      0.096    

# Betadisper analysis
dispersion <- betadisper(d=distance_matrix, group=metdata$Type, type="centroid")
dispersion
#        Homogeneity of multivariate dispersions
#Call: betadisper(d = distance_matrix, group = metdata$Type, type = "centroid")
#No. of Positive Eigenvalues: 28
#No. of Negative Eigenvalues: 0
#Average distance to centroid:
#                Bayo                Negro       Pinto_Saltillo Soil_other_varieties  Soil_Pinto_Saltillo 
#               26.47                24.70                16.58                15.54                15.54 
#Eigenvalues for PCoA axes:
#(Showing 8 of 28 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#4618.7 1496.3 1133.4 1094.8  926.9  868.4  770.9  717.8

permutest(dispersion, permutations=999, pairwise=TRUE)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     4 558.12  139.53 11.225    999  0.001 ***
#Residuals 24 298.33   12.43                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Pairwise comparisons:
#(Observed p-value below diagonal, permuted p-value above diagonal)
#                           Bayo      Negro Pinto_Saltillo Soil_other_varieties Soil_Pinto_Saltillo
#Bayo                            3.0300e-01     6.0000e-03           2.0000e-03               0.001
#Negro                2.8009e-01                9.0000e-03           8.0000e-03               0.001
#Pinto_Saltillo       2.1893e-03 5.0332e-03                          7.1200e-01               0.629
#Soil_other_varieties 2.7063e-05 1.9947e-03     7.3791e-01                                    0.986
#Soil_Pinto_Saltillo  1.8477e-07 4.4383e-05     6.3052e-01           9.9333e-01                    


