library("ggvenn")
library("ggplot2")
library("gridExtra")

TableAll <- read.table(file="VennData_Tab_Pinto_Negro_Bayo_Genus.txt", header=TRUE, sep="\t", dec=".")

# Check for empty values
sum(is.na(TableAll))
sum(TableAll == "")
#sum(TableAll$DS_PS_Genus == "")  # Cantidad de valores vacíos en col1
#sum(TableAll$AN_PS_Genus == "")  # Cantidad de valores vacíos en col2

# Replace NA and empty values ??with zero
TableAll[is.na(TableAll) | TableAll == ""] <- 0
View(TableAll)
# Check again for empty values
sum(is.na(TableAll))  # Should return 0
sum(TableAll == "")   # Should return 0

vennPS <- list(
  DESeq2_PS = TableAll$DS_PS_Genus[TableAll$DS_PS_Genus != 0],
  ANCOMBC_PS = TableAll$AN_PS_Genus[TableAll$AN_PS_Genus != 0]
)    # ggvenn(vennPS, fill_color = c("red", "blue"), digits = 2)
vennPS_plot <- ggvenn(vennPS, fill_color=c("white", "white"), set_name_color="black", set_name_size=5, text_color="black", text_size=5) + theme(plot.margin = unit(c(1, 1, 1, 1),"cm"))
vennPS_plot

# Calculate intersections and unique values
only_A <- setdiff(TableAll$DS_PS_Genus, TableAll$AN_PS_Genus)
only_B <- setdiff(TableAll$AN_PS_Genus, TableAll$DS_PS_Genus)
A_B <- intersect(TableAll$DS_PS_Genus, TableAll$AN_PS_Genus)

# Show results
cat("Only DS_PS_Genus:", length(only_A), "\n")
cat("Only AN_PS_Genus:", length(only_B), "\n")
cat("Intersection:", length(A_B), "\n")

# Save each set in separate files
#writeLines(paste(only_A, collapse = "\n"), "VennData_Pinto_onlyDS.txt")
#writeLines(paste(only_B, collapse = "\n"), "VennData_Pinto_onlyAN.txt")
#writeLines(paste(A_B, collapse = "\n"), "VennData_Pinto_DS_and_AN.txt")

vennNe <- list(
  DESeq2_Negro = TableAll$DS_Ne_Genus[TableAll$DS_Ne_Genus != 0],
  ANCOMBC_Negro = TableAll$AN_Ne_Genus[TableAll$AN_Ne_Genus != 0]
)
vennNe_plot <- ggvenn(vennNe, fill_color=c("white", "white"), set_name_color="black", set_name_size=5, text_color="black", text_size=5) + theme(plot.margin = unit(c(1, 1, 1, 1),"cm"))
vennNe_plot

# Calculate intersections and unique values
only_An <- setdiff(TableAll$DS_Ne_Genus, TableAll$AN_Ne_Genus)
only_Bn <- setdiff(TableAll$AN_Ne_Genus, TableAll$DS_Ne_Genus)
An_Bn <- intersect(TableAll$DS_Ne_Genus, TableAll$AN_Ne_Genus)

# Show results
cat("Only DS_Ne_Genus:", length(only_An), "\n")
cat("Only AN_Ne_Genus:", length(only_Bn), "\n")
cat("Intersection:", length(An_Bn), "\n")

# Save each set in separate files
#writeLines(paste(only_An, collapse = "\n"), "VennData_Negro_onlyDS.txt")
#writeLines(paste(only_Bn, collapse = "\n"), "VennData_Negro_onlyAN.txt")
#writeLines(paste(An_Bn, collapse = "\n"), "VennData_Negro_DS_and_AN.txt")

vennBa <- list(
  DESeq2_Bayo = TableAll$DS_Ba_Genus[TableAll$DS_Ba_Genus != 0],
  ANCOMBC_Bayo = TableAll$AN_Ba_Genus[TableAll$AN_Ba_Genus != 0]
)
vennBa_plot <- ggvenn(vennBa, fill_color=c("white", "white"), set_name_color="black", set_name_size=5, text_color="black", text_size=5) + theme(plot.margin = unit(c(1, 1, 1, 1),"cm"))
vennBa_plot

# Calculate intersections and unique values
only_Ab <- setdiff(TableAll$DS_Ba_Genus, TableAll$AN_Ba_Genus)
only_Bb <- setdiff(TableAll$AN_Ba_Genus, TableAll$DS_Ba_Genus)
Ab_Bb <- intersect(TableAll$DS_Ba_Genus, TableAll$AN_Ba_Genus)

# Show results
cat("Only DS_Ba_Genus:", length(only_Ab), "\n")
cat("Only AN_Ba_Genus:", length(only_Bb), "\n")
cat("Intersection:", length(Ab_Bb), "\n")

# Save each set in separate files
#writeLines(paste(only_Ab, collapse = "\n"), "VennData_Bayo_onlyDS.txt")
#writeLines(paste(only_Bb, collapse = "\n"), "VennData_Bayo_onlyAN.txt")
#writeLines(paste(Ab_Bb, collapse = "\n"), "VennData_Bayo_DS_and_AN.txt")

# Create plot
grid.arrange(vennPS_plot, vennNe_plot, vennBa_plot, nrow=3)
pdf("VennData_Diagrams.pdf", width = 8, height = 11, family="Times")
grid.arrange(vennPS_plot, vennNe_plot, vennBa_plot, nrow=3)
dev.off()


