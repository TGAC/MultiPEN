# It requieres R Packages:
# clusterProfiler, https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# BBmisc
# KEGG.db
# org.Hs.eg.db

# Check if packages are installed, otherwise, install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("BBmisc", "clusterProfiler")

# TEST:
# file with example data: 'ExampleOutputs/MultiPEN-Rankings_lambda0.1.txt'
# Run in terminal:
# Rscript enrichmentGO.R '/path-to-file/MultiPEN-Rankings_lambda0.0001.txt' output-folder/

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "output_MultiPEN/enrichment-GO/"
}

dataFile <- args[1]
cat(sprintf("Loading data from file: %s\n", dataFile))

outputDir <- args[2]
if(!dir.exists(outputDir)){
  cat(sprintf("Creating output directory: %s\n", outputDir))
  dir.create(outputDir)
}

# Load file with data for analysis
data <- read.table(dataFile, header=TRUE, sep = "\t", stringsAsFactors=FALSE)
cat(sprintf("Number of features: %i\n",nrow(data)))


library(clusterProfiler)
library(BBmisc)
library(GO.db)

D <- sortByCol(data, 'ranking')
D <- D[D[,2]!=0,]
D <- D[,c(1,2,3)]  # Only interested in the first three columns [name, weight, ranking]
entrez<-bitr(D$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE)
ranked<-merge(D,entrez,by.x='name',by.y='SYMBOL')
ranked <- sortByCol(ranked, 'weight', asc = F)
geneList <- ranked$weight
names(geneList) <- ranked$ENTREZID

cat(sprintf("Performing over-representation analysis (with KEGG) ...  "))
cat(sprintf("Results saved to folder: %s\n", outputDir))

#Enrichment with KEGG
enrichment_kegg <- enrichKEGG(ranked$ENTREZID, organism = 'hsa', keyType = "kegg")
results <- summary(enrichment_kegg)
head(results)


# modify column names for consistency for valid MATLAB identifiers
# change: p.adjust to pAdjust
aux <- colnames(results)
aux[c(5,6,7)] <- c("pValue", "pAdjust", "qValue")
colnames(results) <- aux

#write results (table) to file: 
fileName <- paste(outputDir, "enrichment-KEGG.txt", sep = "")
cat(sprintf("writing results to file: %s\n", fileName))
write.table(results, fileName, sep = '\t', row.names = FALSE)

# save plot to file:
fileName <- paste(outputDir, 'enrichment-KEGG.pdf', sep = "")
pdf(fileName)
barplot(enrichment_kegg, showCategory=20)


## Gene Set Enrichment with KEGG
kk2 <- gseKEGG(geneList, organism = 'hsa', keyType = "kegg")
results <- summary(kk2)
results


# modify column names for consistency for valid MATLAB identifiers
# change: p.adjust to pAdjust
aux <- colnames(results)
aux[c(6,7,8)] <- c("pValue", "pAdjust", "qValue")
colnames(results) <- aux

#write results (table) to file: 
fileName <- paste(outputDir, "gse-KEGG.txt", sep = "")
cat(sprintf("writing results to file: %s\n", fileName))
write.table(results, fileName, sep = '\t', row.names = FALSE)

# save plot to file:
fileName <- paste(outputDir, 'gse-KEGG.pdf', sep = "")
pdf(fileName)
library("DOSE")
barplot(results$setSize, horiz = T, legend.text = T)