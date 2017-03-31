# It requieres R Packages:
# clusterProfiler, https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# BBmisc
# GO.db
# org.Hs.eg.db

# Check if packages are installed, otherwise, install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("BBmisc", "GO.db", "org.Hs.eg.db", "clusterProfiler")

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
D <- D[,c(1,3)]  # Only interested in the first three columns [name, weight, ranking]
entrez<-bitr(D$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE)
ranked<-merge(D,entrez,by.x='name',by.y='SYMBOL')

cat(sprintf("Performing over-representation analysis (enrichGO) ...  "))
cat(sprintf("Results saved to folder: %s\n", outputDir))

#Enrichment for subontology BP (Biological Process)
subclassOnt <- "BP"
enrichment_BP <- enrichGO(ranked$ENTREZID, OrgDb="org.Hs.eg.db", ont=subclassOnt, readable=TRUE)
enrichmentSummary_BP <- summary(enrichment_BP)
head(enrichmentSummary_BP)
if(nrow(enrichmentSummary_BP)>0){
  aux <- cbind(enrichmentSummary_BP, "BP")
  colnames(aux)[10]<- 'subontology'
  enrichmentSummary_BP <- aux        
}  
#add to results: enrichment for BP category
results <- enrichmentSummary_BP

#Enrichment for subontology MF (Molecular Function)
subclassOnt <- "MF"
enrichment_MF <- enrichGO(ranked$ENTREZID, OrgDb="org.Hs.eg.db", ont=subclassOnt, readable=TRUE)
enrichmentSummary_MF <- summary(enrichment_MF)
head(enrichmentSummary_MF)
if(nrow(enrichmentSummary_MF)>0){
  aux <- cbind(enrichmentSummary_MF, "MF")
  colnames(aux)[10]<- 'subontology'
  enrichmentSummary_MF <- aux       
  #add to results: enrichment for MF category
  results <- rbind(results, enrichmentSummary_MF)
}  


#Enrichment for subclass CC (Cellular Component)
subclassOnt <- "CC"
enrichment_CC <- enrichGO(ranked$ENTREZID, OrgDb="org.Hs.eg.db", ont=subclassOnt, readable=TRUE)
enrichmentSummary_CC <- summary(enrichment_CC)
head(enrichmentSummary_CC)
if(nrow(enrichmentSummary_CC)>0){
  aux <- cbind(enrichmentSummary_CC, "CC")
  colnames(aux)[10]<- 'subontology'
  enrichmentSummary_CC <- aux
  #add to results: enrichment for CC category
  results <- rbind(results, enrichmentSummary_CC)
}

# modify column names for consistency for valid MATLAB identifiers
# change: p.adjust to pAdjust
aux <- colnames(results)
aux[c(5,6,7)] <- c("pValue", "pAdjust", "qValue")
colnames(results) <- aux

#write results to file: 
fileName <- paste(outputDir, "enrichment-GO.txt", sep = "")
cat(sprintf("writing results to file: %s\n", fileName))
write.table(results, fileName, sep = '\t', row.names = FALSE)


fileName <- paste(outputDir, 'enrichment-GO_BP.pdf', sep = "")
pdf(fileName)
barplot(enrichment_BP, showCategory=20)
dev.off()

fileName <- paste(outputDir, 'enrichment-GO_MF.pdf', sep = "")
pdf(fileName)
barplot(enrichment_MF, drop=TRUE, showCategory=20)
dev.off()

fileName <- paste(outputDir, 'enrichment-GO_CC.pdf', sep = "")
pdf(fileName)
barplot(enrichment_CC, showCategory=20)
dev.off()
