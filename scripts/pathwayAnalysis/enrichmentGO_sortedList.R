# enrichment go 
# features are already sorted by log2FoldChange

# It requieres R Packages:
# clusterProfiler, https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# BBmisc
# GO.db
# org.Hs.eg.db

# TEST:
# with file: 
# path = '~/Documents/Projects/EscapePilot/2017-01-31_Gene-Expression_and_Metabolome/VotingAlgorithm/'
# vote = 3
# UpDown = 'upregulated'   or UpDown = 'downregulated'
# dataFile = paste(path, 'FS-voting-lambda0000001_vote-', vote, '_', UpDown, '.txt', sep = "")
# outputDir = '~/Documents/Projects/EscapePilot/2017-01-31_Gene-Expression_and_Metabolome/VotingAlgorithm/enrichGO/'
# Run in terminal:
#  Rscript name-R-script.R 'file-name_List-of-genes' output-directory
#  Rscript enrichmentGO_sortedList.R '~/Documents/Projects/EscapePilot/2017-01-31_Gene-Expression_and_Metabolome/VotingAlgorithm/FS-voting-lambda0000001_vote-3_downregulated.txt' output-folder/

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

D = data
entrez<-bitr(D$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE)
ranked = entrez  

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
fileName <- paste(outputDir, "enrichment-GO-Voting-", vote, "_", UpDown, ".txt", sep = "")
cat(sprintf("writing results to file: %s\n", fileName))
write.table(results, fileName, sep = '\t', row.names = FALSE)


fileName <- paste(outputDir, 'enrichment-GO_BP-Voting-', vote, '_', UpDown, '.pdf', sep = "")
pdf(fileName)
barplot(enrichment_BP, showCategory=20)
dev.off()

fileName <- paste(outputDir, 'enrichment-GO_MF-Voting-', vote, '_', UpDown, '.pdf', sep = "")
pdf(fileName)
barplot(enrichment_MF, drop=TRUE, showCategory=20)
dev.off()

fileName <- paste(outputDir, 'enrichment-GO_CC-Voting-', vote, '_', UpDown, '.pdf', sep = "")
pdf(fileName)
barplot(enrichment_CC, showCategory=20)
dev.off()
