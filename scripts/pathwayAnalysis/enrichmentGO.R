  args = commandArgs(trailingOnly=TRUE)

  # test if there is at least one argument: if not, return an error
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    args[2] = "output/"
  }
  
  dataFile <- args[1]
  cat(sprintf("Loading data from file: %s\n", dataFile))

  outputDir <- args[2]
  if(!dir.exists('output')){
    cat(sprintf("Creating output directory: %s\n", outputDir))
    dir.create("output/")
  }


#GSEAnalysis <- function(dataFile){
  # example data:
  # dataFile <- '~/Documents/Projects/multipen/2016-08-16_MultiPEN_v002/testing_MultiPEN/ExampleOutputs/MultiPEN-Rankings_lambda0.1.txt'
  getwd()
  #path <- '/usr/users/ga004/troncosp/apps/multipen/pathwayAnalysis/'
  #setwd(path)
  #dataFile <- "MultiPEN-Rankings_lambda0.1.txt"
  data <- read.table(dataFile, header=TRUE, sep = "\t", stringsAsFactors=FALSE)
  cat(sprintf("Number of features: %i\n",nrow(data)))


  library(clusterProfiler)
  library(BBmisc)
  library(GO.db)

  D <- sortByCol(data, 'ranking')
  D <- D[D[,2]!=0,]
  D <- D[,c(1,3)]  # Only interested in the first three columns [name, ranking]
  entrez<-bitr(D$name, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db",drop = FALSE)
  ranked<-merge(D,entrez,by.x='name',by.y='SYMBOL')

  cat(sprintf("Performing over-representation analysis (enrichGO) ...  "))
  cat(sprintf("Results saved to folder: %s\n", outputDir))
  
  #Enrichment for subontology BP (Biological Processes)
  subclassOnt <- "BP"
  enrichment <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary <- summary(enrichment)
  if(nrow(enrichmentSummary)>0){
    aux <- cbind(enrichmentSummary, "BP")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary <- aux
    head(summary(enrichment))
    fileName <- paste(outputDir, 'enrichment_BP.pdf', sep = "")
    pdf(fileName)
    barplot(enrichment, drop=TRUE, showCategory=20)
    dev.off()
  }  
  #add to results: enrichment for BP category
  results <- enrichmentSummary
  
  #Enrichment for subontology MF (Molecular Function)
  subclassOnt <- "MF"
  enrichment <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary <- summary(enrichment)
  if(nrow(enrichmentSummary)>0){
    aux <- cbind(enrichmentSummary, "MF")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary <- aux
    head(summary(enrichment))
    fileName <- paste(outputDir, 'enrichment_MF.pdf', sep = "")
    pdf(fileName)
    barplot(enrichment, drop=TRUE, showCategory=20)
    dev.off()
    #add to results: enrichment for BP category
    results <- rbind(results, enrichmentSummary)
  }  
  
  
  #Enrichment for subclass CC (Cellular Cellular Component)
  subclassOnt <- "CC"
  enrichment <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary <- summary(enrichment)
  if(nrow(enrichmentSummary)>0){
    aux <- cbind(enrichmentSummary, "CC")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary <- aux
    head(summary(enrichment))
    fileName <- paste(outputDir, 'enrichment_CC.pdf', sep = "")
    pdf(fileName)
    barplot(enrichment, drop=TRUE, showCategory=20)
    dev.off()
    #add to results: enrichment for BP category
    results <- rbind(results, enrichmentSummary)
  }
  
  
  if(nrow(results)>0){
    #write results to file: 
    fileName <- paste(outputDir, "enrichGO_subcategories.txt", sep = "")
    cat(sprintf("writing results to file: %s\n", fileName))
    write.table(results, fileName, sep = '\t', row.names = FALSE)
  }
  