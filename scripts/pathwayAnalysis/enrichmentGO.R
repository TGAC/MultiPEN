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
  enrichment_BP <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary_BP <- summary(enrichment_BP)
  head(summary(enrichment_BP))
  if(nrow(enrichmentSummary_BP)>0){
    aux <- cbind(enrichmentSummary_BP, "BP")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary_BP <- aux        
  }  
  #add to results: enrichment for BP category
  results <- enrichmentSummary_BP
  
  #Enrichment for subontology MF (Molecular Function)
  subclassOnt <- "MF"
  enrichment_MF <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary_MF <- summary(enrichment_MF)
  head(summary(enrichment_MF)) 
  if(nrow(enrichmentSummary_MF)>0){
    aux <- cbind(enrichmentSummary_MF, "MF")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary_MF <- aux       
    #add to results: enrichment for BP category
    results <- rbind(results, enrichmentSummary_MF)
  }  
  
  
  #Enrichment for subclass CC (Cellular Cellular Component)
  subclassOnt <- "CC"
  enrichment_CC <- enrichGO(ranked$ENTREZID, organism='human', ont=subclassOnt, readable=TRUE)
  enrichmentSummary_CC <- summary(enrichment_CC)
  head(summary(enrichment_CC))
  if(nrow(enrichmentSummary_CC)>0){
    aux <- cbind(enrichmentSummary_CC, "CC")
    colnames(aux)[10]<- 'subontology'
    enrichmentSummary_CC <- aux
    #add to results: enrichment for BP category
    results <- rbind(results, enrichmentSummary_CC)
  }
  
  
  #if(nrow(results)>0){
    #write results to file: 
    fileName <- paste(outputDir, "enrichGO_subcategories.txt", sep = "")
    cat(sprintf("writing results to file: %s\n", fileName))
    write.table(results, fileName, sep = '\t', row.names = FALSE)
  #}
  
  fileName <- paste(outputDir, 'enrichmentGO.pdf', sep = "")
  pdf(fileName)
  par(mfcol=c(1,3))
  barplot(enrichment_BP, drop=TRUE, showCategory=20)
  barplot(enrichment_MF, drop=TRUE, showCategory=20)
  barplot(enrichment_CC, drop=TRUE, showCategory=20)
  dev.off()
  