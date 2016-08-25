# Read data from study on
# Fatty Liver Disease Tissue (Multi-Omics)

# Read transcriptomics data
pathInput <- "~/Documents/Public-Datasets/2016-8-24_Fatty-Liver-Disease_Gene-Expression_Metabolomics/2016-08-22_Data-Analysis/MultiPEN_input-data/"
fileName <- paste(pathInput, 'illuminaProbes.txt', sep = "")
data <- read.delim(fileName, sep = '\t', stringsAsFactors = FALSE)

# Install illuminaHumanv4 package
#source("https://bioconductor.org/biocLite.R")
#biocLite("illuminaHumanv4.db")

library("illuminaHumanv4.db")

# Pick one of the following:
# Map between Manufacturer IDs and ensembl gene accession number
x <- illuminaHumanv4ENSEMBL
# ## Map between Manufacturer IDs and Genes
# x <- illuminaHumanv4GENENAME
# ##  Map between Manufacturer Identifiers and Entrez Gene
# x <- illuminaHumanv4ENTREZID
# ## Map between Manufacturer IDs and Genes
# x <- illuminaHumanv4GENENAME


# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the GENENAME for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}


## Create table with two columns: entrezID and probeID
data[1:15,]
entrezID <- as.character(xx)
probeID <- names(xx)
annot <- data.frame(entrezID = entrezID, probeID = probeID)



####  MAP ENTREZID to Gene Names
library("org.Hs.eg.db")

## Bimap interface:
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

gSymbols <- as.character(xx)
aux <- names(xx)
annot2 <- data.frame(gSymbols = gSymbols, entrezID = aux)

annot_full <- merge(annot, annot2, by.x = "entrezID", by.y = "entrezID", all = FALSE)


#  Save table with columns:  entrezID, probeID, gSymbols
fileName <- paste(pathInput, "annot_entrezID_probeID_gSymbols.txt", sep = "")
write.table(annot_full, fileName, sep = "\t", row.names = FALSE)


# # Next try with:
# # Map between Manufacturer Identifiers and Gene Symbols
# x <- illuminaHumanv4SYMBOL
# # Get the probe identifiers that are mapped to a gene name
# mapped_probes <- mappedkeys(x)
# # Convert to a list
# xx <- as.list(x[mapped_probes])
# if(length(xx) > 0) {
#   # Get the GENENAME for the first five probes
#   xx[1:5]
#   # Get the first one
#   xx[[1]]
# }