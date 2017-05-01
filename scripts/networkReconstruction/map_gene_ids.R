#map list of gene IDs
library("mygene")
args <- commandArgs(trailingOnly = TRUE)

#options
file<-args[1]
#file<-"genes.txt"
n<-substr(file,1,nchar(file)-4)

#recognised ID types
#symbol; entrezgene; ensemblgene; ensemblprotein; ensembltranscript
#for full list of fields see: http://mygene.info/doc/query_service.html#available-fields
from <- args[2]
#from <- "symbol"
to <- args[3]
#to <- "ensembl.gene"

#define species NCBI taxonomy ID
spec<-args[4]
#spec<-"9606"

genes <- read.delim(file, header=F, sep = "\t", stringsAsFactors=FALSE,fill=T)
ids<-as.list(genes)

out <- queryMany(ids, scopes=from, fields=to, species=spec)
m<-out@listData$ensembl
p<-lapply(m,function(x) unlist(x))
new<-data.frame(V1 = rep(out$query, sapply(p, length)), V2 = unlist(p))
colnames(new)<-c(from,to)
write.table(new, paste(n,from,to,"IDmapped.txt", sep='-'),sep="\t",row.names=F,quote = FALSE)



