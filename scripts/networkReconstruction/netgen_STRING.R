###################################################
##find interactions for small number of proteins/genes
library(STRINGdb)
args <- commandArgs(trailingOnly = TRUE)

##options
#input file with gene names
file<-args[1]
file<-"test.txt"
n<-substr(file,1,nchar(file)-4)

#species identifier, uniprot taxonomy
sp<-as.numeric(args[2])
sp<-as.numeric("9606")
#oryzias latipes: 8090
#danio rerio: 7955
#homo sapiens: 9606
#mus musculus: 10090

#STRING db version
sdb<-as.numeric(args[3])

##map content of STRING db
#load STRING
string_db <- STRINGdb$new( version=sdb, species=sp, score_threshold=0, input_directory="" )
#read gene list
v <- read.delim(file, header=T,sep = "\t", stringsAsFactors=FALSE)
#map
mapped <- string_db$map( v,  "V1", removeUnmappedRows = TRUE )
# get the best 200 hits
#hits = mapped$STRING_id[1:200]
# plot the STRING network png
#string_db$plot_network( hits )
#string_db$get_png( hits )

#get interactions 
inter<-string_db$get_interactions(mapped$STRING_id)
#annotate source and target nodes
p<-paste(sp,".",sep="")
from <- gsub(p,"",inter$from)
to <- gsub(p,"",inter$to)
#from <- string_db$add_proteins_description(from)
#colnames(from) <- c("STRING_id_1", "preferred_name_1", "protein_size_1", "annotation_1")
#to <- string_db$add_proteins_description(to)
#colnames(to) <- c("STRING_id_2", "preferred_name_2", "protein_size_2", "annotation_2")
#from_gs <- getBM(attributes = c("mgi_symbol"), filters = "ensembl_peptide_id", values = from, mart = mouse, uniqueRows = FALSE)

#from_gs<- queryMany(from, scopes="ensemblprotein", fields="symbol", species="mouse")
#to_gs<- queryMany(to, scopes="ensemblprotein", fields="symbol", species="mouse")
rest<-inter[3:16]
interact<-cbind(from,to,inter[3:16])
net<-cbind(from,to,inter[16]/1000)

#interact<-cbind(from,to,rest)
write.table(interact,paste("interactome.",args[1],sep="_"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(net,paste("SI_network",args[1],sep="_"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

