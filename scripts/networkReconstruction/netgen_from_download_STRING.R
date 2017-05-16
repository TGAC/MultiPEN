#generation of PPI network with files downloaded from STRING server
args <- commandArgs(trailingOnly = TRUE)
#species
spec<-args[3]
#spec<-"9606"

#read data
pp <- read.delim(args[1], header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
#pp <- read.delim("STRING_prot2prot.txt", header=F, sep = " ", stringsAsFactors=F,fill=T)
colnames(pp)<-c("protS","protT","score")
gids <- read.delim(args[2], header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
#gids <- read.delim("HGNC-9606.protein.aliases.v10.txt", header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
colnames(gids)<-c("string_protein_id","HGNC","sourcedb")

#prot 2 prot
#hgnc
p_mapped <- merge(pp,gids,by.x='protS',by.y='string_protein_id',all.x=T)
p_mapped2 <- merge(p_mapped,gids,by.x='protT',by.y='string_protein_id',all.x=T)
p2p_hgnc<-subset(p_mapped2, select=c("HGNC.x", "HGNC.y", "score"))
colnames(p2p_hgnc)<-c("source","target","score")
p2p_hgnc$score<-p2p_hgnc$score/1000
write.table(p2p_hgnc,"PPI_STRING_hgnc.txt",sep="\t",row.names=F,quote = FALSE)
