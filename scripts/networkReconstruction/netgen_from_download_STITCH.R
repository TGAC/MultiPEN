#generation of PPI network with
args <- commandArgs(trailingOnly = TRUE)
spec<-args[4]
#spec<-"9606"

#read data
chch <- read.delim(args[1], header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
#chch <- read.delim("test1.txt", header=F, sep = "\t", stringsAsFactors=FALSE,fill=T)
colnames(chch)<-c("chemS","chemT","score")
chp <- read.delim(args[2], header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
#chp <- read.delim("test2.txt", header=F, sep = "\t", stringsAsFactors=FALSE,fill=T)
colnames(chp)<-c("chemS","protT","score")
ids <- read.delim(args[3], header=T, sep = "\t", stringsAsFactors=FALSE,fill=T)
#ids <- read.delim("CHEBI-chemical.aliases.v5.0.tsv", header=F, sep = "\t", stringsAsFactors=FALSE,fill=T)
colnames(ids)<-c("flat","stereo","ChEBI","source")

#map ChEBI IDs
#chem 2 chem
ch1_mapped <- merge(chch,ids,by.x='chemS',by.y='stereo',all.x=T)
ch2_mapped <- merge(ch1_mapped,ids,by.x='chemT',by.y='stereo',all.x=T)
DF <- within(ch2_mapped, source <- ifelse(!is.na(ChEBI.x),ChEBI.x,chemS))
DF2 <- within(DF,target <- ifelse(!is.na(ChEBI.y),ChEBI.y,chemT))
ch2ch<-subset(DF2, select=c("source", "target", "score"))
write.table(ch2ch,"ch2ch_ChEBI.txt",sep="\t",row.names=F,quote = FALSE)

#chem 2 prot
ch_mapped <- merge(chp,ids,by.x='chemS',by.y='stereo',all.x=T)
DF <- within(ch_mapped, source <- ifelse(!is.na(ChEBI),ChEBI,chemS))
ch2p<-subset(DF, select=c("source", "protT", "score"))
p<-paste(spec,".",sep="")
ch2p$protT <- gsub(p,"",ch2p$protT)
write.table(ch2p,"ch2p_ChEBI.txt",sep="\t",row.names=F,quote = FALSE)

