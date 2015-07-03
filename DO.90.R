library(GenomicRanges)
library(sqldf)
library(data.table)

gene <- read.delim("gene_location_info.txt")
gene$chrom <- as.character(gene$chrom)
gene[which(gene$chrom == 23),"chrom"] = "X"
gene[which(gene$chrom == 24),"chrom"] = "Y"

gene.gr <- with(gene, GRanges(chrom, IRanges(start, end, names=Gene)))

get_copy_number <- function (sample) {
  file = paste("OutputFolder/Results/", sample, "/FastCallResults_", sample, ".txt", sep="")
  # FastCallResults <- read.delim("OutputFolder/Results/05427188T/FastCallResults_05427188T.txt")
  FastCallResults <- read.delim(file=file)
  FastCallResults.gr <- with(FastCallResults,
                                       GRanges(Chromosome, IRanges(Start, End, names=CN)))
  
  overlaps.all.genes <- findOverlaps(gene.gr, FastCallResults.gr, type="within", select="all")
  
  fastcall.subset <- FastCallResults[overlaps.all.genes@subjectHits, 
                          c("Chromosome", "Start", "End", "Segment", "CNF", "CN", "Call", "ProbCall")]
  row.names(fastcall.subset) <- NULL
  
  gene.subset <- gene[overlaps.all.genes@queryHits, 
                      c("Gene", "chrom", "start", "end")]
  row.names(gene.subset) <- NULL
  
  rm(overlaps.all.genes)
  
  CNV = data.table(cbind(gene.subset, fastcall.subset)) # data.table function in data.table library
  return(CNV)
}

samples = c("05427188T", "07945734T", "12879057T", "17519732T", "24588763T", "24929834T", "25913317T", "25938404T", "28573150T", "30076971T", "30261085T", "30657633T")

CNV = NULL
for (sample in samples) {
  print (sample)
  CNV = rbind(CNV, cbind(sample=sample, as.data.frame(get_copy_number(sample))))
  
  #head(cbind(sample=sample, CNV))
  #head(cbind(sample=sample, CNV[,c("Gene", "Segment", "CNF", "CN", "Call", "ProbCall")]))
  
}

library(reshape2)
library(dplyr)

## CNV.summary.call
CNV.summary.Call <- as.data.frame(dcast(CNV, Gene ~ sample, value.var="Call", fill=0))

#rowsum <- apply(CNV.summary.Call[,-1], 1, function(x) length(which(!is.na(x))))
rowsum <- rowSums(CNV.summary.Call[,-1])
CNV.summary.Call <- cbind(CNV.summary.Call, rowsum=rowsum)
CNV.summary.Call <- CNV.summary.Call[order(-CNV.summary.Call$rowsum),]
#CNV.summary.Call <- arrange(CNV.summary.Call, desc(rowsum))

write.csv(CNV.summary.Call, "CNV.summary.Call.csv", quote=F, row.names=F)

## CNV.summary.CN
CNV.summary.CN <- as.data.frame(dcast(CNV, Gene ~ sample, value.var="CN", fill=2))

#rowsum <- apply(CNV.summary.CN[,-1], 1, function(x) length(which(!is.na(x))))
rowsum <- rowSums(CNV.summary.CN[,-1])
CNV.summary.CN <- cbind(CNV.summary.CN, rowsum=rowsum)
CNV.summary.CN <- CNV.summary.CN[order(-CNV.summary.Call$rowsum),]
#CNV.summary.CN <- arrange(CNV.summary.CN, desc(rowsum))

write.csv(CNV.summary.CN, "CNV.summary.CN.csv", quote=F, row.names=F)



mutated.gene <- c("MUC2", "KRT8", "PIK3CA", "TNRC6A", "ABHD17A", "AHNAK2", "ARHGAP35", "DCD27", "DCSTAMP", "FLG", "HERC2", "MST1", "MYO5B", "PRKRIR", "RECQL4", "SPOP", "XIRP2")
View(CNV.summary.Call[which(CNV.summary.Call$Gene %in% mutated.gene),])
