#!/usr/bin/env Rscript
SS50M <- cbind(read.delim("SureSelect_Human_All_Exon_50M.bed", header=F), V4=0)
#SS50M$V1 = gsub("chr", "", SS50M$V1)
write.table(SS50M, "SS50M.bed", row.names=F, col.names=F, quote=F, sep="\t")
