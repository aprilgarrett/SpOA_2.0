af <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_allele_freqs_v2.txt", header=TRUE)
dat3 <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_variants_v2.txt", header=TRUE)

library(stringr)

#A1 <- dat3[11:ncol(dat3)]
#A1[] <- lapply(A1, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,3]) })
#A2 <- dat3[11:ncol(dat3)]
#A2[] <- lapply(A2, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,4]) })

DP1 <- dat3[5]
DP1[] <- lapply(DP1, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,3]) })
names(DP1) <- "DP1"

DP2 <- dat3[5]
DP2[] <- lapply(DP2, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,4]) })
names(DP2) <- "DP2"

DPtotal <- dat3[5]
DPtotal[] <- lapply(DPtotal, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,2]) })
names(DPtotal) <- "DPtotal"

alldp <- as.data.frame(c(DP1, DP2, DPtotal))

write.table(file="spoa_varscan_v2_depths.txt", alldp, sep="\t", row.names=FALSE,quote=FALSE)