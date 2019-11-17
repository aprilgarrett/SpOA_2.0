# go_format.R

# Rscript ~/RESEARCH/SpOA/scripts/go_format_spoa.R 2> ~/RESEARCH/SpOA/log_out/go_format_spoa.stderr_$(date +"%F_%R").txt 1> ~/RESEARCH/SpOA/log_out/go_format_spoa.stdout_$(date +"%F_%R").txt

# all variants
dat <- read.delim("~/RESEARCH/SpOA/combined_data/analysis/cmh.master.sort.out.2", header=TRUE, stringsAsFactor=FALSE, sep= "\t")
dat$logP_75 <- -log10(dat$D7_7_5_S_selection_pval)
dat$logP_70 <- -log10(dat$D7_7_0_S_selection_pval)
#dat$sig_pH75 <- gsub(" TRUE", "TRUE",dat$D7_7_5_S_SIG)
#dat$sig_pH70 <- gsub(" TRUE", "TRUE",dat$D7_7_0_S_SIG)

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
  a <- dat[which(dat$SPU_1 == gene[i]),]
  out[i,] <- a[which.min(as.numeric(a$D7_7_5_S_selection_pval)),]
}

write.table(file="~/RESEARCH/SpOA/combined_data/analysis/go_enrichment/cmh.pH7_5S_master_GO.out", out, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
  a <- dat[which(dat$SPU_1 == gene[i]),]
  out[i,] <- a[which.min(as.numeric(a$D7_7_0_S_selection_pval)),]
}

write.table(file="~/RESEARCH/SpOA/combined_data/analysis/go_enrichment/cmh.pH7_0S_master_GO.out", out, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")


#### overlap between the two ph
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
out$overlap_sig <- FALSE

for(i in 1:length(gene)){
  a <- dat[which(dat$SPU_1 == gene[i]),]
  if(sum(a$D7_7_5_S_SIG == TRUE) >= 1 & sum(a$D7_7_0_S_SIG == TRUE) >= 1){
    out[i,] <- a[which.min(as.numeric(a$D7_7_0_S_selection_pval)),]
    out$overlap_sig[i] <- TRUE
  }
  if(sum(a$D7_7_5_S_SIG == TRUE) >= 1 & sum(a$D7_7_0_S_SIG == TRUE) == 0){
    out[i,] <- a[which.min(as.numeric(a$D7_7_5_S_selection_pval)),]
    out$overlap_sig[i] <- FALSE
  }
  if(sum(a$D7_7_5_S_SIG == TRUE) == 0 & sum(a$D7_7_0_S_SIG == TRUE) >= 1){
    out[i,] <- a[which.min(as.numeric(a$D7_7_0_S_selection_pval)),]
    out$overlap_sig[i] <- FALSE
  }
  if(sum(a$D7_7_5_S_SIG == TRUE) == 0 & sum(a$D7_7_0_S_SIG == TRUE) == 0){
    out[i,] <- a[which.min(as.numeric(a$D7_7_5_S_selection_pval)),]
    out$overlap_sig[i] <- FALSE
  }
}

write.table(file="~/RESEARCH/SpOA/combined_data/analysis/go_enrichment/cmh.overlap_7_5S_7_0S_master_GO.out", out, col.names=TRUE,
            row.names=FALSE, quote=FALSE,sep="\t")