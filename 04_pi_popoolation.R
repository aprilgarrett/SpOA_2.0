# Groups:
D1_control.pi
D7_control.pi
D7_7_5S.pi
D7_7_0S.pi
D7_7_5V.pi
D7_7_0V.pi
D7_8_1V.pi


# pi columns are: chrom, window position (mean val), number of snps in window, fraction of site with sufficient cov, pi measure for window

D1_control <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D1_control.pi", header=FALSE,stringsAsFactors=FALSE)
D1_control <- D1_control[which(D1_control$V5 != "na"),]
colnames(D1_control) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D1_control$pi))
#0.01383994

D7_control <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_control.pi", header=FALSE,stringsAsFactors=FALSE)
D7_control <- D7_control[which(D7_control$V5 != "na"),]
colnames(D7_control) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_control$pi))
#0.01433202

D7_7_5S <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.pi", header=FALSE,stringsAsFactors=FALSE)
D7_7_5S <- D7_7_5S[which(D7_7_5S$V5 != "na"),]
colnames(D7_7_5S) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_7_5S$pi))
#0.01355107

D7_7_0S <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.pi", header=FALSE,stringsAsFactors=FALSE)
D7_7_0S <- D7_7_0S[which(D7_7_0S$V5 != "na"),]
colnames(D7_7_0S) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_7_0S$pi))
#0.01385727

D7_7_5V <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_7_5V.pi", header=FALSE,stringsAsFactors=FALSE)
D7_7_5V <- D7_7_5V[which(D7_7_5V$V5 != "na"),]
colnames(D7_7_5V) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_7_5V$pi))
#0.01403871

D7_7_0V <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_7_0V.pi", header=FALSE,stringsAsFactors=FALSE)
D7_7_0V <- D7_7_0V[which(D7_7_0V$V5 != "na"),]
colnames(D7_7_0V) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_7_0V$pi))
#0.0133569

D7_8_1V <- read.table("~/RESEARCH/SpOA/combined_data/mapped/D7_8_1V.pi", header=FALSE,stringsAsFactors=FALSE)
D7_8_1V <- D7_8_1V[which(D7_8_1V$V5 != "na"),]
colnames(D7_8_1V) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D7_8_1V$pi))
#0.01425859
