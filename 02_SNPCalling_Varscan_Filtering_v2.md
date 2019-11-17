--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------

######################################### FILTERING SNP_ALL_OUT WITH 09_filter_variants.R #########################################
## filtering raw variants

#install.packages('stringr')
library(stringr)

dat <- read.table("~/RESEARCH/SpOA/combined_data/analysis/snp_all_out", stringsAsFactors=FALSE, skip=1)

> nrow(dat)
[1] 19529443

datnm <- read.table("~/RESEARCH/SpOA/combined_data/analysis/snp_all_out", stringsAsFactors=FALSE, nrows=1)

pops <- c(
                "D1_7_0_S_02", "D1_7_0_S_03", "D1_7_0_S_04", "D1_7_0_S_05",
                "D1_7_0_S_08", "D1_7_5_S_01", "D1_7_5_S_03", "D1_7_5_S_04",
                "D1_7_5_S_05", "D1_7_5_S_08", "D1_8_1_S_01", "D1_8_1_S_03",
                "D1_8_1_S_05", "D1_8_1_S_08", "D7_7_0_S_25", "D7_7_0_S_26",
                "D7_7_0_S_27", "D7_7_0_S_28", "D7_7_0_S_29", "D7_7_0_S_30",
                "D7_7_0_V_02", "D7_7_0_V_03", "D7_7_0_V_04", "D7_7_0_V_05",
                "D7_7_0_V_06", "D7_7_5_S_31", "D7_7_5_S_32", "D7_7_5_S_33",
                "D7_7_5_S_36", "D7_7_5_V_07", "D7_7_5_V_08", "D7_7_5_V_09",
                "D7_7_5_V_10", "D7_7_5_V_12", "D7_8_1_S_13", "D7_8_1_S_14",
                "D7_8_1_S_15", "D7_8_1_S_18", "D7_8_1_V_19", "D7_8_1_V_20",
                "D7_8_1_V_21", "D7_8_1_V_22", "D7_8_1_V_23", "D7_8_1_V_24")

colnames(dat) <- c(datnm[1,1:10], pops)

# first, remove all where the number of samples not covered/called is 0
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1) #87403



--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
### NO LONGER DOING THIS, REID SAYS THERE IS AN ISSUE WITH TRYING TO CALL HOMOZYGOTES THIS WAY, TO INSTEAD
### FILTER BY ALLELE FREQUENCIES

# now, filter to include only variable sites
#    # For 44 samples, will set this to 41 (min of 3/44 variable sites - maf)
#    # so homozygotes both need to be < 41
#dat1 <- dat1[which(dat1$SamplesRef < 41 & dat1$SamplesHom < 41 ),]
#nrow(dat1)
#[1] 59,645


--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------


dat2 <- dat1

# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3***

# filt: 21942


dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
[1] 85,899



# Reid NOTES:
#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 40
# from the manual: Also, VarScan reports variants on a biallelic basis.
    #That is, for a given SNP call, the "reads1" column is the number of
    #reference-supporting reads (RD), and the "reads2" column is the number of
    #variant-supporting reads (AD).
    #There may be additional reads at that position showing other bases (SNP or indel variants).
    #If these other variants meet the calling criteria, they will be reported in
    #their own line. If not, it may look like you have "missing" reads.
# columns for each call are:
    #consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.



keep <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(keep) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])
        # sum up reads
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

#low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 40)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 36,318

#low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 35)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 54,435

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 30)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 84,395   # USING THIS ONE


dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 84,395

########## Checking out updates in Reid's code: ##############################################


##################################################################################################################
# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]
nrow(dat4)
# [1] 411702
dat3 <- dat4
##################################################################################################################


# here calculate allele freqs
# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
af <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(af) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        af[,grep(i_pop, colnames(af))] <- maj/(maj+ minor)

    }

##################################################################################################################
    # get rid of invariant sites
    dat4 <- dat3[(which(rowSums(af) > 0)),]
    dat3 <- dat4
    af <- af[(which(rowSums(af) > 0)),]
    nrow(dat3)
    #[1] 411379

    af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))

##################################################################################################################


afct.maf <- (sapply(af,function(x)
          ifelse(x > 0.5, (1-x), x)))

low_maf <- (apply(afct.maf, 1, function(x) {ifelse((length(which(x > 0.025)) < 4), FALSE, TRUE)}))
sum(low_maf)
# [1] 54,427



dat4 <- dat3[low_maf,]
nrow(dat4)
#[1] *54,427* SNPs

af_f <- af[low_maf,]

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af_f))

colnames(af.out) <- c("SNP", colnames(af_f))



# save filtered genotypes


write.table(file="~/RESEARCH/SpOA/combined_data/analysis/filtered_variants_v2.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/RESEARCH/SpOA/combined_data/analysis/filtered_allele_freqs_v2.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)
