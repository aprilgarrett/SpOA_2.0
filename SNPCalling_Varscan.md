# Remaking vcf with Reid's script (from UrchinMS)

tail -n +2 /data/OASV2/RG_7601_Probe.extended.bed | cut -f 1 | sort | uniq| grep -v "#scaffold" | xargs -I {} -n 1 -P 10 sh -c "/data/programs/samtools-1.4.1/samtools mpileup -u --skip-indels -d 10000  -t DP,AD,INFO/AD -f /data/OASV2/Spur_3.1.LinearScaffold.fa -r {} ~/RESEARCH/SpOA/combined_data/mapped/complete.merged.sorted.bam | /data/programs/bcftools/bin/bcftools call -mv - > ~/RESEARCH/SpOA/combined_data/variants/temp/tmp.{}.vcf"


### PUTTING THIS ON HOLD FOR NOW

--------------------------------------------------------------------------------------------------------
# USING REID'S SCRIPTS FROM HIS A. TONSA WORK

08_varscan_all.sh (very liberal in SNP calls, need additional filtering) --> going to try making VCF from this and continue with urchinMS code, if that
doesn't work then will go back to this and continue to 09_filter_variants.R after making a snp_all_out file from Varscan.
09_filter_variants.R
snp_analysis.R


##!/bin/bash -l
#java -jar ~/bin/VarScan.v2.3.9.jar mpileup2snp ~/reciprocal_t/analysis/merged.pileup --min-coverage 30 --min-reads 2 --min-avg-qual 20 --min-var-freq 0.01 --p-value 0.1 > ~/reciprocal_t/analysis/snp_out


• From Reid's 08_varscan_all.sh

# varscan_v1.sh
# Apparently varscan doesn't read one input
```
#!/bin/bash -l

cd ~/RESEARCH/SpOA/combined_data/mapped

/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f ~/RESEARCH/OASV2/Spur_3.1.LinearScaffold.fa ~/RESEARCH/SpOA/combined_data/mapped/complete.merged.sorted.bam | java -jar /users/r/b/rbrennan/bin/VarScan.v2.3.9.jar mpileup2snp --output-vcf 1 --min-coverage 30 --min-reads 10 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/RESEARCH/SpOA/combined_data/variants/varscan_spoa_v1.vcf

# 48 files that Reid used

# /users/r/b/rbrennan/bin/samtools-1.6/samtools ### Version Reid used
```

### The vcf file created above only recognized 1 sample --> need to input individual bam files, so need to merge and sort data and Newdata runs together so that I have 44 samples total (Merging_data_NewData_1-4.sh files).
--------------------------------------------------------------------------------------------------------

• From Reid's 08_varscan_all.sh (RE-DO)

```
#!/bin/bash -l

cd ~/RESEARCH/SpOA/combined_data/mapped

/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f ~/RESEARCH/OASV2/Spur_3.1.LinearScaffold.fa SpOA_DNA_D1_7_0_S_02.sorted.bam SpOA_DNA_D1_7_0_S_03.sorted.bam SpOA_DNA_D1_7_0_S_04.sorted.bam SpOA_DNA_D1_7_0_S_05.sorted.bam SpOA_DNA_D1_7_0_S_08.sorted.bam SpOA_DNA_D1_7_5_S_01.sorted.bam SpOA_DNA_D1_7_5_S_03.sorted.bam SpOA_DNA_D1_7_5_S_04.sorted.bam SpOA_DNA_D1_7_5_S_05.sorted.bam SpOA_DNA_D1_7_5_S_08.sorted.bam SpOA_DNA_D1_8_1_S_01.sorted.bam SpOA_DNA_D1_8_1_S_03.sorted.bam SpOA_DNA_D1_8_1_S_05.sorted.bam SpOA_DNA_D1_8_1_S_08.sorted.bam SpOA_DNA_D7_7_0_S_25.sorted.bam SpOA_DNA_D7_7_0_S_26.sorted.bam SpOA_DNA_D7_7_0_S_27.sorted.bam SpOA_DNA_D7_7_0_S_28.sorted.bam SpOA_DNA_D7_7_0_S_29.sorted.bam SpOA_DNA_D7_7_0_S_30.sorted.bam SpOA_DNA_D7_7_0_V_02.sorted.bam SpOA_DNA_D7_7_0_V_03.sorted.bam SpOA_DNA_D7_7_0_V_04.sorted.bam SpOA_DNA_D7_7_0_V_05.sorted.bam SpOA_DNA_D7_7_0_V_06.sorted.bam SpOA_DNA_D7_7_5_S_31.sorted.bam SpOA_DNA_D7_7_5_S_32.sorted.bam SpOA_DNA_D7_7_5_S_33.sorted.bam SpOA_DNA_D7_7_5_S_36.sorted.bam SpOA_DNA_D7_7_5_V_07.sorted.bam SpOA_DNA_D7_7_5_V_08.sorted.bam SpOA_DNA_D7_7_5_V_09.sorted.bam SpOA_DNA_D7_7_5_V_10.sorted.bam SpOA_DNA_D7_7_5_V_12.sorted.bam SpOA_DNA_D7_8_1_S_13.sorted.bam SpOA_DNA_D7_8_1_S_14.sorted.bam SpOA_DNA_D7_8_1_S_15.sorted.bam SpOA_DNA_D7_8_1_S_18.sorted.bam SpOA_DNA_D7_8_1_V_19.sorted.bam SpOA_DNA_D7_8_1_V_20.sorted.bam SpOA_DNA_D7_8_1_V_21.sorted.bam SpOA_DNA_D7_8_1_V_22.sorted.bam SpOA_DNA_D7_8_1_V_23.sorted.bam SpOA_DNA_D7_8_1_V_24.sorted.bam | java -jar /users/r/b/rbrennan/bin/VarScan.v2.3.9.jar mpileup2snp --output-vcf 1 --min-coverage 30 --min-reads 10 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/RESEARCH/SpOA/combined_data/variants/varscan_spoa_v2.vcf

```

--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
### NOTE :  Having issues with VCF, going to try just using Reid's scripts from: 09_filter_variants.R & snp_analysis.R

```
#!/bin/bash -l

cd ~/RESEARCH/SpOA/combined_data/mapped

/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup -Q 20 -B --max-depth 3000 --skip-indels -f ~/RESEARCH/OASV2/Spur_3.1.LinearScaffold.fa SpOA_DNA_D1_7_0_S_02.sorted.bam SpOA_DNA_D1_7_0_S_03.sorted.bam SpOA_DNA_D1_7_0_S_04.sorted.bam SpOA_DNA_D1_7_0_S_05.sorted.bam SpOA_DNA_D1_7_0_S_08.sorted.bam SpOA_DNA_D1_7_5_S_01.sorted.bam SpOA_DNA_D1_7_5_S_03.sorted.bam SpOA_DNA_D1_7_5_S_04.sorted.bam SpOA_DNA_D1_7_5_S_05.sorted.bam SpOA_DNA_D1_7_5_S_08.sorted.bam SpOA_DNA_D1_8_1_S_01.sorted.bam SpOA_DNA_D1_8_1_S_03.sorted.bam SpOA_DNA_D1_8_1_S_05.sorted.bam SpOA_DNA_D1_8_1_S_08.sorted.bam SpOA_DNA_D7_7_0_S_25.sorted.bam SpOA_DNA_D7_7_0_S_26.sorted.bam SpOA_DNA_D7_7_0_S_27.sorted.bam SpOA_DNA_D7_7_0_S_28.sorted.bam SpOA_DNA_D7_7_0_S_29.sorted.bam SpOA_DNA_D7_7_0_S_30.sorted.bam SpOA_DNA_D7_7_0_V_02.sorted.bam SpOA_DNA_D7_7_0_V_03.sorted.bam SpOA_DNA_D7_7_0_V_04.sorted.bam SpOA_DNA_D7_7_0_V_05.sorted.bam SpOA_DNA_D7_7_0_V_06.sorted.bam SpOA_DNA_D7_7_5_S_31.sorted.bam SpOA_DNA_D7_7_5_S_32.sorted.bam SpOA_DNA_D7_7_5_S_33.sorted.bam SpOA_DNA_D7_7_5_S_36.sorted.bam SpOA_DNA_D7_7_5_V_07.sorted.bam SpOA_DNA_D7_7_5_V_08.sorted.bam SpOA_DNA_D7_7_5_V_09.sorted.bam SpOA_DNA_D7_7_5_V_10.sorted.bam SpOA_DNA_D7_7_5_V_12.sorted.bam SpOA_DNA_D7_8_1_S_13.sorted.bam SpOA_DNA_D7_8_1_S_14.sorted.bam SpOA_DNA_D7_8_1_S_15.sorted.bam SpOA_DNA_D7_8_1_S_18.sorted.bam SpOA_DNA_D7_8_1_V_19.sorted.bam SpOA_DNA_D7_8_1_V_20.sorted.bam SpOA_DNA_D7_8_1_V_21.sorted.bam SpOA_DNA_D7_8_1_V_22.sorted.bam SpOA_DNA_D7_8_1_V_23.sorted.bam SpOA_DNA_D7_8_1_V_24.sorted.bam | java -jar /users/r/b/rbrennan/bin/VarScan.v2.3.9.jar mpileup2snp --mpileup 1 --min-coverage 30 --min-reads 10 --min-avg-qual 20 --min-var-freq 0.01 --variants --p-value 0.1 > ~/RESEARCH/SpOA/combined_data/analysis/snp_all_out

```

--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------

######################################### FILTERING SNP_ALL_OUT WITH 09_filter_variants.R #########################################
## filtering raw variants


library(stringr)

dat <- read.table("~/RESEARCH/SpOA/combined_data/analysis/snp_all_out", stringsAsFactors=FALSE, skip=1)

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

# now, filter to include only variable sites
    # For 44 samples, will set this to 41 (min of 3/44 variable sites - maf)
    # so homozygotes both need to be < 41
dat1 <- dat1[which(dat1$SamplesRef < 41 & dat1$SamplesHom < 41 ),]
nrow(dat1)
#[1] 59,645
dat2 <- dat1

# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3

dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
[1] 58,630

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

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 41)) > 0), FALSE, TRUE)}))
sum(low_cv)
#[1] 22491
dat4 <- dat3[low_cv,]
nrow(dat4)
dat3 <- dat4
nrow(dat3)
# [1] 22491

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

afct.maf <- (sapply(af,function(x)  
          ifelse(x > 0.5, (1-x), x)))

low_maf <- (apply(afct.maf, 1, function(x) {ifelse((length(which(x > 0.025)) < 4), FALSE, TRUE)}))
sum(low_maf)
# [1] 18089

dat4 <- dat3[low_maf,]
nrow(dat4)
#[1] 18089

af_f <- af[low_maf,]

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af_f))

colnames(af.out) <- c("SNP", colnames(af_f))



# save filtered genotypes


write.table(file="~/RESEARCH/SpOA/combined_data/analysis/filtered_variants.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/RESEARCH/SpOA/combined_data/analysis/filtered_allele_freqs.txt", af.out, sep="\t",
              row.names=FALSE, quote=FALSE)
