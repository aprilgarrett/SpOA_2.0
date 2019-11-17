##########################################################################
##########################################################################
# CMH RE-DO (WITH GROUP) - 8/26/19
##########################################################################
##########################################################################

af <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_allele_freqs_v2.txt", header=TRUE)
dat3 <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_variants_v2.txt", header=TRUE)

# need count data for this. 
# get from dat3

library(stringr)

A1 <- dat3[11:ncol(dat3)]
A1[] <- lapply(A1, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,3]) })
A2 <- dat3[11:ncol(dat3)]
A2[] <- lapply(A2, function(x){as.numeric(str_split_fixed(x, ":", n=6)[,4]) })


##########################################################################

#install.packages('stringr')
library(stringr)

#install.packages("tidyr")
library(tidyr)
library(dplyr)


# get unique names of each treatment
treats <- as.data.frame(str_split_fixed(colnames(A1),n = 5,pattern = "_"))[,1:4] %>% unite("names", V1,V2,V3,V4)

uniq_treats <- unique(treats$names)
uniq_treats

# making data frames to save loop outputs to
A1.df <- data.frame(matrix(nrow=54427,ncol=0))
dim(A1.df)
A2.df <- data.frame(matrix(nrow=54427,ncol=0))
dim(A2.df)


# FOR LOOP TO RANDOMLY SAMPLE 
for(i in uniq_treats) {
  
  tmpdfA1 <- A1[,grep(i, (colnames(A1)))]
  tmpdfA2 <- A2[,grep(i, (colnames(A2)))]
  
  # sample randomly
  tmp_rs <- sample(1:ncol(tmpdfA1), size=4, replace=FALSE)
  
  rs_A1<-tmpdfA1[,tmp_rs] 
  
  rs_A2<-tmpdfA2[,tmp_rs]
  
  # pull each out from your tmp files above and save to new data frame
  A1.df <- cbind(A1.df,rs_A1)
  A2.df <- cbind(A2.df,rs_A2)
  
}

A1 <- A1.df
A2 <- A2.df

################################################################################################
### (1) cmh for D1-8.1S vs. D7-7.5S

D7_7_5_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$replicate <- substr(sub_all$ind, 1, 11)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_0_V|D7_7_5_V|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_7_5_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
  
}


D7_7_5_S_selection_pval[which(is.na(D7_7_5_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-7.5S cmh done")

cut_off <- (0.001)

af$D7_7_5_S_selection_pval <- D7_7_5_S_selection_pval
length(which(af$D7_7_5_S_selection_pval < cut_off))

################################################################################################

### (2) cmh for D1-8.1S vs. D7-7.0S

D7_7_0_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$replicate <- substr(sub_all$ind, 1, 11)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_5_S|D7_7_0_V|D7_7_5_V|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_7_0_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D7_7_0_S_selection_pval[which(is.na(D7_7_0_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-7.0S cmh done")

cut_off <- (0.001)

af$D7_7_0_S_selection_pval <- D7_7_0_S_selection_pval

length(which(af$D7_7_0_S_selection_pval < cut_off))


################################################################################################

### (3) cmh for D1-8.1S vs. D7-7.5V

D7_7_5_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_0_V|D7_7_5_S|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_7_5_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D7_7_5_V_selection_pval[which(is.na(D7_7_5_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-7.5V cmh done")

cut_off <- (0.001)

af$D7_7_5_V_selection_pval <- D7_7_5_V_selection_pval

length(which(af$D7_7_5_V_selection_pval < cut_off))

################################################################################################

### (4) cmh for D1-8.1S vs. D7-7.0V

D7_7_0_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_5_V|D7_7_5_S|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_7_0_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D7_7_0_V_selection_pval[which(is.na(D7_7_0_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-7.0V cmh done")

cut_off <- (0.001)

af$D7_7_0_V_selection_pval <- D7_7_0_V_selection_pval

length(which(af$D7_7_0_V_selection_pval < cut_off))

################################################################################################

### (5) cmh for D1-8.1S vs. D7-8.1V

D7_8_1_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_5_V|D7_7_5_S|D7_8_1_S|D7_7_0_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_8_1_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D7_8_1_V_selection_pval[which(is.na(D7_8_1_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-8.1V cmh done")

cut_off <- (0.001)

af$D7_8_1_V_selection_pval <- D7_8_1_V_selection_pval

length(which(af$D7_8_1_V_selection_pval < cut_off))

################################################################################################

### (6) cmh for D1-8.1S vs. D7-8.1S

D7_8_1_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_5_V|D7_7_5_S|D7_8_1_V|D7_7_0_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
                           paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D7_8_1_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D7_8_1_S_selection_pval[which(is.na(D7_8_1_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-8.1S cmh done")

cut_off <- (0.001)

af$D7_8_1_S_selection_pval <- D7_8_1_S_selection_pval

length(which(af$D7_8_1_S_selection_pval < cut_off))

################################################################################################

### (7) cmh for D1-8.1S vs. D1-7.5S

D1_7_5_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D7_7_5_S|D7_7_0_S|D7_7_0_V|D7_7_5_V|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D1_7_5_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D1_7_5_S_selection_pval[which(is.na(D1_7_5_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D1-7.5S cmh done")

cut_off <- (0.001)

af$D1_7_5_S_selection_pval <- D1_7_5_S_selection_pval

length(which(af$D1_7_5_S_selection_pval < cut_off))

################################################################################

### (8) cmh for D1-8.1S vs. D1-7.0S

D1_7_0_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D7_7_0_S|D1_7_5_S|D7_7_5_S|D7_7_0_V|D7_7_5_V|D7_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  D1_7_0_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


D1_7_0_S_selection_pval[which(is.na(D1_7_0_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D1-7.0S cmh done")

cut_off <- (0.001)

af$D1_7_0_S_selection_pval <- D1_7_0_S_selection_pval

length(which(af$D1_7_0_S_selection_pval < cut_off))


################################################################################################
################################################################################################
############################### DAY 7 COMPARISONS TO D7-8.1S ###################################
################################################################################################
################################################################################################


### (9) cmh for D7-8.1S vs. D7-7.5S

DAY7_D7_7_5_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_0_V|D7_7_5_V|D1_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  DAY7_D7_7_5_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


DAY7_D7_7_5_S_selection_pval[which(is.na(DAY7_D7_7_5_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D7-8.1S vs. D7-7.5S cmh done")

cut_off <- (0.001)

af$DAY7_D7_7_5_S_selection_pval <- DAY7_D7_7_5_S_selection_pval

length(which(af$DAY7_D7_7_5_S_selection_pval < cut_off))

################################################################################

### (10) cmh for D7-8.1S vs. D7-7.0S

DAY7_D7_7_0_S_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_5_S|D7_7_0_V|D7_7_5_V|D1_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  DAY7_D7_7_0_S_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


DAY7_D7_7_0_S_selection_pval[which(is.na(DAY7_D7_7_0_S_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D1-8.1S vs. D7-7.0S cmh done")

cut_off <- (0.001)

af$DAY7_D7_7_0_S_selection_pval <- DAY7_D7_7_0_S_selection_pval

length(which(af$DAY7_D7_7_0_S_selection_pval < cut_off))

################################################################################################

### (11) cmh for D7-8.1S vs. D7-7.5V

DAY7_D7_7_5_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_0_V|D7_7_5_S|D1_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  DAY7_D7_7_5_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


DAY7_D7_7_5_V_selection_pval[which(is.na(DAY7_D7_7_5_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D7-8.1S vs. D7-7.5V cmh done")

cut_off <- (0.001)

af$DAY7_D7_7_5_V_selection_pval <- DAY7_D7_7_5_V_selection_pval

length(which(af$DAY7_D7_7_5_V_selection_pval < cut_off))

################################################################################################

### (12) cmh for D7-8.1S vs. D7-7.0V

DAY7_D7_7_0_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  sub_all$pH <- substr(sub_all$ind, 4,6)
  #sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_5_V|D7_7_5_S|D1_8_1_S|D7_8_1_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
                           paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  DAY7_D7_7_0_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}


DAY7_D7_7_0_V_selection_pval[which(is.na(DAY7_D7_7_0_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D7-8.1S vs. D7-7.0V cmh done")

cut_off <- (0.001)

af$DAY7_D7_7_0_V_selection_pval <- DAY7_D7_7_0_V_selection_pval

length(which(af$DAY7_D7_7_0_V_selection_pval < cut_off))

################################################################################################

### (13) cmh for D7-8.1S vs. D7-8.1V

DAY7_D7_8_1_V_selection_pval <- c() # empty vector to save output

for(i in 1:nrow(A1)){
  
  #pull out snp
  sub_A1 <- A1[i,]
  sub_A2 <- A2[i,]
  #transform to data frame
  sub_A1 <- stack(sub_A1)
  sub_A2 <- stack(sub_A2)
  # add ID
  sub_A1$allele <- rep("ac1", nrow(sub_A1))
  sub_A2$allele <- rep("ac2", nrow(sub_A2))
  # add col ID
  colnames(sub_A1) <- c("count", "ind", "allele")
  colnames(sub_A2) <- c("count", "ind", "allele")
  # combine all
  sub_all <- rbind(sub_A1, sub_A2)
  
  # add ids
  #sub_all$pH <- substr(sub_all$ind, 4,6)
  sub_all$condition <- substr(sub_all$ind, 8,8)
  #sub_all$day <- substr(sub_all$ind, 1,2)
  sub_all$group <- substr(sub_all$ind, 1,8)
  
  #only using control D1, removing:
  sub_all <- sub_all[grep("D1_7_0_S|D1_7_5_S|D7_7_0_S|D7_7_5_V|D7_7_5_S|D1_8_1_S|D7_7_0_V", sub_all$ind, invert=TRUE),]
  sub_all$replicate <- ave(paste(sub_all$condition, sub_all$allele, sep=":"),
                           paste(sub_all$condition, sub_all$allele, sep=":"), FUN=seq_along)
  
  # make table for cmh
  Data.xtabs = xtabs(count ~ allele + group + replicate,
                     data=sub_all)
  
  # cmh test using counts for each allele
  test <- mantelhaen.test(Data.xtabs)
  
  # add pvalue to output
  DAY7_D7_8_1_V_selection_pval[i] <- test$p.value
  if (i%%5000 == 0){print(i)}
}

DAY7_D7_8_1_V_selection_pval[which(is.na(DAY7_D7_8_1_V_selection_pval))] <- 1 # bc some are invariant due to removing the unused group

print("D7-8.1S vs. D7-8.1V cmh done")

cut_off <- (0.001)

af$DAY7_D7_8_1_V_selection_pval <- DAY7_D7_8_1_V_selection_pval

length(which(af$DAY7_D7_8_1_V_selection_pval < cut_off))

################################################################################################

head(af)

write.table(file="cmh.out.txt", af, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
write.csv(file="cmh.out.csv", af)


