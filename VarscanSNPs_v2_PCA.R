######################################### snp_analysis.R #########################################
## pca, etc for full data set, f1, f2, f3.


install.packages("ggplot2")
install.packages("reshape")
install.packages("data.table")
install.packages("gridExtra")
install.packages("scales")


library(stringr) # yes
library(ggplot2)
library(reshape)
library(data.table)
library(gridExtra)
library(scales)

af <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_allele_freqs_v2.txt", header=TRUE)
dat3 <- read.table("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/filtered_variants_v2.txt", header=TRUE)

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


freqs <- t(af[,2:ncol(af)])
colnames(freqs) <- c(paste(dat3$Chrom, dat3$Position, sep=":"))


# apply PCA - scale. = TRUE is highly
# advisable, but default is FALSE.

pcaResult <- prcomp(freqs, scale=TRUE)
#round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)
#summary(pcaResult)






#########################################################################################################

####
##
## plot pca: pH vs. Day
##
####

# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(pH=substr(pops, 4,6), condition=substr(pops, 8,8),
                    day=substr(pops, 1,2),
                    PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

head(dat.p)

a <- ggplot(dat.p, aes(PC1, PC2, fill=pH, shape=day)) +
  geom_point(size=4, color="black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('brown3', 'green', 'blue3'))+
  scale_fill_manual(values=c('brown3', 'green', 'blue3'))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23,23), 
                                             size=c(5,5,5), fill=c('brown3', 'green', 'blue3'))))+
  ggtitle("pH vs. Day")

a

png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_pHvsDay.png", res=300, height=5, width=8, units="in")

a

dev.off()

#########################################################################################################

####
##
## plot pca: pH vs. Condition
##
####

b <- ggplot(dat.p, aes(PC1, PC2, fill=pH, shape=condition)) +
  geom_point(size=4, color="black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('brown3', 'green', 'blue3'))+
  scale_fill_manual(values=c('brown3', 'green', 'blue3'))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23,23),
                                             size=c(5,5,5), fill=c('brown3', 'green', 'blue3'))))+
  ggtitle("pH vs. Condition")

b


png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_pHvsCondition.png", res=300, height=5, width=8, units="in")

b

dev.off()


#########################################################################################################

####
##
## plot pca: Day vs. Condition
##
####

c <- ggplot(dat.p, aes(PC1, PC2, fill=day, shape=condition)) +
  geom_point(size=4, color="black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('blue3', 'brown3'))+
  scale_fill_manual(values=c('blue3', 'brown3'))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23),
                                             size=c(5,5), fill=c('blue3', 'brown3'))))+
  ggtitle("Day vs. Condition")
  

c


png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_DayvsCondition.png", res=300, height=5, width=8, units="in")

c

dev.off()


#########################################################################################################

####
##
## plot pca: Day vs. pH vs. Condition
##
####

d <- ggplot(dat.p, aes(PC1, PC2, fill=pH, shape=condition, size=day)) +
  geom_point(color="black", alpha=0.75) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('brown3', 'green', 'blue3'))+
  scale_fill_manual(values=c('brown3', 'green', 'blue3'))+
  scale_size_manual(values=c(6,4))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23,23), 
                                             size=c(5,5,5), fill=c('brown3', 'green', 'blue3'))))+
  ggtitle("pH vs. Condition vs. Day")

d

png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_pHvsConditionvsDay_alpha.png", res=300, height=5, width=8, units="in")

d

dev.off()


#########################################################################################################

####
##
## plot pca: Day vs. pH vs. Condition (only 8.1S for D1)
##
####

freqs2 <- freqs[11:44,]
rownames(freqs2)
nrow(freqs2)
pcaResult2 <- prcomp(freqs2, scale=TRUE)

pops <- c(
  "D1_8_1_S_01", "D1_8_1_S_03",
  "D1_8_1_S_05", "D1_8_1_S_08", "D7_7_0_S_25", "D7_7_0_S_26",
  "D7_7_0_S_27", "D7_7_0_S_28", "D7_7_0_S_29", "D7_7_0_S_30",
  "D7_7_0_V_02", "D7_7_0_V_03", "D7_7_0_V_04", "D7_7_0_V_05",
  "D7_7_0_V_06", "D7_7_5_S_31", "D7_7_5_S_32", "D7_7_5_S_33",
  "D7_7_5_S_36", "D7_7_5_V_07", "D7_7_5_V_08", "D7_7_5_V_09",
  "D7_7_5_V_10", "D7_7_5_V_12", "D7_8_1_S_13", "D7_8_1_S_14",
  "D7_8_1_S_15", "D7_8_1_S_18", "D7_8_1_V_19", "D7_8_1_V_20",
  "D7_8_1_V_21", "D7_8_1_V_22", "D7_8_1_V_23", "D7_8_1_V_24")


percentVar <- round(pcaResult2$sdev^2/sum(pcaResult2$sdev^2)*100, digits=2)

dat.p <- data.frame(pH=substr(pops, 4,6), condition=substr(pops, 8,8),
                    day=substr(pops, 1,2),
                    PC1 = pcaResult2$x[,1],  PC2= pcaResult2$x[,2])



e <- ggplot(dat.p, aes(PC1, PC2, fill=pH, shape=condition, size=day)) +
  geom_point(color="black", alpha=0.75) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('brown3', 'green', 'blue3'))+
  scale_fill_manual(values=c('brown3', 'green', 'blue3'))+
  scale_size_manual(values=c(6,4))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23,23), 
                                             size=c(5,5,5), fill=c('brown3', 'green', 'blue3'))))+
  ggtitle("pH vs. Condition vs. Day")

e

png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_pHvsConditionvsDay_alpha_D1controlonly.png", res=300, height=5, width=8, units="in")

e

dev.off()


#########################################################################################################

####
##
## plot pca: pH vs. Condition (no day)
##
####


freqs3 <- freqs[15:44,]
rownames(freqs3)
nrow(freqs3)
pcaResult3 <- prcomp(freqs3, scale=TRUE)

pops <- c(
  "D7_7_0_S_25", "D7_7_0_S_26",
  "D7_7_0_S_27", "D7_7_0_S_28", "D7_7_0_S_29", "D7_7_0_S_30",
  "D7_7_0_V_02", "D7_7_0_V_03", "D7_7_0_V_04", "D7_7_0_V_05",
  "D7_7_0_V_06", "D7_7_5_S_31", "D7_7_5_S_32", "D7_7_5_S_33",
  "D7_7_5_S_36", "D7_7_5_V_07", "D7_7_5_V_08", "D7_7_5_V_09",
  "D7_7_5_V_10", "D7_7_5_V_12", "D7_8_1_S_13", "D7_8_1_S_14",
  "D7_8_1_S_15", "D7_8_1_S_18", "D7_8_1_V_19", "D7_8_1_V_20",
  "D7_8_1_V_21", "D7_8_1_V_22", "D7_8_1_V_23", "D7_8_1_V_24")


percentVar <- round(pcaResult3$sdev^2/sum(pcaResult3$sdev^2)*100, digits=2)

dat.p <- data.frame(pH=substr(pops, 4,6), condition=substr(pops, 8,8),
                    day=substr(pops, 1,2),
                    PC1 = pcaResult3$x[,1],  PC2= pcaResult3$x[,2])





f <- ggplot(dat.p, aes(PC1, PC2, fill=pH, shape=condition)) +
  geom_point(size=4, color="black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  #ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c(21,22))+
  #scale_color_manual(values=c('brown3', 'green', 'blue3'))+
  scale_fill_manual(values=c('brown3', 'green', 'blue3'))+
  #theme(legend.position = c(0.88,0.17))+
  theme(legend.text=element_text(size=10),legend.title=element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(fill=guide_legend(override.aes=list(shape=c(23,23,23),
                                             size=c(5,5,5), fill=c('brown3', 'green', 'blue3'))))+
  ggtitle("pH vs. Condition: Day 7 only")

f


png("/Users/aprilgarrett/Desktop/MastersThesis/GenomicAnalyses/figures/pca_pHvsCondition_D7only.png", res=300, height=5, width=8, units="in")

f

dev.off()


