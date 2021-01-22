#####################################################
#                                                   #
#   Script for identifying diploids and triploids,  #
#   among GBS data.                                 #
#                                                   # 
#   Author: Roos Goessen                            #
#   Date: 2020-12-14                                #
#                                                   #                     
#####################################################

#Note: make sure to only use bi-allelic snps
#Note: you can modify the Ratios if you wish to look for other ploidy levels

# load packages
library(tidyverse); library(vcfR)
library(forcats); library(ggplot2)

# decide depth cutoff 
# the higher the cut-off, the more clear will be the prediction, however if the cut-off is too high the will not be enough SNPs. 
depth=16

#######################
#                     #
#   Parse vfc file    #
#                     # 
#######################

# read files
input_vcf <- ""
yourvcf <- read.vcfR(input_vcf)

# make sure the order of the VCF aspects (GT,DP etc.) is the same for your own VCF file
# gather data from vcffile
VCF_file <- data.frame(yourvcf@gt) %>% gather(sample, gt, -FORMAT) %>% separate(gt, into = c("GT","DP","AD","GQ","GL"), sep = ":") %>% separate(AD, into = c("count_REF", "count_ALT"), sep = ",")

# take depth reference and divide it by the total depth (ref + alt). 
VCF_file$HET <- as.numeric(VCF_file$count_REF)/(as.numeric(VCF_file$count_REF)+as.numeric(VCF_file$count_ALT))

# filter out snps with frequencies below 0.1 and over 0.9 to only look at heterozygous SNPs 
VCF_file <- VCF_file %>% filter(HET > 0.1) %>% filter(HET < 0.9)

# only SNPs with coverage above certain cut-off depth, to avoid background noise. Also to define triploids correctly
VCF_file <- VCF_file %>% filter((as.numeric(VCF_file$count_REF)+as.numeric(VCF_file$count_ALT))>depth)

#########################
#                       #
#   Calculate Ratios    #
#                       # 
#########################

# Define Range (+-6)
# Range_A (33.3):  
# Range_B (50):
# Range_C (66.7):

VCF_file.mutated <- VCF_file %>% mutate(CLASS = case_when(
  HET > 0.273 & HET < 0.393 ~ "RANGE_A",
  HET > 0.440 & HET < 0.560 ~ "RANGE_B",
  HET > 0.607 & HET < 0.727 ~ "RANGE_C",
  TRUE ~ "Other"))

# count number of SNPs per range, per sample
Freq_per_sample <- as.data.frame(unclass(table(VCF_file.mutated$sample,VCF_file.mutated$CLASS)))

# total number of snps counted
Freq_per_sample$TOTAL <- Freq_per_sample$RANGE_A+Freq_per_sample$RANGE_C+Freq_per_sample$RANGE_B+Freq_per_sample$Other

# Calculation of ratio A+C/B
Freq_per_sample$RATIO <- (Freq_per_sample$RANGE_A+Freq_per_sample$RANGE_C)/Freq_per_sample$RANGE_B

# Calculation of ratio C/B
Freq_per_sample$RATIO_high <- Freq_per_sample$RANGE_C/Freq_per_sample$RANGE_B

# add names sample
Freq_per_sample$Sample <- rownames(Freq_per_sample)


#########################
#                       #
#   Visualize Ratios    #
#                       # 
#########################

# plot for count per ratio value
a <- ggplot(Freq_per_sample, aes(x = RATIO)) + geom_dotplot(binwidth = 1/30) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
a

# plot total number snps vs ratio
plot(Freq_per_sample$RATIO~Freq_per_sample$TOTAL,xlab="Number SNPs (het 0.1-0.9),depth>16", 
     ylab= "Ratio")
abline(h = 1.5) # you can choose ratio cut-off yourself

# Choose ratio based on plots
CHOSEN_RATIO <- 1.5

#########################
#                       #
#       Save data       #
#                       # 
#########################

# add new ploidy based on script 
Freq_per_sample$ploidy_script <- ifelse(Freq_per_sample$RATIO<CHOSEN_RATIO, "diploid", "triploid")

# save data
write.csv(Freq_per_sample, paste("Ratios_ploidy_",depth,".csv", sep="",collapse=NULL),row.names=FALSE)

