library(tidyverse)
library(bigsnpr)
library(readr)


###background

###we use bigsnpr for two purposes:
###(1) to thin dataset based on LD and then run PCAadapt for outlier detection using the best practice;
###(2) detect long-range LD. 

###Output:
##list of related and non-related individuals: second-degrees or above; 
##list of prunning SNPs after removing outlier SNPs; 
##list of long-range LD if any. 


###Input: bed/bim/fam files generated from PLINK
###It doesn't handle missing values, use snp_fastInput() or snp_fastInputSimple() to impute missing values of genotyped variant. 

#download plink executable. only need to do this once.
#plink <- download_plink2("./plink/")

rm(list=ls())
plink <- "./scripts/1_neutral_vs_outlier/plink2.exe"
bedfile <- './data/plink/lwf.west.pop.bed'


popmap <- read.table("./data/plink/lwf.west.pop.fam", header=F)
popmap <- popmap %>% 
  select(sample=V2, pop=V1)

head(popmap)

###check the data
(obj.bed <- bed(bedfile))


###Get the input ready
# Read from bed/bim/fam, it will create new files.
snp_readBed(bedfile) #only need to do this once
#this creates two new output: X.rds and X.bk
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("./data/plink/lwf.west.pop.rds") 
#load the "bigSNP" object


# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")

###Get alases for useful slots
G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
#CHR <- parse_number(CHR) #bigsnpr only takes chrs as integers
POS <- obj.bigSNP$map$physical.pos
# Check and make sure there is no missing data
big_counts(G, ind.col = 1:10)

G2 <- snp_fastImputeSimple(G) #imput missing data
big_counts(G2, ind.col = 1:10)

###Automatic processure to clump and remove long-range LD
svd <- snp_autoSVD(G2, CHR, POS, 
                   min.mac=3,
                   max.iter=10,
                   ncores=nb_cores(),
                   roll.size=0) #set roll.size=0 otherwise it will say it is too big. so annoying!

#if it says "XXX" already exists, run it again. 

subset=attr(svd, "subset")

thinned <- obj.bigSNP$map %>% 
  rownames_to_column("index") %>% 
  filter(index %in% subset) %>% 
  pull(marker.ID)

write.table(thinned, "./data/neutral_outlier/lwf.pcadapt.west.thinned.snps", quote = F, row.names = F)
length(thinned)





