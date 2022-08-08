####Adapted from: https://github.com/ksamuk/gene_flow_linkage/tree/master/shared_functions
#### sample N SNPs (number of outliers) and extract LD for these SNPs, and return the mean among these SNPs.

calculate.emp.metrics <- function(ld, ol){
  
  ld.ol <- ld %>% 
    filter(SNP_A %in% ol$snp & SNP_B %in% ol$snp)
  
  #ld
  mean <- mean(ld.ol$R2)
  
  return(mean)
}