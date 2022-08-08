####Adapted from: https://github.com/ksamuk/gene_flow_linkage/tree/master/shared_functions
#### sample N SNPs (number of outliers) and extract LD for these SNPs, and return the mean among these SNPs.

calculate.null.metrics <- function(ld, ol){
  
  num.outliers <- nrow(ol)
  
  snp.list <- c(as.character(ld$SNP_A), as.character(ld$SNP_B))
  snp.list <- as.data.frame(unique(snp.list))
  colnames(snp.list) <- "snp"
  
  
  snp.sample <- snp.list %>% 
    sample_n(num.outliers) %>% 
    pull(snp)
  
  ld.sample <- ld %>% 
    filter(SNP_A %in% snp.sample & SNP_B %in% snp.sample)
  
  
  #ld
  mean <- mean(ld.sample$R2)
  
  return(mean)
}