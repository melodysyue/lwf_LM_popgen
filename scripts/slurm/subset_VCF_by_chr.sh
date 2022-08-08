#module load samtools
#subset vcf by chrs
while IFS= read -r line; do
  tabix -h lwf.nomixed.snps.N829.sorted.chr.vcf.gz $line > ../lostruct/$line.vcf;
done < lwf.chr.names
#-h : keep header in the vcf files;
