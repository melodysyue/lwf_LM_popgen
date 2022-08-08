#module load samtools

while IFS= read -r CHR; do

# print progress
echo ${CHR}

# set input and output
OUTFILE=/scratch-lustre/wlarson/yue/LWF/lostruct/geno_chrs/${CHR}.geno
VCF=/scratch-lustre/wlarson/yue/LWF/lostruct/vcf_chrs/${CHR}.vcf
HEAD=/scratch-lustre/wlarson/yue/LWF/lostruct/geno_chrs/${CHR}.header

# initialise file
echo -ne "chromosome\tpos\t" > $HEAD

# echo sample names and also deal with running line
# head line: chromosome pos each sample name
bcftools query -l $VCF | tr "\n" "\t" >> $HEAD
awk '{ if (NR%2) {printf $0} else {print $0} }' FS='\t' $HEAD > $OUTFILE

# output IUPAC genotypes and replace missing genotype ./. or . with N
#GT: genotype, encoded as alleles valued separated by eitehr / or |;
#O for ref; 1 for alt;
#for diploid: 0/1 or 1|0 etc;
#for haploid: only one allele should be given;
#missing genotype: ./. or .

bcftools query -e 'STRLEN(REF)>1' -f '%CHROM\t%POS[\t%IUPACGT]\n' $VCF | sed 's/\.\/\./N/g' | sed 's/\./N/g' | sed 's/N1/\.1/g' >> $OUTFILE
gzip $OUTFILE

done < lwf.chr.names 
