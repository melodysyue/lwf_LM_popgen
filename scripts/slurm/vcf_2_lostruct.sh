#module load samtools

while IFS= read -r CHR; do

# print progress
echo ${CHR}

plink --vcf ./vcf_chrs/${CHR}.vcf --recode12 --chr-set 40 no-xy --out ./lostruct_input/${CHR} 
plink --file ./lostruct_input/${CHR} --recodeA --chr-set 40 no-xy --out ./lostruct_input/${CHR}
cat ./lostruct_input/${CHR}.raw | sed 1d | cut -d' ' -f7- > ./lostruct_input/${CHR}.geno
cut -f1,4 ./lostruct_input/${CHR}.map > ./lostruct_input/${CHR}.loci
done < lwf.chr.names 
