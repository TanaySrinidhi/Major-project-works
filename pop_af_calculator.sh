Populations=(
     "barbados_afrocarribean"
     "colombian"
     "mexican_la"
     "peruvian"
     "puerto_rican"
     "southwest_us_afroamerican"
     "utah_us_european"
     "bangladesh_bengali"
     "british"
     "chinese_dai"
     "chinese_han_beijing"
     "chinese_han_south"
     "esan_nigeria"
     "finland_finnish"
     "gambia_mandinka"
     "gujarati_houston"
     "italy_toscani"
     "japanese_tokyo"
     "kenya_luhya"
     "nigerian_yoruba"
     "punjabi_lahore"
     "sierra_leone_mende"
     "spain_iberian"
     "srilanka_tamil_uk"
     "telugu_uk"
     "vietnam_kinh"
)

#For-loop to calculate allele frequency for NTD-linked SNPs in specific populations using superpopulation data
for population in "${Populations[@]}"; do
    bcftools view "$population"_lct_active_snps_only.vcf -Oz -o "$population"_lct_active_snps_only.vcf.gz
    # Index
    bcftools index "$population"_lct_active_snps_only.vcf.gz

    # Calculate population-specific AF
    bcftools +fill-tags "$population"_lct_active_snps_only.vcf.gz -- -t AF | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    > "$population"_lct_allele_frequencies.txt
done  
