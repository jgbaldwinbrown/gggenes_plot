#!/bin/bash
set -e

cat _breed_BlackHomer_time_36_bit_Bitted_replicate_R1_breed_Feral_time_36_bit_Bitted_replicate_R1_Color_Low_Mid__subsetted_input.vcf | \
awkf '$1!~/^#/ {print $1, $2-1, $2, "snp"}' \
> snps.bed

cat _breed_WhiteHomer_time_36_bit_Unbitted_replicate_All_breed_WhiteHomer_time_36_bit_Bitted_replicate_All_Color_Low_Mid___multiplot_selec_plfmt_bedified.bed | \
awkf '{print $1, $2, $3, $4, "selec"}' \
> data.bed

gunzip -c louse_annotation_0.1.1.gff.gz | \
mawk -F "\t" -v OFS="\t" '$1=="chr1"' | \
gzip -c > louse_annotation_0.1.1_chr1.gff.gz
