#!/bin/sh
# gets unique positions from VEP output file

dir="/dcs04/arking/data/active/projects/mito_rare-variant/annotations/qualifying_variants/protein_coding/"
#dir="/dcs04/arking/data/active/projects/mito_rare-variant/annotations/qualifying_variants/protein_coding_phred18/"
out_dir="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/qual.var/"
echo -e "directory: ${dir}"
#out_prefix="variants_in_cds"
#out_prefix="variants_in_cds_LoF2"
#out_prefix="variants_in_cds_CADDphred18_LoF"
#out_prefix="variants_in_cds_CADDphred_18"
#out_prefix="variants_in_cds_CADDphred_18_impact_high"
out_prefix="variants_in_cds_synonymous"
mkdir -p ${out_dir}/${out_prefix}/
echo -e "out directory: ${out_dir}/${out_prefix}/"

#for i in {1..22}
for i in X
do

echo -e "chr: ${i}"
# note: this code "| awk '$11 ~ /\// {print}'" removes any synonymous variants
#cat ${dir}/ukb23156_c${i}_v1_cadd_coding.txt | grep -v '#' | awk '$11 ~ /\// {print}' | cut -f2 | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos
#cat ${dir}/ukb23156_c${i}_v1_coding_phred18.txt | grep -v '#' | awk '$11 ~ /\// {print}' | cut -f2 | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos
#cat ${dir}/ukb23156_c${i}_v1_coding_phred18.txt | grep -v '#' | awk '$11 ~ /\// {print}' | grep 'IMPACT=HIGH' | cut -f2 | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos
#cat ${dir}/ukb23156_c${i}_v1_cadd_coding.txt | grep -v '#' | awk '$7 !~ /synonymous/ && $7 !~ /missense/ && $11 ~ /\// {print$2}' | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos
#cat ${dir}/ukb23156_c${i}_v1_coding_phred18.txt | grep -v '#' | awk '$7 !~ /synonymous/ && $7 !~ /missense/ && $11 ~ /\// {print$2}' | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos
cat ${dir}/ukb23156_c${i}_v1_cadd_coding.txt | grep -v '#' | awk '$7 ~ /synonymous/ {print$2}' | tr ':' '\t' | tr '-' '\t' | cut -f2 | sort -n | uniq > ${out_dir}/${out_prefix}/${out_prefix}.chr${i}.uniq_pos

done




