#!/bin/sh
# Run `ukb_plinkBED2gds.R` on each of chr1-22 of UKB 450k Exome data
for c in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
  if [ ! -f "ukb23149_c${c}_b0_v1.bed" ]; then
      dx cd /Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ interim\ 450k\ release/
      dx download ukb23149_c${c}_b0_v1.*
  fi

  ~/bin/R-4.1.2/bin/Rscript ukb_plinkBED2gds.R --bed_dir ./ --output_dir ./ --num.cores 10 ${c} 2>&1 | tee UKBexomeQQFE_chr${c}.gds.log
  dx cd /Data/gds/
  dx upload UKBexomeQQFE_chr${c}.gds*
  rm ukb23149_c${c}_b0_v1.*
done


# create symbolic links to exome plink files (combined common/rare)
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  echo "chr: ${c}"
  ln -s /home/dnanexus/ref2/UKB\ -mtDNA_2022/Data/merged/combined_chr${c}.bed combined_chr${c}.bed
  ln -s /home/dnanexus/ref2/UKB\ -mtDNA_2022/Data/merged/combined_chr${c}.bim combined_chr${c}.bim
  ln -s /home/dnanexus/ref2/UKB\ -mtDNA_2022/Data/merged/combined_chr${c}.fam combined_chr${c}.fam
done