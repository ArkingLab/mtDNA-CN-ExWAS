dx cd project-G7KB5zQJz55qG0158ZXb2J5p:/vkp/TFAM
dx download TFAM_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt
dx download TFAM_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt
dx download TFAM_mtDNA_CN_n180k_MAC_gteq6_variants.FDR10pct.txt

less TFAM_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt | cut -f1 | tail -n +2 > TFAM_mtDNA_CN_n180k_MAC_lt6_variants.txt
plink --bfile ~/ref/UKB\ -mtDNA_2022/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ interim\ 450k\ release/ukb23149_c10_b0_v1 --extract TFAM_mtDNA_CN_n180k_MAC_lt6_variants.txt --recode AD --out TFAM_mac_lt6_variants

less TFAM_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt | cut -f1 | tail -n +2 > TFAM_mtDNA_CN_n180k_MAC_gteq6_variants.txt
plink --bfile ~/ref/UKB\ -mtDNA_2022/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ interim\ 450k\ release/ukb23149_c10_b0_v1 --extract TFAM_mtDNA_CN_n180k_MAC_gteq6_variants.txt --recode vkp_snapshot_jan27.snapshotAD --out TFAM_mac_gteq6_variants

rm *.nosex *.log
