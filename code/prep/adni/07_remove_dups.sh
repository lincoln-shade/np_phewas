
plink --bfile data/mega/all_cohorts_notqced --remove data/mega/dup_ids.txt --make-bed --out data/mega/all_cohorts_np

cat data/mega/all_cohorts_np.fam | awk '$1=="ADNI"{print $1, $2}' > data/adni/ids_qced_not_dups.txt

# plink --bfile data/mega/all_cohorts_np --keep data/adni/ids_qced_not_dups.txt --make-bed --out data/adni/adni_np
