
plink --bfile tmp/adni_2 --keep data/adni/ids_qced.txt --make-bed --out tmp/adni_3

python code/prep/merge_plink_filesets.py -b data/act/act_rosmap_adc -b1 tmp/adni_3 -o data/mega/all_cohorts_notqced

king -b data/mega/all_cohorts_notqced.bed --duplicate --prefix data/mega/duplicates_bt_cohorts