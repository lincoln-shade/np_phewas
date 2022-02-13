
python code/prep/merge_plink_filesets.py -b data/rosmap/rosmap_adc -b1 data/act/act_np -o data/act/act_rosmap_adc

king -b data/act/act_rosmap_adc.bed --related --degree 2 --prefix data/act/act_rosmap_adc