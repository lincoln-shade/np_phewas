#!/bin/bash


seq 1 22 | sudo xargs -P 22 -n 1 bash -c 'bcftools merge data/imputed/ADC_NHW/ADC*/chr$0.dose.vcf.gz --force-samples -Oz -o /data_global/ADC_imputed/chr$0_all_adcs.dose.vcf.gz'
