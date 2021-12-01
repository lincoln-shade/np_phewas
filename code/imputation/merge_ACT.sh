#!/bin/bash


seq 1 22 | sudo xargs -P 22 -n 1 bash -c 'bcftools merge --force-samples -Oz -o /data_global/ACT_imputed/chr$0.dose.vcf.gz data/imputed/ACT_NHW/ACT*/chr$0.dose.vcf.gz'
