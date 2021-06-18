#!/bin/bash


CHR=$1

cat /data_global/UKY_ADC/TOPMedImputation/preQC/chr$CHR.alladcs.bim | xargs -n 6 bash -c 'grep -E -w -m1 "$3[ \t].*$4[ \t].*$5|$3[ \t].*$5[ \t].*$4" /data_global/code/rs/hg38/ucsc.hg38.snp151.chr'$CHR'.linker | (printf $1; printf " ";cut -f 3 -d " " | tr -d "\n";printf "\n") | awk '\''$2!="" {print}'\'''  > /data_global/UKY_ADC/TOPMedImputation/linker/chr$CHR.rs.linker





