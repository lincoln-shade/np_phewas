cat output/gwas/mega/polmm/atheroscler_polmm_results.txt | \
  ~/software/locuszoom/bin/locuszoom \
    --metal - \
    --markercol SNPID \
    --pvalcol pval.spa \
    --chr 13 \
    --refsnp rs2000660 \
    --flank 200000 \
    --source 1000G_Nov2014 \
    --build hg38 \
    --pop EUR \
    --delim space \
    --no-date \
    ymax=10 \
    legend="none" \
    geneFontSize=1.18 \
    axisTextSize=1.4 \
    refsnpTextSize=1.2 \
    showRecomb=FALSE \
    rfrows=2	\
    --prefix doc/ \
    --rundir ~/software/locuszoom