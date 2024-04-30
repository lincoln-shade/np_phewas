# cat output/gwas/mega/polmm/atheroscler_polmm_results.txt | \
#   ~/software/locuszoom/bin/locuszoom \
#     --metal - \
#     --markercol SNPID \
#     --pvalcol pval.spa \
#     --chr 13 \
#     --refsnp rs2000660 \
#     --flank 500000 \
#     --source 1000G_March2002 \
#     --build hg19 \
#     --pop EUR \
#     --delim space \
#     --no-date \
#     ymax=10 \
#     legend="none" \
#     geneFontSize=1.18 \
#     axisTextSize=1.4 \
#     refsnpTextSize=1.2 \
#     showRecomb=FALSE \
#     rfrows=2	\
#     --prefix doc/ 

cat "output/gwas/metal/results/caa_apoe1.csv" | \
  ~/software/locuszoom/bin/locuszoom \
    --metal - \
    --markercol MarkerName \
    --pvalcol P.value \
    --chr 19 \
    --refsnp rs4803778 \
    --flank 500000 \
    --source 1000G_March2012 \
    --build hg19 \
    --pop EUR \
    --delim space \
    --no-date \
    ymax=13 \
    legend="none" \
    geneFontSize=1 \
    axisTextSize=1 \
    refsnpTextSize=1.2 \
    showRecomb=FALSE \
    signifLine="7.3" \
    --prefix doc/