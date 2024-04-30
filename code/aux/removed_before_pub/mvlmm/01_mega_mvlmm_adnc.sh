cat output/gwas/mega/mega_adnc_pca.PC1.glm.linear | \
  awk '$14 < 5e-4 {print $3}' \
  > tmp/gemma_snps.tmp

cat data/mega/mega_np_ord.pheno | \
  awk '$3!=-1 && $4!=-1 && $5!=-1 && $10!=-1 && NR!=1 {
       print$1,$2,0,0,0,$3,$4,$5,$10
  }' \
  > tmp/mega_mvlmm_adnc.txt

plink2 \
  --pfile data/mega/mega_np \
  --keep tmp/mega_mvlmm_adnc.txt  \
  --make-bed \
  --out tmp/mega_mvlmm_adnc


ldak5.2.linux \
  --bfile tmp/mega_mvlmm_adnc \
  --calc-kins-direct tmp/mega_mvlmm_adnc \
  --kinship-gz YES  \
  --ignore-weights YES  \
  --power -0.25 \
  --hwe-stand NO  \
  --max-threads 32

Rscript --vanilla code/mvlmm/01a_format_grm_covar.txt

~/software/GEMMA-0.98.5/bin/gemma \
  -bfile tmp/mega_mvlmm_adnc  \
  -k tmp/mega_mvlmm_adnc_grm.txt  \
  -km 2 \
  -n 1 2 3 4  \
  -lmm  \
  -c tmp/mega_mvlmm_adnc.covar  \
  -outdir ./tmp/ \
  -o mega_mvlmm_adnc \
  -snps tmp/gemma_snps.tmp

