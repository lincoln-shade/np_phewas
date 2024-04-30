library(data.table)
library(ordinal)
coloc <- fread("data/coloc_results.txt")
coloc <- coloc[!duplicated(coloc)]
coloc <- coloc[QTL_Source == "ROSMAP"]
coloc80 <- coloc[QTL_Type == "mQTL" & PPH4 >= 0.8]
npe <- coloc80[, names(table(GWAS_Phenotype))]
load("data/cpg.RData")
np_cpgs <- coloc80[, QTL_Phenotype]
cpg_cols <- c("IID", "specimenID", "batch", "conversion.efficiency", np_cpgs)
cpg <- cpg[, ..cpg_cols]

np <- fread("data/rosmap_data.csv")
np[, IID := as.character(projid)]

dt <- merge(np, cpg, "IID")



run_regression <- function(npe, cpg) {
    f <- as.formula(paste0(
        "as.ordered(", npe, ") ~ ", cpg, " + scale(age_death) + msex + factor(study) + scale(pmi) + factor(apoe) + batch + scale(conversion.efficiency.x)"))
    m <- clm(f, data=dt)
    return(m)
}

setnames(dt, c("arteriol_scler", "cvda_4gp2", "braaksc", "plaq_d", 
               "ci_num2_gct", "hspath_typ", "tdp_st4", "dlbdx", "ci_num2_mct"),
         npe)

npe1 <- "arteriol"
cpg1 <- "cg00438541"

m1 <- run_regression(npe1, cpg1)

coloc80_ord <- coloc80[GWAS_Phenotype != "diffuse_abeta"]
mlist <- vector("list", nrow(coloc80_ord))

for (i in 1:length(mlist)) {
    npei <- coloc80_ord$GWAS_Phenotype[i]
    cpgi <- coloc80_ord$QTL_Phenotype[i]
    mlist[[i]] <- run_regression(npei, cpgi)
}


