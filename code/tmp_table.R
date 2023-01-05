
library(data.table)
library(magrittr)
library(flextable)
library(table1)

pord <- fread("data/mega/mega_np_ord.pheno", na.strings = '-1')
pbin <- fread("data/mega/mega_np.pheno", na.strings = '-1')
pheno <- merge(pord, pbin, c("FID", "IID"))
covar <- fread("data/mega/mega_np.covar")
covar_apoe <- fread("data/mega/mega_np_apoe.covar", stringsAsFactors = TRUE)
covar <- merge(covar, 
               covar_apoe[, .(FID, IID, apoe)], 
               c("FID", "IID"), 
               all.x = TRUE)
covar <- covar[, .(FID, IID, msex, age_death, apoe)]
covar[, Sex := factor(msex, labels = c("Female", "Male"))]

dt <- merge(pheno, covar, by = c("FID", "IID"))
setnames(dt, "age_death", "Age of Death")
dt[, Study := factor(FID, labels = c("NACC", "ROSMAP", "ACT", "ADNI"))]
dt[, CAA := ordered(caa_ord, 
                      labels = c("None", "Mild", "Moderate", "Severe"))]

# Display apoe variable as number of APOE e4 alleles for table
dt[apoe %in% c("apoe22", "apoe23", "apoe33"), 
   `APOE e4 alleles` := 0]
dt[apoe %in% c("apoe24", "apoe34"), 
   `APOE e4 alleles` := 1]
dt[apoe %in% c("apoe44"), 
   `APOE e4 alleles` := 2]
dt[, `APOE e4 alleles` := ordered(`APOE e4 alleles`)]

t1dt <- 
  as.data.frame(
    table1(~ Sex + `Age of Death` + `APOE e4 alleles` + CAA | Study,
           data = dt
    )
  ) 

t1 <- flextable(t1dt) %>% 
  autofit() %>% 
  padding(padding = 2, part = "all") %>%
  add_header_lines("Table 1: Participant Demographics")

t1
