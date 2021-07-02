#=============================
# calculate power for PheWAS
#=============================

library(genpwr)

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=c(800), Case.Rate=c(0.3), k=NULL,
                  MAF=seq(0.05, 0.5, 0.05), OR=c(4.5),Alpha=5e-8,
                  True.Model=c("Additive"), 
                  Test.Model=c("Additive"))
pw
