library(data.table)
library(ggplot2)
library(stringi)
dt <- fread("data/rosmap_data.csv")
dt[, caa_4gp := ordered(caa_4gp, 
                        labels = c("None", "Mild", "Moderate", "Severe"))]
dt[, apoe := factor(apoe, 
                    levels = c("apoe33", "apoe34", "apoe44", 
                               "apoe23", "apoe22"))]
cpgs <- colnames(dt)[grep("cg", colnames(dt))]
# regress cpg sites on covariates
get_residuals <- function(cpg, 
                          covars=c("conversion.efficiency", "msex", "apoe",
                                   "age_death", "pmi", "study")) {
    f <- as.formula(paste0(cpg, " ~ ", paste(covars, collapse = "+")))
    m <- dt[, lm(f, .SD)]
    res <- dt[, residuals(lm(f, .SD))]
    res <- data.table(as.integer(names(res)), res)
    setnames(res, colnames(res), c("row_num", "residual"))
}

for (i in 1:length(cpgs)) {
    res <- get_residuals(cpgs[i])
    dt[res[["row_num"]], (paste0(cpgs[i], "_res")) := scale(res[["residual"]])]
    rm(res)
}

plot_data <- melt(dt[!is.na(caa_4gp), 
                     .(caa_4gp, cg04401876_res, cg09555818_res,
                       cg10169327_res, cg13119609_res)],
                  id.vars = "caa_4gp",
                  variable.factor = FALSE)

plot_data[, variable := stri_replace_first_fixed(variable, "_res", "")]
p1 <- plot_data[value > -4, ggplot(.SD, aes(caa_4gp, value, fill=variable, alpha=0.8))] +
    facet_wrap("variable", nrow = 1) +
    geom_boxplot(aes(alpha = 1), show.legend = FALSE, outlier.alpha = 0) +
    geom_jitter(aes(alpha = 0.8), width = 0.15, show.legend = FALSE) +
    labs(y = "Methylation at CpG site",
         x = "CAA Pathology Severity",
         fill = "CAA Severity",
         title = "Methylation at CpG sites vs. cerebral amyloid angiopathy pathology") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11),
          strip.text.x = element_text(size = 14))
p1
