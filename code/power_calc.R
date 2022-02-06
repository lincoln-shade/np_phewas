#=============================
# calculate power for PheWAS
#=============================
#=============================
# calculate power for PheWAS
#=============================

library(genpwr)
library(ggplot2)
library(data.table)
library(scales)

# NACC + ROSMAP + ACT + ADNI approx. N
n_case <- 2479
n_control <- 3878
n <- n_case + n_control
case_rate <- n_case / n
n_tests <- 1000000
alpha <- 0.05 / n_tests
maf <- seq(0.01, 0.5, 0.01)
or <- seq(1.3, 1.7, 0.1)
pw <- genpwr.calc(calc = "power", 
                  model = "logistic", 
                  ge.interaction = NULL,
                  N=n, 
                  Case.Rate=case_rate, 
                  k=NULL,
                  MAF= maf, 
                  OR= or,
                  Alpha= alpha,
                  True.Model=c("Additive"), 
                  Test.Model=c("Additive"))
pw <- as.data.table(pw)
colnames(pw)[9] <- "Power"

pw_plot <- pw[, ggplot(.SD, aes(MAF, Power, color=factor(OR)))] + 
  # geom_smooth(se=FALSE, size = 1.5, span = 1) + 
  geom_line(size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.justification = c(1, 0),
        legend.position = c(0.8, 0.2),) +
  xlab("Minor Allele Frequency") +
  ylab("Power") +
  labs(color = "Odds Ratio") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_discrete(labels = as.character(or))

print(pw_plot)
