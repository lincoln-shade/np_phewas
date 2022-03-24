pacman::p_load(data.table, magrittr, ggplot2)

np <- readRDS("data/np_cl.Rds")
np <- np[, .(ASC, HS, LATE, NPMICRO)
         ][
           !(is.na(ASC) & is.na(HS) & is.na(LATE) & is.na(NACCMICR))
         ]


bar_data <- data.table(NPE = factor(rep(c("HS", "LATE", "NACCMICR"), 2)),
                       NPE_yes = factor(rep(c("Yes", "No"), 3))
                       )
setorder(bar_data, NPE, NPE_yes)

calc_percent_asc <- function(NPE, yes) {
  np[get(NPE) == yes & ASC == 1, .N] / (np[get(NPE) == yes & ASC == 1, .N] + np[get(NPE) == yes & ASC == 0, .N]) * 100
}

bar_data[, percent_asc := c(calc_percent_asc("HS", 0),
                            calc_percent_asc("HS", 1),
                            calc_percent_asc("LATE", 0),
                            calc_percent_asc("LATE", 1),
                            calc_percent_asc("NACCMICR", 0),
                            calc_percent_asc("NACCMICR", 1)
                            )
         ]

# plot
bar_data[, ggplot(.SD, aes(NPE, fill = NPE_yes, weight = percent_asc))] + 
  geom_bar(position = "dodge") +
  theme_minimal() + 
  labs(x = element_blank(),
       y = "% with Moderate or Severe B-ASC",
       fill = "NPE Present") +
  scale_x_discrete(labels = c("HS", "TDP-43", "Micro Infarcts"))
