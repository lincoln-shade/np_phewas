library(forestplot)

#list.files(path = "/home/kau232/locuszoom_test/lincoln/Forest_plot/")

cil = function(beta, se) {
  cil = exp(beta - 1.96 * se)
}

ciu = function(beta, se) {
  ciu = exp(beta + 1.96 * se)
}
# NACC 
# apoe33 -0.163  0.056
# apoe34 -0.248 0.055
# apoe44 -0.190112   0.103950
# apoe24 -0.56 0.246
# apoe23 -0.4 0.172
# Create a data frame for forestplot
forest_data <- data.frame(
  labeltext = c("APOE 3/3", "APOE 3/4", "APOE 4/4", "APOE 2/4", "APOE 2/3"),
  mean = c(exp(-0.163), exp(-0.248), exp(-0.19), exp(-0.56), exp(-0.4)),
  lower = c(cil(-0.163, 0.056), cil(-0.248, 0.055), cil(-0.190112, 0.103950), cil(-0.56, 0.246), cil(-0.4, 0.172)),
  upper = c(ciu(-0.163, 0.056), ciu(-0.248, 0.055), ciu(-0.190112, 0.103950), ciu(-0.56, 0.246), ciu(-0.4, 0.172))
)

# Create a forest plot with labels
forestplot(
  mean = forest_data$mean,
  lower = forest_data$lower,
  upper = forest_data$upper,
  labeltext = forest_data$labeltext,
  title = "NACC: Forest Plot of APOE-stratified rs7247551 assocation with CAA",
  xlab = "Odds Ratio with 95% CI",
  col = fpColors(box = "blue", line = "darkblue", summary = "royalblue"),
  legend.pos = "bottom"  # Adjust legend position
)

# ROSMAP 
# apoe33 -0.23656    0.10231
# apoe34_44_24 -0.38474    0.15768
# apoe23 -0.42770    0.28593
# Create a data frame for forestplot
forest_data <- data.frame(
  labeltext = c("APOE 3/3", "APOE 4 carrier", "APOE 2/3"),
  mean = c(exp(-0.23656), exp(-0.38474), exp(-0.42770)),
  lower = c(cil(-0.23656, 0.10231), cil(-0.38474, 0.15768), cil(-0.42770, 0.28593)),
  upper = c(ciu(-0.23656, 0.10231), ciu(-0.38474, 0.15768), ciu(-0.42770, 0.28593))
)

# Create a forest plot with labels
forestplot(
  mean = forest_data$mean,
  lower = forest_data$lower,
  upper = forest_data$upper,
  labeltext = forest_data$labeltext,
  title = "ROSMAP: Forest Plot of APOE-stratified rs7247551 assocation with CAA",
  xlab = "Odds Ratio with 95% CI",
  col = fpColors(box = "blue", line = "darkblue", summary = "royalblue"),
  legend.pos = "bottom"  # Adjust legend position
)

# ACT 
# apoe33 -0.12013    0.14586
# apoe34_44_24 -0.07230    0.20464
# apoe23 0.36319    0.62281
# Create a data frame for forestplot
forest_data <- data.frame(
  labeltext = c("APOE 3/3", "APOE 4 carrier", "APOE 2/3"),
  mean = c(exp(-0.12013), exp(-0.07230), exp(0.36319)),
  lower = c(cil(-0.12013, 0.14586), cil(-0.07230, 0.20464), cil(0.36319, 0.62281)),
  upper = c(ciu(-0.12013, 0.14586), ciu(-0.07230, 0.20464), ciu(0.36319, 0.62281))
)

# Create a forest plot with labels
forestplot(
  mean = forest_data$mean,
  lower = forest_data$lower,
  upper = forest_data$upper,
  labeltext = forest_data$labeltext,
  title = "ACT: Forest Plot of APOE-stratified rs7247551 assocation with CAA",
  xlab = "Odds Ratio with 95% CI",
  col = fpColors(box = "blue", line = "darkblue", summary = "royalblue"),
  legend.pos = "bottom"  # Adjust legend position
)


# Create a forest plot with ggplot2
ggplot(ggplot_data, aes(x = mean, y = labeltext)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "blue") +
  labs(title = "CAA results for rs7247551, by study",
       x = "Odds Ratio with 95% CI",
       y = "Study") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.6, 1.14)) +  # Adjust x-axis limits
  theme(axis.title.x = element_text(hjust = 1))  # Adjust x-axis label position
# Convert "Study" to a factor with custom order
ggplot_data$labeltext <- factor(ggplot_data$labeltext, levels = c("Pooled", "ACT", "ROSMAP", "NACC"))

ggplot_data$labeltext <- factor(ggplot_data$labeltext, levels = c("Pooled", "ACT", "ROSMAP", "NACC"))



# Create a forest plot with ggplot2
forest_p1 <- ggplot(ggplot_data, aes(x = mean, y = labeltext)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "blue") +
  labs(title = "CAA results for rs7247551, by study",
       x = "Odds Ratio with 95% CI",
       y = "Study") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.6, 1.14)) +  # Adjust x-axis limits
  theme(axis.title.x = element_text(hjust = 0.5),  # Center x-axis label
        plot.title = element_text(hjust = 0.5)) +  # Center the plot title
  geom_vline(xintercept = 1, linetype = "dotted", color = "blue") +  # Add a dotted line in Y-axis
  theme(plot.background = element_rect(fill = "white"))  # Set the background color


#save
ggsave("/home/kau232/locuszoom_test/lincoln/Forest_plot/Forest_plot.png", plot = forest_p1, width = 10, height = 6)

