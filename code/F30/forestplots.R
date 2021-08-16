#========================================
# Forest plot from original B-ASC study
#========================================

pacman::p_load(data.table, magrittr, forestplot)

elovl4_vals <- structure(list(
  mean = c(NA, NA, 1.45, 1.06, 4.75, 1.37),
  lower = c(NA, NA, 1.27, 0.84, 1.42, 1.23),
  upper = c(NA, NA, 1.66, 1.31, 15.91, 1.53),
  .Names  = c("mean", "lower", "upper"),
  row.names = c(NA, 6L),
  class = "data.frame"
))

elovl4_text <- cbind(
  c("", "Study", "NACC", "ROSMAP", "ADNI", "Summary"),
  c("", "Variant", "rs2603462", "rs2603462", "rs2603462", "rs2603462"),
  c("", "OR", "1.45", "1.06", "4.75", "1.37")
)

# save 615x30
forestplot(elovl4_text, 
           elovl4_vals$mean,
           elovl4_vals$lower,
           elovl4_vals$upper,
           new_page = FALSE,
           is.summary = c(TRUE,TRUE,rep(FALSE,3),TRUE),
           clip = c(0.3, 3),
           # xlog = FALSE, 
           xticks = c(0, 1.0, 2.0, 3),
           grid = 1,
           col = fpColors(box = c("royalblue"),
                          line = "darkblue",
                          summary = c("royalblue")))

sorcs3_vals <- structure(list(
  mean = c(NA, NA, 1.57, 1.61, 6.49, 1.60),
  lower = c(NA, NA, 1.29, 1.14, 2.07, 1.35),
  upper = c(NA, NA, 1.92, 2.26, 8.80, 1.89),
  .Names  = c("mean", "lower", "upper"),
  row.names = c(NA, 6L),
  class = "data.frame"
))

sorcs3_text <- cbind(
  c("", "Study", "NACC", "ROSMAP", "ADNI", "Summary"),
  c("", "Variant", "rs7902929", "rs7902929", "rs7902929", "rs7902929"),
  c("", "OR", "1.57", "1.61", "26.49", "1.60")
)

# save 615x30
forestplot(sorcs3_text, 
           sorcs3_vals$mean,
           sorcs3_vals$lower,
           sorcs3_vals$upper,
           new_page = FALSE,
           is.summary = c(TRUE,TRUE,rep(FALSE,3),TRUE),
           clip = c(0.3, 3),
           # xlog = FALSE, 
           xticks = c(0, 1.0, 2.0, 3),
           grid = 1,
           col = fpColors(box = c("royalblue"),
                          line = "darkblue",
                          summary = c("royalblue")))
