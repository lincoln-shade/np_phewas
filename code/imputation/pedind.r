options(stringsAsFactors = FALSE)
list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(dplyr)

argv <- commandArgs(TRUE)

sample.1000g = read.table("/data_global/1000g/integrated_call_samples_v3.20130502.ALL.panel", h = T)

if (argv[1] != "ALL") {
  
  sample.1000g = sample.1000g %>% filter(super_pop == argv[1]) %>% mutate(code = as.numeric(factor(pop)) + 2) %>% mutate(gender = ifelse(gender=="male", 1, 2)) %>% select(sample, pop, code, gender) %>% rename(pop_c = pop)
 
} else {
  
  sample.1000g = sample.1000g %>% mutate(code = as.numeric(factor(super_pop)) + 2) %>% mutate(gender = ifelse(gender=="male", 1, 2)) %>% select(sample, super_pop, code, gender) %>% rename(pop_c = super_pop)
  
}

fam = read.table(paste(argv[2], ".fam", sep = ""), h = F)

fam = fam %>% rename(sample = V2) %>% left_join(sample.1000g, by = "sample") %>% mutate(code = ifelse(is.na(code), 1, code))  %>% select(V1, sample, V3, V4, V5, code, pop_c)
system(paste("mv ", argv[2], ".fam ", argv[2], ".fam_delete", sep = ""))

write.table(fam, paste(argv[2], ".fam", sep = ""), quote = F, row.names = F, col.names = F)
