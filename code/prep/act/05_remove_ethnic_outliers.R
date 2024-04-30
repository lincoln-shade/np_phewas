##----------------------------------------------
## create PCA plot of 1KG and adc NP subjects
##----------------------------------------------


library(data.table)
library(magrittr)
library(ggplot2)
cargs <- commandArgs(trailingOnly = TRUE)
merged_eigenvec_file <- cargs[1]
study <- cargs[2]
study_ids <- fread(paste0("tmp/", study, "_qc1.tmp.fam"), header = F)

# 1000g (otg) data population data
otg <- fread("/data_global/1000g/integrated_call_samples_v3.20130502.ALL.panel") %>% 
  .[, sample, super_pop]

# PCA data
merged <- fread(merged_eigenvec_file, header = F) %>% 
  .[, 2:4] %>% 
  setnames(., c("V2", "V3", "V4"), c("sample", "PC1", "PC2")) %>% 
  .[, PC1.norm := scale(PC1)] %>% 
  .[, PC2.norm := scale(PC2)] 

# merge
otg_merged <- merge(merged, otg, by = c("sample"), all = TRUE)

# reverse PCs if needed to make European pop top right of plot & label samples
if (otg_merged[super_pop == "EUR", mean(PC1.norm)] < 0) {
  otg_merged[, `:=`(PC1 = -PC1, PC1.norm = -PC1.norm)]
  }

if (otg_merged[super_pop == "EUR", mean(PC2.norm)] < 0) {
  otg_merged[, `:=`(PC2 = -PC2, PC2.norm = -PC2.norm)]
  }
#otg_merged[is.na(super_pop), super_pop := "adc"]

#ggplot of first 2 normalized PCs
pca.plot <- ggplot(otg_merged, 
                   aes(PC1.norm, 
                       PC2.norm, 
                       color = as.factor(super_pop))) +
  geom_point() +
  ggtitle(paste0("PCA of 1000 Genomes and ", study, " participants")) 

# create objects for mean values for EUR PC1 and PC2 
# and create selection criteria based off them

#######################
## set circle radius ##
#######################
radius <- as.numeric(cargs[3])
#######################

eur_pc1_mean <- otg_merged[super_pop == "EUR", mean(PC1.norm)]
eur_pc2_mean <- otg_merged[super_pop == "EUR", mean(PC2.norm)]

CircleFun <- function(center = c(eur_pc1_mean,
                                 eur_pc2_mean), 
                      r = 1, 
                      npoints = 100){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

##-------------------
## Make plot
##-------------------
pca.plot2 <- ggplot(otg_merged, aes(PC1.norm, 
                                    PC2.norm, 
                                    color = as.factor(super_pop))) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_point() +
  #ggtitle("PCA of 1000 Genomes and adc Participants") +
  xlab("PC1") +
  ylab("PC2") +
  geom_path(data = CircleFun(r=radius), aes(x, y), color = "red")

pca.plot2

InsideCircle <- function(x, y, r, center = c(0, 0)) {
  ifelse((x - center[1])^2 + (y - center[2])^2 <= r, TRUE, FALSE)
}

otg_merged <- 
  otg_merged[, include := 
               ifelse(
                 is.na(super_pop) == T & 
                   InsideCircle(PC1.norm, 
                                PC2.norm, 
                                radius, 
                                c(eur_pc1_mean, 
                                  eur_pc2_mean)
                   ) == T, 
                 TRUE, 
                 FALSE)
             ]

keep <- otg_merged[include == T, sample]
study_ids <- study_ids[V2 %in% keep, .(V1, V2)]

ggsave(filename = paste0("doc/", study, "_1000g_pca.png"), 
       pca.plot2, 
       units="in", 
       width=7, 
       height=7)

fwrite(study_ids, 
            file = paste0("tmp/", study, "_qced.txt"), 
            col.names = F, 
            sep = " ", 
            quote = F)

