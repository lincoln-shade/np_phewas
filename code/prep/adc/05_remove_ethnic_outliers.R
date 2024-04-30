##----------------------------------------------
## create PCA plot of 1KG and adc NP subjects
##----------------------------------------------

library(pacman)
p_load(data.table, magrittr, ggplot2)
# groundhog.library("ggplot2", day)

# 1000g (otg) data population data
otg <- fread(
  "/data_global/1000g/integrated_call_samples_v3.20130502.ALL.panel"
  ) %>% 
  .[, sample, super_pop]

# PCA data
merged <- fread("tmp/adc_1000g_merged_pca.eigenvec", header = F) %>% 
  .[, 2:4] %>% 
  setnames(., c("V2", "V3", "V4"), c("sample", "PC1", "PC2")) %>% 
  .[, PC1.norm := scale(PC1)] %>% 
  .[, PC2.norm := scale(PC2)] 

# merge
otg.merged <- merge(merged, otg, by = c("sample"), all = T)

# reverse PCs if needed to make European pop top right of plot & label samples
if (otg.merged[super_pop == "EUR", mean(PC1.norm)] < 0) {
  otg.merged[, `:=`(PC1 = -PC1, PC1.norm = -PC1.norm)]
  }
if (otg.merged[super_pop == "EUR", mean(PC2.norm)] < 0) {
  otg.merged[, `:=`(PC2 = -PC2, PC2.norm = -PC2.norm)]
  }
#otg.merged[is.na(super_pop), super_pop := "adc"]

#ggplot of first 2 normalized PCs
pca.plot <- ggplot(otg.merged, aes(PC1.norm, PC2.norm, color = as.factor(super_pop))) +
  geom_point() +
  ggtitle("PCA of 1000 Genomes and adc participants") 

# create objects for mean values for EUR PC1 and PC2 and create selection criteria based off them

#######################
## set circle radius ##
#######################
radius <- 0.38  #######
#######################

CircleFun <- function(center = c(otg.merged[super_pop == "EUR", 
                                            mean(PC1.norm)],
                                 otg.merged[super_pop == "EUR", 
                                            mean(PC2.norm)]), 
                      r = 1, npoints = 100){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

##-------------------
## Make plot
##-------------------
pca.plot2 <- ggplot(otg.merged, 
                    aes(PC1.norm, PC2.norm, 
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

 
otg.merged[, include := 
             ifelse(
               is.na(super_pop) == T & 
                 InsideCircle(
                   PC1.norm, 
                   PC2.norm, 
                   radius, 
                   c(
                     otg.merged[super_pop == "EUR", mean(PC1.norm)], 
                     otg.merged[super_pop == "EUR", mean(PC2.norm)]
                   )
                 ) == T, 
               TRUE, 
               FALSE
             )
]

adc <- otg.merged[include == T, ]
adc_ids <- fread("tmp/np_ids.tmp", header = F)
adc_ids <- adc_ids[V2 %in% adc$sample]

ggsave(filename = "doc/adc_1000g_pca.png", 
       pca.plot2, units="in", width=7, height=7)

write.table(adc_ids, "data/adc/ids_qced.txt", 
            col.names = F, row.names = F, quote = F)
