options(stringsAsFactors = FALSE)
options(warn=-1)

list.of.packages <- c("ggplot2", "RColorBrewer", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

argv <- commandArgs(TRUE)
dat = read.csv(argv[1], na.strings = c(""))
filename = tools::file_path_sans_ext(argv[1])

legend.order = dat %>% dplyr::select(Superpop) %>% unique() %>% arrange(Superpop)
legend.order = c("Subject", legend.order$Superpop[-nrow(legend.order)])
dat$Superpop[is.na(dat$Superpop)] = "Subject"

dat$Superpop = factor(dat$Superpop, legend.order)


my.cols <- brewer.pal(length(legend.order), "Set1")
my.cols = c("#000000", my.cols[1:(length(unique(dat$Superpop))-1)])

p = ggplot()+
  geom_point(data = dat, aes(x = V1, y = V2, color = Superpop), size = 0.5)+
  scale_color_manual(values = my.cols)+
  theme(legend.key = element_rect(colour = "white", fill = "white"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill = NA, colour = "black"),
        text = element_text(size=16),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(vjust=-0.2, size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_blank())+
  labs(x = "1st component", y = "2nd component")

ggsave(paste(filename, ".ancestry-plot.pdf", sep = ""), plot = p, width=10, height=8)
write.csv(dat, argv[1], quote = F, row.names = F)
