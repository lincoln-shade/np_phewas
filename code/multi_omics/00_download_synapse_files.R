library(synapser)
library(synapserutils)

synLogin()
files <- syncFromSynapse('syn3157275',
                         path="data/synapse/DNA_methylation",
                         ifcollision="keep.local")

clin_file <- syncFromSynapse('syn3191087', 
                         path="data/synapse", 
                         ifcollision="keep.local")
