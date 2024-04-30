library(synapser)
library(synapserutils)

synLogin()
out_path <- "data/synapse/rna_seq/"
syncFromSynapse('syn21088596', path=out_path, ifcollision="keep.local")
syncFromSynapse('syn21323366', path=out_path, ifcollision="keep.local")
syncFromSynapse('syn3505732', path=out_path, ifcollision="keep.local")
syncFromSynapse('syn3505724', path=out_path, ifcollision="keep.local")
