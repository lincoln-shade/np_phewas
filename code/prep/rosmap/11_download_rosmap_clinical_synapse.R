library(synapser)
library(synapserutils)

synapser::synLogin()

# download rosmap clinical data set
syncFromSynapse('syn3191087', 
                path="data/synapse", 
                ifcollision="keep.local")