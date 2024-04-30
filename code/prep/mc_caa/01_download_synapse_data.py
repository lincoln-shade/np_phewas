import synapseclient
import synapseutils
syn = synapseclient.Synapse()
syn.login()


outdir = "data/mc_caa/"
# files = synapseutils.syncFromSynapse(syn, 'syn17010931', path = outdir + "Metadata/", ifcollision="keep.local") 
# files = synapseutils.syncFromSynapse(syn, 'syn21547862', path = outdir + "Genetic_Variants/", ifcollision="keep.local") 
# files = synapseutils.syncFromSynapse(syn, 'syn10850933', path = outdir + "Gene_Expression/", ifcollision="keep.local") 
