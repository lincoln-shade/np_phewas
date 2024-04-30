import synapseclient 
import synapseutils
syn = synapseclient.Synapse()
syn.login()


outdir = "data/msbb/"

# files = synapseutils.syncFromSynapse(syn, 'syn21447661', path = outdir, ifcollision="keep.local") 
# files = synapseutils.syncFromSynapse(syn, 'syn51138026', path = outdir, ifcollision="keep.local") 
# files = synapseutils.syncFromSynapse(syn, 'syn22157601', path = outdir, ifcollision="keep.local")
# files = synapseutils.syncFromSynapse(syn, 'syn22447899', path = outdir, ifcollision="keep.local")
# files = synapseutils.syncFromSynapse(syn, 'syn6101474', path = outdir, ifcollision="keep.local")
files = synapseutils.syncFromSynapse(syn, 'syn21893059', path = outdir, ifcollision="keep.local")
