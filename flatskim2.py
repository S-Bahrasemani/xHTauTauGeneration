import ROOT

import ROOT
ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')

if not ROOT.xAOD.Init().isSuccess():
    raise RuntimeError



from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel, TreeBuffer
from rootpy.extern.argparse import ArgumentParser
from rootpy.io import root_open
from rootpy import stl

from flat.models import *
from flat.filters import *
from flat.objects import define_objects
import flat.branches as branches
from flat import log; log = log[__name__]
parser = ArgumentParser()
parser.add_argument('input', type=str, help='input file name')
parser.add_argument('output', type=str, help='output file name')
parser.add_argument('--tree-name', type=str, default='truth', help='Input tree name')
args = parser.parse_args()

output = root_open(args.output, 'recreate')
output.cd()
model = get_model()
outtree = Tree(name='Tree', model=model)
outtree.define_object(name='tau1', prefix='tau1_')
outtree.define_object(name='tau2', prefix='tau2_')
outtree.define_object(name='jet1', prefix='jet1_')
outtree.define_object(name='jet2', prefix='jet2_')
outtree.define_object(name='higgs', prefix='resonance_')
outtree.define_object(name='met', prefix='MET_')

log.info(model)
def mc_weight_count(event):
    return event.mc_event_weight
count_funcs = {
    'mc_weight': mc_weight_count,
}
event_filters = EventFilterList([
        Higgs(count_funcs=count_funcs),
        TrueTaus(count_funcs=count_funcs),
        ClassifyDecay(
            count_funcs=count_funcs,
            tree=outtree),
        TrueJets(count_funcs=count_funcs),
        ])

files = [args.input]

# peek at first tree to determine which branches to exclude
# with root_open(files[0]) as test_file:
#     test_tree = test_file.Get(args.tree_name)
#     ignore_branches = test_tree.glob(branches.REMOVE)

# chain = TreeChain(
#     args.tree_name, 
#     files=files, 
#     filters=event_filters,
#     cache=True,
#     cache_size=50000000,
#     learn_entries=100,
#     ignore_branches=ignore_branches)

print 30 * '*'
ch = ROOT.TChain('CollectionTree')
ch.Add(files[0])
tree = ROOT.xAOD.MakeTransientTree(ch, ROOT.xAOD.TEvent.kBranchAccess)
for i in xrange(tree.GetEntries()):
    tree.GetEntry(i)
    print 20 * '-'
    print '######  EVENT {0} ######'.format(i) 
    print tree.EventInfo.runNumber()
    
    muons = tree.Muons
    if muons.size() > 0:
        print 'First muon pT = ', muons[0].pt(),  muons[0].e(),  muons[0].phi(),  muons[0].eta()

    if not muons.size() > 1: 
        continue
    print muons.size()
    
    tau1 , mu2 = muons[0], muons[1]
    # outtree.tau1_pt = mu1.pt()
    # outtree.tau2_pt = mu2.pt()
    # outtree.tau1 = mu1

    TrueTau.set_full(outtree.tau1, tau1)


    outtree.fill(reset=True)
output.cd()
outtree.FlushBaskets()
outtree.Write()
    
#    TrueTauBlock.set(outtree, mu1, mu2)



# define_objects(chain)


# for event in chain:

#     outtree.runnumber = event.RunNumber
#     outtree.evtnumber = event.EventNumber
#     outtree.weight = event.mc_event_weight
#     # sort taus and jets in decreasing order by pT
#     event.taus.sort(key=lambda tau: tau.decay.fourvect_vis.Pt(), reverse=True)
#     event.jets.sort(key=lambda jet: jet.pt, reverse=True)

#     # Set variables describing the two taus 
#     # and the ditau system
#     tau1, tau2 = event.taus
#     FourMomentum.set(outtree.higgs, event.higgs[0])
#     outtree.Fill(reset=-1)
   
#     # MET = tau1.decay.fourvect_missing + tau2.decay.fourvect_missing
#     jets = list(event.jets)
#     outtree.numJets= len(jets)
#     if len(jets) >=2:
#         jet1, jet2 = jets[:2]
#         TrueTauBlock.set(outtree, tau1, tau2, jet1, jet2)
#         TrueJetBlock.set(outtree, jet1, jet2)
        
#         jet1 = jets[0]
#         TrueTauBlock.set(outtree, tau1, tau2, jet1)
#         TrueJetBlock.set(outtree, jet1, jet2)

#     elif len(jets) == 1:
#         jet1 = jets[0]
#         TrueTauBlock.set(outtree, tau1, tau2, jet1)
#         TrueJetBlock.set(outtree, jet1, None)
    
#     else:
#         TrueTauBlock.set(outtree, tau1, tau2)

#     outtree.sum_pt_full = sum(
#         [outtree.tau1.pt, outtree.tau2.pt] + 
#         [jet.pt for jet in jets])
    

#             # vector sum pT with all selected jets and MET
#     outtree.vector_sum_pt_full = sum(
#         [tau1.fourvect, tau2.fourvect] +
#         [jet.fourvect for jet in jets] + [tau1.decay.fourvect_missing + tau2.decay.fourvect_missing]).Pt()



#     outtree.true_resonance_pt = outtree.higgs.pt
