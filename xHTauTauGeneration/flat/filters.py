import ROOT

from rootpy.tree.filtering import *
from rootpy.extern.hep import pdg

from math import *
from .units import GeV
from . import utils
from . import log; log = log[__name__]
from .taudecay import TauDecay
from .models import FourMomentum

class TrueTaus(EventFilter):

    def __init__(self, **kwargs):
        super(TrueTaus, self).__init__(**kwargs)

    def passes(self, event):
        higgs = event.higgs[0]
        event.taus = [tau for tau in higgs.iter_children() 
                      if abs(tau.pdgId) == pdg.tau and tau.status == 2]
        return len(event.taus) == 2

class TrueJets(EventFilter):

    def __init__(self, **kwargs):
        super(TrueJets, self).__init__(**kwargs)
        
    def passes(self, event):
        # The number of anti kt R = 0.4 truth jets with pT>25 GeV, not
        # originating from the decay products of the Higgs boson.
        # Start from the AntiKt4Truth collection. Reject any jet with pT<25
        # GeV. Reject any jet withing dR < 0.4 of any electron, tau, photon or
        # parton (directly) produced in the Higgs decay.
        event.jets.select(lambda jet: jet.pt >= 25 * GeV and \
                              not any([tau for tau in event.taus if
                                       utils.dR(jet.eta, jet.phi,
                                                tau.decay.fourvect_vis.Eta(), 
                                                tau.decay.fourvect_vis.Phi()) < 0.4]))
        return True


class ClassifyDecay(EventFilter):
    def __init__(self, tree, **kwargs):
        super(ClassifyDecay, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        tau1, tau2 = event.taus
        tau1.decay = TauDecay(tau1)
        tau2.decay = TauDecay(tau2)
        if tau1.decay.hadronic and tau2.decay.hadronic:
            self.tree.hadhad = 1
        elif tau1.decay.hadronic and tau2.decay.leptonic or\
                tau1.decay.leptonic and tau2.decay.hadronic:
            self.tree.lephad = 1
        elif tau1.decay.leptonic and tau2.decay.leptonic:
            self.tree.leplep = 1

        return True
    

class Higgs(EventFilter):

    def __init__(self, **kwargs):
        super(Higgs, self).__init__(**kwargs)
        self.status = (62, 195)
        
    def passes(self, event):
        pt = 0
        higgs = None
        status = self.status
        # find the Higgs
        for mc in event.mc:
            if mc.pdgId == 25 and mc.status in status:
                higgs = mc
                break
        if higgs is None:
            raise RuntimeError("Higgs not found!")
        event.higgs = [higgs]
        return True

