import math
from decorators import cached_property

from rootpy.vector import LorentzVector, Vector3
from rootpy.extern.hep import pdg

from . import log; log = log[__name__]
from .utils import dR, et2pt
from .units import GeV

"""
This module contains "mixin" classes for adding
functionality to Tree objects ("decorating" them).
"""

__all__ = [
    'FourMomentum',
    'FourMomentumMeV',
    'MCParticle',
    'MCTauParticle',
    'MCTauFourMomentum',
]


class MatchedObject(object):

    def __init__(self):
        self.matched = False
        self.matched_dR = 9999.
        #self.matched_collision = False
        self.matched_object = None

    def matches(self, other, thresh=.2):
        return self.dr(other) < thresh

    def dr(self, other):
        return dR(self.eta, self.phi, other.eta, other.phi)

    def dr_vect(self, other):
        return dR(self.eta, self.phi, other.Eta(), other.Phi())

    def angle_vect(self, other):
        return self.fourvect.Angle(other)

    def matches_vect(self, vect, thresh=.2):
        return self.dr_vect(vect) < thresh


class FourMomentum(MatchedObject):

    def __init__(self):
        self.fourvect_boosted = LorentzVector()
        super(FourMomentum, self).__init__()

    @cached_property
    def fourvect(self):
        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.m)
        return vect

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "%s (m: %.3f MeV, pt: %.1f MeV, eta: %.2f, phi: %.2f)" % \
            (self.__class__.__name__,
             self.m,
             self.pt,
             self.eta,
             self.phi)


class FourMomentumMeV(object):

    def __init__(self):
        self.fourvect_boosted = LorentzVector()

    @cached_property
    def fourvect(self):
        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt*GeV, self.eta, self.phi, self.m)
        return vect

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "%s (m: %.3f MeV, pt: %.1f MeV, eta: %.2f, phi: %.2f)" % \
            (self.__class__.__name__,
             self.m,
             self.pt,
             self.eta,
             self.phi)

class MCTauFourMomentum(FourMomentum):

    @cached_property
    def fourvect_vis(self):
        vect = LorentzVector()
        try:
            vect.SetPtEtaPhiM(
                et2pt(self.vis_Et, self.vis_eta, self.vis_m),
                self.eta, self.phi, self.m)
        except ValueError:
            log.warning("DOMAIN ERROR ON TRUTH 4-VECT: "
                        "Et: {0} eta: {1} m: {2}".format(
                            self.vis_Et, self.vis_eta, self.vis_m))
            vect.SetPtEtaPhiM(0, self.eta, self.phi, self.m)
        return vect



class MCParticle(FourMomentum):

    def __init__(self):
        self._particle = pdg.GetParticle(self.pdgId)
        FourMomentum.__init__(self)

    @cached_property
    def num_children(self):
        return len(self.child_index)

    @cached_property
    def num_parents(self):
        return len(self.parent_index)

    def get_child(self, index):
        index = self.child_index[index]
        return getattr(self.tree, self.name)[index]

    def get_parent(self, index):
        index = self.parent_index[index]
        return getattr(self.tree, self.name)[index]

    def iter_children(self):
        try:
            for child in self.child_index:
                yield getattr(self.tree, self.name)[child]
        except GeneratorExit:
            pass

    def iter_parents(self):
        try:
            for parent in self.parent_index:
                yield getattr(self.tree, self.name)[parent]
        except GeneratorExit:
            pass

    def traverse_children(self):
        try:
            for child in self.iter_children():
                yield child
                for desc in child.traverse_children():
                    yield desc
        except GeneratorExit:
            pass

    def traverse_parents(self):
        try:
            for parent in self.iter_parents():
                yield parent
                for ancestor in parent.traverse_parents():
                    yield ancestor
        except GeneratorExit:
            pass

    def is_stable(self):
        return self.status == 1

    @cached_property
    def first_self(self):
        for parent in self.iter_parents():
            if parent.pdgId == self.pdgId:
                return parent.first_self
        return self

    @cached_property
    def last_self(self):
        for child in self.iter_children():
            if child.pdgId == self.pdgId:
                return child.last_self
        return self

    @cached_property
    def final_state(self):
        if self.is_stable():
            return [self]
        return [particle for particle in self.traverse_children()
                if particle.is_stable()]

    @cached_property
    def fourvect(self):
        vect = LorentzVector()
        vect.SetPtEtaPhiM(
                self.pt,
                self.eta,
                self.phi,
                self.m)
        #       self._particle.Mass() * GeV)
        return vect


    def export_graphvis(self, out_file=None):
        def particle_to_str(particle):
            return ('%s\\n'
                    'mass = %.3f MeV\\n'
                    'pt = %.3f GeV\\n'
                    'eta = %.2f\\n'
                    'status = %d') % (
                    particle._particle.GetName(),
                    #particle._particle.Mass() * GeV,
                    particle.m,
                    particle.pt / GeV,
                    particle.eta,
                    particle.status)

        def recurse(particle, parent=None):
            out_file.write('%d [label="%s"] ;\n' % (
                particle.barcode, particle_to_str(particle)))

            if parent is not None:
                # Add edge to parent
                out_file.write('%d -> %d ;\n' % (
                    parent.barcode,
                    particle.barcode))

            # recurse on children
            for child in particle.iter_children():
                recurse(child, particle)

        close_file = True
        if out_file is None:
            out_file = open('event.dot', 'w')
        elif isinstance(out_file, basestring):
            out_file = open(out_file, 'w')
        else:
            close_file = False

        out_file.write('digraph Tree {\n')
        out_file.write('size="7.5,10" ;\n')
        out_file.write('orientation=landscape ;\n')
        recurse(self, None)
        out_file.write('}')

        if close_file:
            out_file.close()
        else:
            return out_file

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return ("%s ("
                "status: %d, "
                "m: %.3f MeV, pt: %.1f GeV, eta: %.2f, phi: %.2f" %
            (self._particle.GetName(),
             self.status,
             #self._particle.Mass() * GeV,
             self.m,
             self.pt / GeV,
             self.eta, self.phi))

class MCTauParticle(MCParticle):
    
    def __init__(self):
        super(MCTauParticle, self).__init__()
        self.decay = None



