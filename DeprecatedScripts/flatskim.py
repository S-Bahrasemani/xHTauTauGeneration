from ROOT import gROOT

gROOT.SetBatch( True )

import ROOT
import os
import sys
import math
import random
import glob

ROOT.gROOT.SetBatch( True )

from math import sqrt
from array import array

if not os.path.exists( 'loadvec_cc.so' ):
    hndl = open( 'loadvec.cc', 'w' )
    hndl.write( '#include "TLorentzVector.h"\n' )
    hndl.write( '#include "TH1F.h"\n' )
    hndl.write( '#ifdef __CINT__\n' )
    hndl.write( '#pragma link C++ class vector<vector<int> >+;\n' )
    hndl.write( '#pragma link C++ class vector<vector<double> >+;\n' )
    hndl.write( '#pragma link C++ class vector<TLorentzVector>+;\n' )
    hndl.write( '#endif\n' )
    hndl.close()

ROOT.gROOT.ProcessLine( '.L loadvec.cc+' )

#HiggsMass=91
#ifilenames='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/ntup.EVNT_300000.root'
#ofilename='/cluster/data03/sbahrase/BrtStudies/PracticeDesk/ntup.EVNT_300000_F.root'
#ifilenames = i = 0

i =0
ifilenames = sys.argv[1]
ofilename = sys.argv[2] 


print sys.argv


itree = ROOT.TChain( 'truth', 'truth' )
ifiles = glob.glob( ifilenames )
for fn in ifiles:
    sys.stdout.write( '\t' + fn + '\n' )
    itree.Add( fn )
del ifiles
sys.stdout.flush()

nevts = itree.GetEntries()
sys.stdout.write( 'Will run over ' + str(nevts) + ' events ...\n' )
sys.stdout.flush()

itree.SetBranchStatus( '*', 0 )
itree.SetBranchStatus( 'RunNumber', 1 )
itree.SetBranchStatus( 'EventNumber', 1 )
itree.SetBranchStatus( 'mc_*', 1 )
itree.SetBranchStatus( 'mcevt_*', 1 )
itree.SetBranchStatus( 'MET_Truth_Int_*', 1 )
# itree.SetBranchStatus( 'trueTau_*', 1 )
itree.SetBranchStatus( 'tau_*', 1 )

ofile = ROOT.TFile( ofilename, 'recreate' )
ofile.cd()
otree = ROOT.TTree( 'Tree', 'Tree' )
        
#----------------------------------------
        
def delta_phi(phi1,phi2):
    dphi = phi2-phi1
    while dphi >  math.pi: dphi -= 2*math.pi
    while dphi < -math.pi: dphi += 2*math.pi
    return dphi

def tau_relative_momentum_resolution(pt,eta):
    # Takes pt in GeV, returns the relative resolution (sigma of a gaussian) on the tau momentum.
    # Numbers taken from https://cds.cern.ch/record/1544036/files/ATLAS-CONF-2013-044.pdf
    # corresponding to multi-prong taus. The resolutions for 1-prong taus are better, so I took the worst scenario.
    p = pt*math.cosh(eta)
    if abs(eta) < 0.3:
        a = 0.85
        b = 0.0
        c = 0.01
    elif abs(eta) < 0.8:
        a = 0.79
        b = 0.0
        c = 0.03
    elif abs(eta) < 1.3:
        a = 1.02
        b = 0.0
        c = 0.04
    elif abs(eta) < 1.6:
        a = 1.24
        b = 0.0
        c = 0.06
    else:
        a = 1.27
        b = 0.0
        c = 0.01
    return math.sqrt(pow(a/math.sqrt(p),2)+pow(b/p,2)+pow(c,2))
    
def electron_relative_energy_resolution(E,eta):
    # Takes E in GeV, returns the relative resolution on the electron energy.
    # Numbers taken from https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection/tags/ElectronPhotonFourMomentumCorrection-00-00-15/data/resolutionFit_electron.root. Thanks to Quentin!
    eta_bins = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.37,1.52,1.56,1.62,1.72,1.82,1.92,2.02,2.12,2.22,2.32,2.4]
    sampl = [0.0922468, 0.0946236, 0.0986818, 0.121886, 0.155332, 0.180168, 0.23041, 0.539857, 0.500733, 0.315278, 0.420421, 0.242754, 0.216817, 0.197961, 0.205227, 0.202158, 0.196951, 0.245199]
    noise = [0.264571, 0.213702, 0.256527, 0.397013, 0.700537, 0.87247, 0.967603, 0.05, 3.3815, 2.45919, 2.39302, 1.06174, 0.05, 0.05, 0.05, 0.274303, 0.0500002, 0.05]
    const = [0.00442843, 0.00356228, 0.00340941, 0.002, 0.002, 0.002, 0.002, 0.01, 0.002, 0.002, 0.002, 0.002, 0.00637377, 0.00272288, 0.00228653, 0.002, 0.002, 0.00273704]
    eta = abs(eta)
    if eta >= eta_bins[-1]:
        return 1
    etabin = 0
    for i in range(len(eta_bins)-1):
        if eta_bins[i] <= eta < eta_bins[i+1]:
            etabin = i
            break
    a = noise[etabin]
    b = sampl[etabin]
    c = const[etabin]
    return math.sqrt(pow(a/E,2)+pow(b/math.sqrt(E),2)+pow(c,2))
    
def met_energy_resolution(SumET):
    # Takes E in GeV, returns the relative resolution on the electron energy.
    # Numbers for c taken from https://cds.cern.ch/record/1335395/files/ATL-COM-PHYS-2011-263.pdf, table 14, 
    return 0.5*math.sqrt(SumET)

#----------------------------------------
                
br_runnumber = array( 'i', [0] )
br_evtnumber = array( 'i', [0] )

br_weight = array( 'd', [0] )
                
br_bad_event = array( 'i', [0] )
                
br_n_photons = array( 'i', [0] )

br_boson_pt    = array( 'd', [0] )
br_boson_eta   = array( 'd', [0] )
br_boson_phi   = array( 'd', [0] )
br_boson_E     = array( 'd', [0] )
br_boson_m     = array( 'd', [0] )
br_boson_mtrue = array( 'd', [0] )

br_taup_pt     = array( 'd', [0] )
br_taup_eta    = array( 'd', [0] )
br_taup_phi    = array( 'd', [0] )
br_taup_E      = array( 'd', [0] )
br_taup_m      = array( 'd', [0] )
br_taup_mtrue  = array( 'd', [0] )

br_taum_pt     = array( 'd', [0] )
br_taum_eta    = array( 'd', [0] )
br_taum_phi    = array( 'd', [0] )
br_taum_E      = array( 'd', [0] )
br_taum_m      = array( 'd', [0] )
br_taum_mtrue  = array( 'd', [0] )

br_ditau_pt    = array( 'd', [0] )
br_ditau_eta   = array( 'd', [0] )
br_ditau_phi   = array( 'd', [0] )
br_ditau_E     = array( 'd', [0] )
br_ditau_m     = array( 'd', [0] )
br_ditau_mtrue = array( 'd', [0] )

br_tau1_decay_nutau_pt    = array( 'd', [0] )
br_tau1_decay_nutau_eta   = array( 'd', [0] )
br_tau1_decay_nutau_phi   = array( 'd', [0] )
br_tau1_decay_nutau_E     = array( 'd', [0] )
br_tau1_decay_nutau_m     = array( 'd', [0] )
br_tau1_decay_nutau_mtrue = array( 'd', [0] )

br_tau1_decay_lep_pt      = array( 'd', [0] )
br_tau1_decay_lep_eta     = array( 'd', [0] )
br_tau1_decay_lep_phi     = array( 'd', [0] )
br_tau1_decay_lep_E       = array( 'd', [0] )
br_tau1_decay_lep_m       = array( 'd', [0] )
br_tau1_decay_lep_mtrue   = array( 'd', [0] )

br_tau1_decay_nulep_pt    = array( 'd', [0] )
br_tau1_decay_nulep_eta   = array( 'd', [0] )
br_tau1_decay_nulep_phi   = array( 'd', [0] )
br_tau1_decay_nulep_E     = array( 'd', [0] )
br_tau1_decay_nulep_m     = array( 'd', [0] )
br_tau1_decay_nulep_mtrue = array( 'd', [0] )

br_tau2_decay_nutau_pt    = array( 'd', [0] )
br_tau2_decay_nutau_eta   = array( 'd', [0] )
br_tau2_decay_nutau_phi   = array( 'd', [0] )
br_tau2_decay_nutau_E     = array( 'd', [0] )
br_tau2_decay_nutau_m     = array( 'd', [0] )
br_tau2_decay_nutau_mtrue = array( 'd', [0] )

br_tau2_decay_lep_pt      = array( 'd', [0] )
br_tau2_decay_lep_eta     = array( 'd', [0] )
br_tau2_decay_lep_phi     = array( 'd', [0] )
br_tau2_decay_lep_E       = array( 'd', [0] )
br_tau2_decay_lep_m       = array( 'd', [0] )
br_tau2_decay_lep_mtrue   = array( 'd', [0] )

br_tau2_decay_nulep_pt    = array( 'd', [0] )
br_tau2_decay_nulep_eta   = array( 'd', [0] )
br_tau2_decay_nulep_phi   = array( 'd', [0] )
br_tau2_decay_nulep_E     = array( 'd', [0] )
br_tau2_decay_nulep_m     = array( 'd', [0] )
br_tau2_decay_nulep_mtrue = array( 'd', [0] )

br_leplep = array( 'i', [0] )
br_lephad = array( 'i', [0] )
br_hadhad = array( 'i', [0] )

br_met_et   = array( 'd', [0] )
br_lep1_pt  = array( 'd', [0] )
br_lep1_eta = array( 'd', [0] )
br_lep2_pt  = array( 'd', [0] )
br_lep2_eta = array( 'd', [0] )
br_transverse_mass_lep1_lep2 = array( 'd', [0] )
br_transverse_mass_lep1_met  = array( 'd', [0] )
br_transverse_mass_lep2_met  = array( 'd', [0] )
br_dphi_lep1_met = array( 'd', [0] )
br_dphi_lep2_met = array( 'd', [0] )
br_dphi_lep_lep  = array( 'd', [0] )
br_deta_lep_lep  = array( 'd', [0] )
br_dR_lep_lep    = array( 'd', [0] )
br_ptsum_lep1_lep2_met = array( 'd', [0] )
br_ptsum_lep1_lep2     = array( 'd', [0] )
br_pttot_lep1_lep2_met = array( 'd', [0] )
br_pttot_lep1_lep2     = array( 'd', [0] )
br_ptdiff_lep1_lep2    = array( 'd', [0] )

br_met_et_sm   = array( 'd', [0] )
br_lep1_pt_sm  = array( 'd', [0] )
br_lep1_eta_sm = array( 'd', [0] )
br_lep2_pt_sm  = array( 'd', [0] )
br_lep2_eta_sm = array( 'd', [0] )
br_transverse_mass_lep1_lep2_sm = array( 'd', [0] )
br_transverse_mass_lep1_met_sm  = array( 'd', [0] )
br_transverse_mass_lep2_met_sm  = array( 'd', [0] )
br_dphi_lep1_met_sm = array( 'd', [0] )
br_dphi_lep2_met_sm = array( 'd', [0] )
br_dphi_lep_lep_sm  = array( 'd', [0] )
br_deta_lep_lep_sm  = array( 'd', [0] )
br_dR_lep_lep_sm    = array( 'd', [0] )
br_ptsum_lep1_lep2_met_sm = array( 'd', [0] )
br_ptsum_lep1_lep2_sm     = array( 'd', [0] )
br_pttot_lep1_lep2_met_sm = array( 'd', [0] )
br_pttot_lep1_lep2_sm     = array( 'd', [0] )
br_ptdiff_lep1_lep2_sm    = array( 'd', [0] )

#----------------------------------------

otree.Branch('runnumber', br_runnumber, 'runnumber/I')
otree.Branch('evtnumber', br_evtnumber, 'evtnumber/I')

otree.Branch( 'weight', br_weight, 'weight/D' )
otree.Branch( 'bad_event', br_bad_event, 'bad_event/I' )

otree.Branch( 'n_photons', br_n_photons, 'n_photons/I' )

otree.Branch( 'boson_pt',    br_boson_pt,    'boson_pt/D'    )
otree.Branch( 'boson_eta',   br_boson_eta,   'boson_eta/D'   )
otree.Branch( 'boson_phi',   br_boson_phi,   'boson_phi/D'   )
otree.Branch( 'boson_E',     br_boson_E,     'boson_E/D'     )
otree.Branch( 'boson_m',     br_boson_m,     'boson_m/D'     )
otree.Branch( 'boson_mtrue', br_boson_mtrue, 'boson_mtrue/D' )

otree.Branch( 'taup_pt',    br_taup_pt,    'taup_pt/D'    )
otree.Branch( 'taup_eta',   br_taup_eta,   'taup_eta/D'   )
otree.Branch( 'taup_phi',   br_taup_phi,   'taup_phi/D'   )
otree.Branch( 'taup_E',     br_taup_E,     'taup_E/D'     )
otree.Branch( 'taup_m',     br_taup_m,     'taup_m/D'     )
otree.Branch( 'taup_mtrue', br_taup_mtrue, 'taup_mtrue/D' )

otree.Branch( 'taum_pt',    br_taum_pt,    'taum_pt/D'    )
otree.Branch( 'taum_eta',   br_taum_eta,   'taum_eta/D'   )
otree.Branch( 'taum_phi',   br_taum_phi,   'taum_phi/D'   )
otree.Branch( 'taum_E',     br_taum_E,     'taum_E/D'     )
otree.Branch( 'taum_m',     br_taum_m,     'taum_m/D'     )
otree.Branch( 'taum_mtrue', br_taum_mtrue, 'taum_mtrue/D' )

otree.Branch( 'ditau_pt',    br_ditau_pt,    'ditau_pt/D'    )
otree.Branch( 'ditau_eta',   br_ditau_eta,   'ditau_eta/D'   )
otree.Branch( 'ditau_phi',   br_ditau_phi,   'ditau_phi/D'   )
otree.Branch( 'ditau_E',     br_ditau_E,     'ditau_E/D'     )
otree.Branch( 'ditau_m',     br_ditau_m,     'ditau_m/D'     )
otree.Branch( 'ditau_mtrue', br_ditau_mtrue, 'ditau_mtrue/D' )

otree.Branch( 'tau1_decay_nutau_pt',    br_tau1_decay_nutau_pt,    'tau1_decay_nutau_pt/D'    )
otree.Branch( 'tau1_decay_nutau_eta',   br_tau1_decay_nutau_eta,   'tau1_decay_nutau_eta/D'   )
otree.Branch( 'tau1_decay_nutau_phi',   br_tau1_decay_nutau_phi,   'tau1_decay_nutau_phi/D'   )
otree.Branch( 'tau1_decay_nutau_E',     br_tau1_decay_nutau_E,     'tau1_decay_nutau_E/D'     )
otree.Branch( 'tau1_decay_nutau_m',     br_tau1_decay_nutau_m,     'tau1_decay_nutau_m/D'     )
otree.Branch( 'tau1_decay_nutau_mtrue', br_tau1_decay_nutau_mtrue, 'tau1_decay_nutau_mtrue/D' )

otree.Branch( 'tau1_decay_lep_pt',      br_tau1_decay_lep_pt,      'tau1_decay_lep_pt/D'      )
otree.Branch( 'tau1_decay_lep_eta',     br_tau1_decay_lep_eta,     'tau1_decay_lep_eta/D'     )
otree.Branch( 'tau1_decay_lep_phi',     br_tau1_decay_lep_phi,     'tau1_decay_lep_phi/D'     )
otree.Branch( 'tau1_decay_lep_E',       br_tau1_decay_lep_E,       'tau1_decay_lep_E/D'       )
otree.Branch( 'tau1_decay_lep_m',       br_tau1_decay_lep_m,       'tau1_decay_lep_m/D'       )
otree.Branch( 'tau1_decay_lep_mtrue',   br_tau1_decay_lep_mtrue,   'tau1_decay_lep_mtrue/D'   )

otree.Branch( 'tau1_decay_nulep_pt',    br_tau1_decay_nulep_pt,    'tau1_decay_nulep_pt/D'    )
otree.Branch( 'tau1_decay_nulep_eta',   br_tau1_decay_nulep_eta,   'tau1_decay_nulep_eta/D'   )
otree.Branch( 'tau1_decay_nulep_phi',   br_tau1_decay_nulep_phi,   'tau1_decay_nulep_phi/D'   )
otree.Branch( 'tau1_decay_nulep_E',     br_tau1_decay_nulep_E,     'tau1_decay_nulep_E/D'     )
otree.Branch( 'tau1_decay_nulep_m',     br_tau1_decay_nulep_m,     'tau1_decay_nulep_m/D'     )
otree.Branch( 'tau1_decay_nulep_mtrue', br_tau1_decay_nulep_mtrue, 'tau1_decay_nulep_mtrue/D' )

otree.Branch( 'tau2_decay_nutau_pt',    br_tau2_decay_nutau_pt,    'tau2_decay_nutau_pt/D'    )
otree.Branch( 'tau2_decay_nutau_eta',   br_tau2_decay_nutau_eta,   'tau2_decay_nutau_eta/D'   )
otree.Branch( 'tau2_decay_nutau_phi',   br_tau2_decay_nutau_phi,   'tau2_decay_nutau_phi/D'   )
otree.Branch( 'tau2_decay_nutau_E',     br_tau2_decay_nutau_E,     'tau2_decay_nutau_E/D'     )
otree.Branch( 'tau2_decay_nutau_m',     br_tau2_decay_nutau_m,     'tau2_decay_nutau_m/D'     )
otree.Branch( 'tau2_decay_nutau_mtrue', br_tau2_decay_nutau_mtrue, 'tau2_decay_nutau_mtrue/D' )

otree.Branch( 'tau2_decay_lep_pt',      br_tau2_decay_lep_pt,      'tau2_decay_lep_pt/D'      )
otree.Branch( 'tau2_decay_lep_eta',     br_tau2_decay_lep_eta,     'tau2_decay_lep_eta/D'     )
otree.Branch( 'tau2_decay_lep_phi',     br_tau2_decay_lep_phi,     'tau2_decay_lep_phi/D'     )
otree.Branch( 'tau2_decay_lep_E',       br_tau2_decay_lep_E,       'tau2_decay_lep_E/D'       )
otree.Branch( 'tau2_decay_lep_m',       br_tau2_decay_lep_m,       'tau2_decay_lep_m/D'       )
otree.Branch( 'tau2_decay_lep_mtrue',   br_tau2_decay_lep_mtrue,   'tau2_decay_lep_mtrue/D'   )

otree.Branch( 'tau2_decay_nulep_pt',    br_tau2_decay_nulep_pt,    'tau2_decay_nulep_pt/D'    )
otree.Branch( 'tau2_decay_nulep_eta',   br_tau2_decay_nulep_eta,   'tau2_decay_nulep_eta/D'   )
otree.Branch( 'tau2_decay_nulep_phi',   br_tau2_decay_nulep_phi,   'tau2_decay_nulep_phi/D'   )
otree.Branch( 'tau2_decay_nulep_E',     br_tau2_decay_nulep_E,     'tau2_decay_nulep_E/D'     )
otree.Branch( 'tau2_decay_nulep_m',     br_tau2_decay_nulep_m,     'tau2_decay_nulep_m/D'     )
otree.Branch( 'tau2_decay_nulep_mtrue', br_tau2_decay_nulep_mtrue, 'tau2_decay_nulep_mtrue/D' )

otree.Branch( 'leplep', br_leplep, 'leplep/I' )
otree.Branch( 'lephad', br_lephad, 'lephad/I' )
otree.Branch( 'hadhad', br_hadhad, 'hadhad/I' )

otree.Branch( 'lep1_pt',  br_lep1_pt,  'lep1_pt/D'  )
otree.Branch( 'lep1_eta', br_lep1_eta, 'lep1_eta/D' )
otree.Branch( 'lep2_pt',  br_lep2_pt,  'lep2_pt/D'  )
otree.Branch( 'lep2_eta', br_lep2_eta, 'lep2_eta/D' )
otree.Branch( 'met_et',   br_met_et,   'met_et/D'   )
otree.Branch( 'transverse_mass_lep1_lep2', br_transverse_mass_lep1_lep2, 'transverse_mass_lep1_lep2/D' )
otree.Branch( 'transverse_mass_lep1_met',  br_transverse_mass_lep1_met,  'transverse_mass_lep1_met/D'  )
otree.Branch( 'transverse_mass_lep2_met',  br_transverse_mass_lep2_met,  'transverse_mass_lep2_met/D'  )
otree.Branch( 'dphi_lep1_met', br_dphi_lep1_met, 'dphi_lep1_met/D' )
otree.Branch( 'dphi_lep2_met', br_dphi_lep2_met, 'dphi_lep2_met/D' )
otree.Branch( 'dphi_lep_lep',  br_dphi_lep_lep,  'dphi_lep_lep/D'  )
otree.Branch( 'deta_lep_lep',  br_deta_lep_lep,  'deta_lep_lep/D'  )
otree.Branch( 'dR_lep_lep',    br_dR_lep_lep,    'dR_lep_lep/D'    )
otree.Branch( 'ptsum_lep1_lep2_met', br_ptsum_lep1_lep2_met, 'ptsum_lep1_lep2_met/D' )
otree.Branch( 'ptsum_lep1_lep2',     br_ptsum_lep1_lep2,     'ptsum_lep1_lep2/D'     )
otree.Branch( 'pttot_lep1_lep2_met', br_pttot_lep1_lep2_met, 'pttot_lep1_lep2_met/D' )
otree.Branch( 'pttot_lep1_lep2',     br_pttot_lep1_lep2,     'pttot_lep1_lep2/D'     )
otree.Branch( 'ptdiff_lep1_lep2',    br_ptdiff_lep1_lep2,    'ptdiff_lep1_lep2/D'    )

otree.Branch( 'lep1_pt_sm',  br_lep1_pt_sm,  'lep1_pt_sm/D'  )
otree.Branch( 'lep1_eta_sm', br_lep1_eta_sm, 'lep1_eta_sm/D' )
otree.Branch( 'lep2_pt_sm',  br_lep2_pt_sm,  'lep2_pt_sm/D'  )
otree.Branch( 'lep2_eta_sm', br_lep2_eta_sm, 'lep2_eta_sm/D' )
otree.Branch( 'met_et_sm',   br_met_et_sm,   'met_et_sm/D'   )
otree.Branch( 'transverse_mass_lep1_lep2_sm', br_transverse_mass_lep1_lep2_sm, 'transverse_mass_lep1_lep2_sm/D' )
otree.Branch( 'transverse_mass_lep1_met_sm',  br_transverse_mass_lep1_met_sm,  'transverse_mass_lep1_met_sm/D'  )
otree.Branch( 'transverse_mass_lep2_met_sm',  br_transverse_mass_lep2_met_sm,  'transverse_mass_lep2_met_sm/D'  )
otree.Branch( 'dphi_lep1_met_sm', br_dphi_lep1_met_sm, 'dphi_lep1_met_sm/D' )
otree.Branch( 'dphi_lep2_met_sm', br_dphi_lep2_met_sm, 'dphi_lep2_met_sm/D' )
otree.Branch( 'dphi_lep_lep_sm',  br_dphi_lep_lep_sm,  'dphi_lep_lep_sm/D'  )
otree.Branch( 'deta_lep_lep_sm',  br_deta_lep_lep_sm,  'deta_lep_lep_sm/D'  )
otree.Branch( 'dR_lep_lep_sm',    br_dR_lep_lep_sm,    'dR_lep_lep_sm/D'    )
otree.Branch( 'ptsum_lep1_lep2_met_sm', br_ptsum_lep1_lep2_met_sm, 'ptsum_lep1_lep2_met_sm/D' )
otree.Branch( 'ptsum_lep1_lep2_sm',     br_ptsum_lep1_lep2_sm,     'ptsum_lep1_lep2_sm/D'     )
otree.Branch( 'pttot_lep1_lep2_met_sm', br_pttot_lep1_lep2_met_sm, 'pttot_lep1_lep2_met_sm/D' )
otree.Branch( 'pttot_lep1_lep2_sm',     br_pttot_lep1_lep2_sm,     'pttot_lep1_lep2_sm/D'     )
otree.Branch( 'ptdiff_lep1_lep2_sm',    br_ptdiff_lep1_lep2_sm,    'ptdiff_lep1_lep2_sm/D'    )

#otree.Branch( 'MT', br_mt, 'MT/D' )
#otree.Branch( 'Ptll', br_ptll, 'Ptll/D' )
#otree.Branch( 'DPhill', br_dphill, 'DPhill/D' )
#otree.Branch( 'Mll', br_mll, 'Mll/D' )
#otree.Branch( 'METRel', br_metrel, 'METRel/D' )
#otree.Branch( 'MET', br_met, 'MET/D' )
#otree.Branch( 'METTrue', br_mettrue, 'METTrue/D' )
#otree.Branch( 'x1', br_x1, 'x1/D' )
#otree.Branch( 'x2', br_x2, 'x2/D' )
#otree.Branch( 'fl1', br_fl1, 'fl1/D' )
#otree.Branch( 'fl2', br_fl2, 'fl2/D' )
#otree.Branch( 'qsqr', br_qsqr, 'qsqr/D' )
#otree.Branch( 'wplusPt', br_wpluspt, 'wplusPt/D' )
#otree.Branch( 'wminusPt', br_wminuspt, 'wminusPt/D' )
#otree.Branch( 'wplusEta', br_wpluseta, 'wplusEta/D' )
#otree.Branch( 'wplusM', br_wplusm, 'wplusM/D' )
#otree.Branch( 'wminusEta', br_wminuseta, 'wminusEta/D' )
#otree.Branch( 'wplusPhi', br_wplusphi, 'wplusPhi/D' )
#otree.Branch( 'wminusPhi', br_wminusphi, 'wminusPhi/D' )
#otree.Branch( 'wminusM', br_wminusm, 'wminusM/D' )

#----------------------------------------
#----------------------------------------

th1_boson_pt   = ROOT.TH1F("boson_pt",  "boson_pt",  1000,0.0,1000.0)
th1_boson_pt_w = ROOT.TH1F("boson_pt_w","boson_pt_w",1000,0.0,1000.0)
th1_ditau_pt   = ROOT.TH1F("ditau_pt",  "ditau_pt",  1000,0.0,1000.0)
th1_ditau_pt_w = ROOT.TH1F("ditau_pt_w","ditau_pt_w",1000,0.0,1000.0)

tlv_boson = ROOT.TLorentzVector()
tlv_taup  = ROOT.TLorentzVector()
tlv_taum  = ROOT.TLorentzVector()
tlv_taup_decay_nutau  = ROOT.TLorentzVector()
tlv_taup_decay_lep    = ROOT.TLorentzVector()
tlv_taup_decay_nulep  = ROOT.TLorentzVector()
tlv_taum_decay_nutau  = ROOT.TLorentzVector()
tlv_taum_decay_lep    = ROOT.TLorentzVector()
tlv_taum_decay_nulep  = ROOT.TLorentzVector()

tlv_lep1 = ROOT.TLorentzVector()
tlv_lep2 = ROOT.TLorentzVector()
tlv_met  = ROOT.TLorentzVector()

#----------------------------------------

for ievt in xrange( nevts ):
    if int(float(ievt) / float(nevts)) % 10 != 0 or ievt == 0:
        print 'Processed %s / %s (%s)' % (ievt, nevts, ievt / nevts)

    # print 100 * "--"
    itree.GetEntry( ievt )

    br_weight[0] = itree.mcevt_weight[0][0]
    br_runnumber[0] = itree.RunNumber
    br_evtnumber[0] = itree.EventNumber
    print itree.RunNumber, itree.EventNumber
    ############################
    ## Get the Higgs
    ############################

    Nb = 0
    ib = -1
    for imc in xrange( itree.mc_n ):
        # ugly but seems to work :-(
        if itree.mc_pdgId[imc] == 25 and len(itree.mc_child_index[imc])>1:
            Nb += 1
            ib = imc
    
    # if Nb != 1:
    #     sys.stdout.write( '********* Found '+str(Nb)+' bosons !!!\n' )

    if ib < 0 or Nb != 1: 
        continue
    

    tlv_boson.SetPtEtaPhiE( itree.mc_pt[ib], itree.mc_eta[ib], itree.mc_phi[ib], itree.mc_E[ib] )

    br_boson_pt[0]    = tlv_boson.Pt()
    br_boson_eta[0]   = tlv_boson.Eta()
    br_boson_phi[0]   = tlv_boson.Phi()
    br_boson_E[0]     = tlv_boson.E()
    br_boson_m[0]     = tlv_boson.M()
    br_boson_mtrue[0] = itree.mc_m[ib]

    th1_boson_pt.Fill(itree.mc_pt[ib]/1000.)
    th1_boson_pt_w.Fill(itree.mc_pt[ib]/1000.,itree.mcevt_weight[0][0])
    
    ######################################
    ## Get the taus from the Higgs decay
    ######################################

    itp = -1
    itm = -1
    Nph = 0

    # for ich in itree.mc_child_index[ib]:
    #     sys.stdout.write( 'H boson child index '+str(ich)+' has pdgId '+str(itree.mc_pdgId[ich])+' and status '+str(itree.mc_status[ich])+'\n' )

    for ich in itree.mc_child_index[ib]:
        if abs(itree.mc_pdgId[ich]) == 15 and itree.mc_status[ich] == 2:
            if itree.mc_pdgId[ich] > 0: 
                itm = ich
            else:
                itp = ich
        #Check for FSR/Loop correction emitted photons
        elif itree.mc_pdgId[ich]==22 and itree.mc_status[ich]==1:
            Nph+=1
    br_n_photons[0]=Nph

    if itp < 0 or itm < 0: continue

    tlv_taup.SetPtEtaPhiE( itree.mc_pt[itp], itree.mc_eta[itp], itree.mc_phi[itp], itree.mc_E[itp] )
    tlv_taum.SetPtEtaPhiE( itree.mc_pt[itm], itree.mc_eta[itm], itree.mc_phi[itm], itree.mc_E[itm] )

    br_taup_pt[0]    = tlv_taup.Pt()
    br_taup_eta[0]   = tlv_taup.Eta()
    br_taup_phi[0]   = tlv_taup.Phi()
    br_taup_E[0]     = tlv_taup.E()
    br_taup_m[0]     = tlv_taup.M()
    br_taup_mtrue[0] = itree.mc_m[itp]

    br_taum_pt[0]    = tlv_taum.Pt()
    br_taum_eta[0]   = tlv_taum.Eta()
    br_taum_phi[0]   = tlv_taum.Phi()
    br_taum_E[0]     = tlv_taum.E()
    br_taum_m[0]     = tlv_taum.M()
    br_taum_mtrue[0] = itree.mc_m[itm]

    br_ditau_pt[0]    = ( tlv_taup + tlv_taum ).Pt()
    br_ditau_eta[0]   = ( tlv_taup + tlv_taum ).Eta()
    br_ditau_phi[0]   = ( tlv_taup + tlv_taum ).Phi()
    br_ditau_E[0]     = ( tlv_taup + tlv_taum ).E()
    br_ditau_m[0]     = ( tlv_taup + tlv_taum ).M()
    br_ditau_mtrue[0] = -9999

    th1_ditau_pt.Fill((tlv_taup + tlv_taum).Pt()/1000.)
    th1_ditau_pt_w.Fill((tlv_taup + tlv_taum).Pt()/1000.,itree.mcevt_weight[0][0])
        
    ################################################################
    ## Get the decay products of the taus that decay leptonically.
    ################################################################

    taup_lep = 0
    taup_had = 0
    tauptaum_leplep = 0
    tauptaum_lephad = 0
    tauptaum_hadlep = 0
    tauptaum_hadhad = 0

    ipnutau = -1
    iplep   = -1
    ipnulep = -1
    iw     = -1
    for ich in itree.mc_child_index[itp]:
        if itree.mc_pdgId[ich] in [-11,12,-13,14,-16] and itree.mc_status[ich] == 1:
            if ipnutau < 0 and itree.mc_pdgId[ich] == -16: ipnutau = ich
            if iplep   < 0 and itree.mc_pdgId[ich] in [-11,-13]: iplep   = ich
            if ipnulep < 0 and itree.mc_pdgId[ich] in [ 12, 14]: ipnulep = ich
        if itree.mc_pdgId[ich] == 24 and itree.mc_status[ich] == 2:
            if iw < 0: iw = ich
    if iw > 0 and ipnutau > 0 and iplep < 0 and ipnulep < 0:
        for ich in itree.mc_child_index[iw]:
            if itree.mc_pdgId[ich] in [-11,12,-13,14] and itree.mc_status[ich] == 1:
                if iplep   < 0 and itree.mc_pdgId[ich] in [-11,-13]: iplep   = ich
                if ipnulep < 0 and itree.mc_pdgId[ich] in [ 12, 14]: ipnulep = ich
    if ipnutau > 0:
        tlv_taup_decay_nutau.SetPtEtaPhiE( itree.mc_pt[ipnutau], itree.mc_eta[ipnutau], itree.mc_phi[ipnutau], itree.mc_E[ipnutau] )
        if iplep > 0 and ipnulep > 0:
            tlv_taup_decay_lep.SetPtEtaPhiE( itree.mc_pt[iplep], itree.mc_eta[iplep], itree.mc_phi[iplep], itree.mc_E[iplep] )
            tlv_taup_decay_nulep.SetPtEtaPhiE( itree.mc_pt[ipnulep], itree.mc_eta[ipnulep], itree.mc_phi[ipnulep], itree.mc_E[ipnulep] )
            taup_lep = 1
        else:
            tlv_taup_decay_lep = tlv_taup - tlv_taup_decay_nutau
            taup_had = 1

    if not (taup_lep or taup_had): continue

    imnutau = -1
    imlep   = -1
    imnulep = -1
    iw     = -1
    for ich in itree.mc_child_index[itm]:
        if itree.mc_pdgId[ich] in [11,-12,13,-14,16] and itree.mc_status[ich] == 1:
            if imnutau < 0 and itree.mc_pdgId[ich] == 16: imnutau = ich
            if imlep   < 0 and itree.mc_pdgId[ich] in [ 11, 13]: imlep   = ich
            if imnulep < 0 and itree.mc_pdgId[ich] in [-12,-14]: imnulep = ich
        if itree.mc_pdgId[ich] == -24 and itree.mc_status[ich] == 2:
            if iw < 0: iw = ich
    if iw > 0 and imnutau > 0 and imlep < 0 and imnulep < 0:
        for ich in itree.mc_child_index[iw]:
            if itree.mc_pdgId[ich] in [11,-12,13,-14] and itree.mc_status[ich] == 1:
                if imlep   < 0 and itree.mc_pdgId[ich] in [ 11, 13]: imlep   = ich
                if imnulep < 0 and itree.mc_pdgId[ich] in [-12,-14]: imnulep = ich
    if imnutau > 0:
        tlv_taum_decay_nutau.SetPtEtaPhiE( itree.mc_pt[imnutau], itree.mc_eta[imnutau], itree.mc_phi[imnutau], itree.mc_E[imnutau] )
        if imlep > 0 and imnulep > 0:
            tlv_taum_decay_lep.SetPtEtaPhiE( itree.mc_pt[imlep], itree.mc_eta[imlep], itree.mc_phi[imlep], itree.mc_E[imlep] )
            tlv_taum_decay_nulep.SetPtEtaPhiE( itree.mc_pt[imnulep], itree.mc_eta[imnulep], itree.mc_phi[imnulep], itree.mc_E[imnulep] )
            if taup_lep == 1:
                tauptaum_leplep = 1
            if taup_had == 1:
                tauptaum_hadlep = 1
        else:
            tlv_taum_decay_lep = tlv_taum - tlv_taum_decay_nutau
            if taup_lep == 1:
                tauptaum_lephad = 1
            if taup_had == 1:
                tauptaum_hadhad = 1
    else:
        continue

    if not (tauptaum_leplep or tauptaum_lephad or tauptaum_hadlep or tauptaum_hadhad): continue

    if ipnutau > 0:
        br_tau1_decay_nutau_pt[0]    = tlv_taup_decay_nutau.Pt()
        br_tau1_decay_nutau_eta[0]   = tlv_taup_decay_nutau.Eta()
        br_tau1_decay_nutau_phi[0]   = tlv_taup_decay_nutau.Phi()
        br_tau1_decay_nutau_E[0]     = tlv_taup_decay_nutau.E()
        br_tau1_decay_nutau_m[0]     = tlv_taup_decay_nutau.M()
        br_tau1_decay_nutau_mtrue[0] = itree.mc_m[ipnutau]
    if imnutau > 0:
        br_tau2_decay_nutau_pt[0]    = tlv_taum_decay_nutau.Pt()
        br_tau2_decay_nutau_eta[0]   = tlv_taum_decay_nutau.Eta()
        br_tau2_decay_nutau_phi[0]   = tlv_taum_decay_nutau.Phi()
        br_tau2_decay_nutau_E[0]     = tlv_taum_decay_nutau.E()
        br_tau2_decay_nutau_m[0]     = tlv_taum_decay_nutau.M()
        br_tau2_decay_nutau_mtrue[0] = itree.mc_m[imnutau]
    if tauptaum_leplep > 0:
        br_tau1_decay_lep_pt[0]      = tlv_taup_decay_lep.Pt()
        br_tau1_decay_lep_eta[0]     = tlv_taup_decay_lep.Eta()
        br_tau1_decay_lep_phi[0]     = tlv_taup_decay_lep.Phi()
        br_tau1_decay_lep_E[0]       = tlv_taup_decay_lep.E()
        br_tau1_decay_lep_m[0]       = tlv_taup_decay_lep.M()
        br_tau1_decay_lep_mtrue[0]   = itree.mc_m[iplep]
        br_tau1_decay_nulep_pt[0]    = tlv_taup_decay_nulep.Pt()
        br_tau1_decay_nulep_eta[0]   = tlv_taup_decay_nulep.Eta()
        br_tau1_decay_nulep_phi[0]   = tlv_taup_decay_nulep.Phi()
        br_tau1_decay_nulep_E[0]     = tlv_taup_decay_nulep.E()
        br_tau1_decay_nulep_m[0]     = tlv_taup_decay_nulep.M()
        br_tau1_decay_nulep_mtrue[0] = itree.mc_m[ipnulep]
        br_tau2_decay_lep_pt[0]      = tlv_taum_decay_lep.Pt()
        br_tau2_decay_lep_eta[0]     = tlv_taum_decay_lep.Eta()
        br_tau2_decay_lep_phi[0]     = tlv_taum_decay_lep.Phi()
        br_tau2_decay_lep_E[0]       = tlv_taum_decay_lep.E()
        br_tau2_decay_lep_m[0]       = tlv_taum_decay_lep.M()
        br_tau2_decay_lep_mtrue[0]   = itree.mc_m[imlep]
        br_tau2_decay_nulep_pt[0]    = tlv_taum_decay_nulep.Pt()
        br_tau2_decay_nulep_eta[0]   = tlv_taum_decay_nulep.Eta()
        br_tau2_decay_nulep_phi[0]   = tlv_taum_decay_nulep.Phi()
        br_tau2_decay_nulep_E[0]     = tlv_taum_decay_nulep.E()
        br_tau2_decay_nulep_m[0]     = tlv_taum_decay_nulep.M()
        br_tau2_decay_nulep_mtrue[0] = itree.mc_m[imnulep]
    if tauptaum_hadhad > 0:
        br_tau1_decay_lep_pt[0]      = tlv_taup_decay_lep.Pt()
        br_tau1_decay_lep_eta[0]     = tlv_taup_decay_lep.Eta()
        br_tau1_decay_lep_phi[0]     = tlv_taup_decay_lep.Phi()
        br_tau1_decay_lep_E[0]       = tlv_taup_decay_lep.E()
        br_tau1_decay_lep_m[0]       = tlv_taup_decay_lep.M()
        br_tau1_decay_lep_mtrue[0]   = itree.mc_m[itp]
        br_tau2_decay_lep_pt[0]      = tlv_taum_decay_lep.Pt()
        br_tau2_decay_lep_eta[0]     = tlv_taum_decay_lep.Eta()
        br_tau2_decay_lep_phi[0]     = tlv_taum_decay_lep.Phi()
        br_tau2_decay_lep_E[0]       = tlv_taum_decay_lep.E()
        br_tau2_decay_lep_m[0]       = tlv_taum_decay_lep.M()
        br_tau2_decay_lep_mtrue[0]   = itree.mc_m[itm]
    if tauptaum_lephad > 0:
        br_tau1_decay_lep_pt[0]      = tlv_taup_decay_lep.Pt()
        br_tau1_decay_lep_eta[0]     = tlv_taup_decay_lep.Eta()
        br_tau1_decay_lep_phi[0]     = tlv_taup_decay_lep.Phi()
        br_tau1_decay_lep_E[0]       = tlv_taup_decay_lep.E()
        br_tau1_decay_lep_m[0]       = tlv_taup_decay_lep.M()
        br_tau1_decay_lep_mtrue[0]   = itree.mc_m[iplep]
        br_tau1_decay_nulep_pt[0]    = tlv_taup_decay_nulep.Pt()
        br_tau1_decay_nulep_eta[0]   = tlv_taup_decay_nulep.Eta()
        br_tau1_decay_nulep_phi[0]   = tlv_taup_decay_nulep.Phi()
        br_tau1_decay_nulep_E[0]     = tlv_taup_decay_nulep.E()
        br_tau1_decay_nulep_m[0]     = tlv_taup_decay_nulep.M()
        br_tau1_decay_nulep_mtrue[0] = itree.mc_m[ipnulep]
        br_tau2_decay_lep_pt[0]      = tlv_taum_decay_lep.Pt()
        br_tau2_decay_lep_eta[0]     = tlv_taum_decay_lep.Eta()
        br_tau2_decay_lep_phi[0]     = tlv_taum_decay_lep.Phi()
        br_tau2_decay_lep_E[0]       = tlv_taum_decay_lep.E()
        br_tau2_decay_lep_m[0]       = tlv_taum_decay_lep.M()
        br_tau2_decay_lep_mtrue[0]   = itree.mc_m[itm]
    if tauptaum_hadlep > 0:
        br_tau1_decay_lep_pt[0]      = tlv_taum_decay_lep.Pt()
        br_tau1_decay_lep_eta[0]     = tlv_taum_decay_lep.Eta()
        br_tau1_decay_lep_phi[0]     = tlv_taum_decay_lep.Phi()
        br_tau1_decay_lep_E[0]       = tlv_taum_decay_lep.E()
        br_tau1_decay_lep_m[0]       = tlv_taum_decay_lep.M()
        br_tau1_decay_lep_mtrue[0]   = itree.mc_m[imlep]
        br_tau1_decay_nulep_pt[0]    = tlv_taum_decay_nulep.Pt()
        br_tau1_decay_nulep_eta[0]   = tlv_taum_decay_nulep.Eta()
        br_tau1_decay_nulep_phi[0]   = tlv_taum_decay_nulep.Phi()
        br_tau1_decay_nulep_E[0]     = tlv_taum_decay_nulep.E()
        br_tau1_decay_nulep_m[0]     = tlv_taum_decay_nulep.M()
        br_tau1_decay_nulep_mtrue[0] = itree.mc_m[imnulep]
        br_tau2_decay_lep_pt[0]      = tlv_taup_decay_lep.Pt()
        br_tau2_decay_lep_eta[0]     = tlv_taup_decay_lep.Eta()
        br_tau2_decay_lep_phi[0]     = tlv_taup_decay_lep.Phi()
        br_tau2_decay_lep_E[0]       = tlv_taup_decay_lep.E()
        br_tau2_decay_lep_m[0]       = tlv_taup_decay_lep.M()
        br_tau2_decay_lep_mtrue[0]   = itree.mc_m[itp]

    tlv_met = tlv_taup_decay_nutau + tlv_taum_decay_nutau
    if tauptaum_leplep > 0 or tauptaum_lephad > 0: tlv_met += tlv_taup_decay_nulep
    if tauptaum_leplep > 0 or tauptaum_hadlep > 0: tlv_met += tlv_taum_decay_nulep

    leplep = tauptaum_leplep
    hadhad = tauptaum_hadhad
    if tauptaum_lephad > 0 or tauptaum_hadlep > 0: lephad = 1
    else:                                          lephad = 0

    br_leplep[0] = leplep
    br_lephad[0] = lephad
    br_hadhad[0] = hadhad

    ##################################################################################################
    ## Get the leading and subleading lepton, where by lepton I mean electron, muon or hadronic-tau.
    ##################################################################################################

    if tlv_taup_decay_lep.Pt() > tlv_taum_decay_lep.Pt():
        tlv_lep1 = tlv_taup_decay_lep
        tlv_lep2 = tlv_taum_decay_lep
    else:
        tlv_lep1 = tlv_taum_decay_lep
        tlv_lep2 = tlv_taup_decay_lep

    if tauptaum_leplep > 0:
        lep1_label = 'emu'
        lep2_label = 'emu'
    elif tauptaum_lephad > 0:
        if tlv_taup_decay_lep.Pt() > tlv_taum_decay_lep.Pt():
            lep1_label = 'emu'
            lep2_label = 'tauhad'
        else:
            lep1_label = 'tauhad'
            lep2_label = 'emu'
    elif tauptaum_hadlep > 0:
        if tlv_taup_decay_lep.Pt() > tlv_taum_decay_lep.Pt():
            lep1_label = 'tauhad'
            lep2_label = 'emu'
        else:
            lep1_label = 'emu'
            lep2_label = 'tauhad'
    elif tauptaum_hadhad > 0:
        lep1_label = 'tauhad'
        lep2_label = 'tauhad'

    met_et   = tlv_met.Pt()
    met_phi  = tlv_met.Phi()
    lep1_pt  = tlv_lep1.Pt()
    lep1_eta = tlv_lep1.Eta()
    lep1_phi = tlv_lep1.Phi()
    lep1_E   = tlv_lep1.E()
    lep2_pt  = tlv_lep2.Pt()
    lep2_eta = tlv_lep2.Eta()
    lep2_phi = tlv_lep2.Phi()
    lep2_E   = tlv_lep2.E()
    dphi_lep1_lep2 = delta_phi(lep1_phi,lep2_phi)
    dphi_lep1_met  = delta_phi(lep1_phi,met_phi )
    dphi_lep2_met  = delta_phi(lep2_phi,met_phi )
    deta_lep1_lep2 = lep2_eta-lep1_eta

    br_met_et[0]  = met_et
    br_lep1_pt[0]  = lep1_pt
    br_lep1_eta[0] = lep1_eta
    br_lep2_pt[0]  = lep2_pt
    br_lep2_eta[0] = lep2_eta
    br_transverse_mass_lep1_lep2[0] = math.sqrt(2*lep1_pt*lep2_pt*(1-math.cos(dphi_lep1_lep2)))
    br_transverse_mass_lep1_met[0]  = math.sqrt(2*lep1_pt*met_et *(1-math.cos(dphi_lep1_met )))
    br_transverse_mass_lep2_met[0]  = math.sqrt(2*lep2_pt*met_et *(1-math.cos(dphi_lep2_met )))
    br_dphi_lep1_met[0] = abs(dphi_lep1_met)
    br_dphi_lep2_met[0] = abs(dphi_lep2_met)
    br_dphi_lep_lep[0]  = abs(dphi_lep1_lep2)
    br_deta_lep_lep[0]  = abs(deta_lep1_lep2)
    br_dR_lep_lep[0]    = math.sqrt(pow(dphi_lep1_lep2,2)+pow(deta_lep1_lep2,2))
    br_ptsum_lep1_lep2_met[0] = lep1_pt + lep2_pt + met_et
    br_ptsum_lep1_lep2[0]     = lep1_pt + lep2_pt
    br_pttot_lep1_lep2_met[0] = ( tlv_lep1 + tlv_lep2 + tlv_met ).Pt() / ( lep1_pt + lep2_pt + met_et )
    br_pttot_lep1_lep2[0]     = ( tlv_lep1 + tlv_lep2 ).Pt() / ( lep1_pt + lep2_pt )
    br_ptdiff_lep1_lep2[0]    = ( tlv_lep1 - tlv_lep2 ).Pt() / ( lep1_pt + lep2_pt )

    smf_met_et  = abs(random.gauss(1.0,met_energy_resolution(tlv_lep1.Pt()/1000.+tlv_lep2.Pt()/1000.)/(met_et/1000.)))
    smf_met_phi = abs(random.gauss(1.0,0.04))
    if   lep1_label == 'emu':
        smf_lep1_eta = abs(random.gauss(1.0,0.02))
        smf_lep1_phi = abs(random.gauss(1.0,0.01))
        smf_lep1_pt  = abs(random.gauss(1.0,electron_relative_energy_resolution(tlv_lep1.E()/1000.,tlv_lep1.Eta())))
    elif lep1_label == 'tauhad':
        smf_lep1_eta = abs(random.gauss(1.0,0.04))
        smf_lep1_phi = abs(random.gauss(1.0,0.02))
        smf_lep1_pt  = abs(random.gauss(1.0,tau_relative_momentum_resolution(tlv_lep1.Pt()/1000.,tlv_lep1.Eta())))
    if   lep2_label == 'emu':
        smf_lep2_eta = abs(random.gauss(1.0,0.02))
        smf_lep2_phi = abs(random.gauss(1.0,0.01))
        smf_lep2_pt  = abs(random.gauss(1.0,electron_relative_energy_resolution(tlv_lep2.E()/1000.,tlv_lep2.Eta())))
    elif lep2_label == 'tauhad':
        smf_lep2_eta = abs(random.gauss(1.0,0.04))
        smf_lep2_phi = abs(random.gauss(1.0,0.02))
        smf_lep2_pt  = abs(random.gauss(1.0,tau_relative_momentum_resolution(tlv_lep2.Pt()/1000.,tlv_lep2.Eta())))
    
    met_et_sm   = met_et   * smf_met_et
    met_phi_sm  = met_phi  * smf_met_phi
    while met_phi_sm >  math.pi: met_phi_sm -= 2*math.pi
    while met_phi_sm < -math.pi: met_phi_sm += 2*math.pi
    lep1_pt_sm  = lep1_pt  * smf_lep1_pt
    lep1_eta_sm = lep1_eta * smf_lep1_eta
    lep1_phi_sm = lep1_phi * smf_lep1_phi
    while lep1_phi_sm >  math.pi: lep1_phi_sm -= 2*math.pi
    while lep1_phi_sm < -math.pi: lep1_phi_sm += 2*math.pi
    lep2_pt_sm  = lep2_pt  * smf_lep2_pt
    lep2_eta_sm = lep2_eta * smf_lep2_eta
    lep2_phi_sm = lep2_phi * smf_lep2_phi
    while lep2_phi_sm >  math.pi: lep2_phi_sm -= 2*math.pi
    while lep2_phi_sm < -math.pi: lep2_phi_sm += 2*math.pi
    dphi_lep1_lep2_sm = delta_phi(lep1_phi_sm,lep2_phi_sm)
    dphi_lep1_met_sm  = delta_phi(lep1_phi_sm,met_phi_sm )
    dphi_lep2_met_sm  = delta_phi(lep2_phi_sm,met_phi_sm )
    deta_lep1_lep2_sm = lep2_eta_sm-lep1_eta_sm

    br_met_et_sm[0]   = met_et_sm
    br_lep1_pt_sm[0]  = lep1_pt_sm
    br_lep1_eta_sm[0] = lep1_eta_sm
    br_lep2_pt_sm[0]  = lep2_pt_sm
    br_lep2_eta_sm[0] = lep2_eta_sm
    br_transverse_mass_lep1_lep2_sm[0] = math.sqrt(2*lep1_pt_sm*lep2_pt_sm*(1-math.cos(dphi_lep1_lep2_sm)))
    br_transverse_mass_lep1_met_sm[0]  = math.sqrt(2*lep1_pt_sm*met_et_sm *(1-math.cos(dphi_lep1_met_sm )))
    br_transverse_mass_lep2_met_sm[0]  = math.sqrt(2*lep2_pt_sm*met_et_sm *(1-math.cos(dphi_lep2_met_sm )))
    br_dphi_lep1_met_sm[0] = abs(dphi_lep1_met_sm )
    br_dphi_lep2_met_sm[0] = abs(dphi_lep2_met_sm )
    br_dphi_lep_lep_sm[0]  = abs(dphi_lep1_lep2_sm)
    br_deta_lep_lep_sm[0]  = abs(deta_lep1_lep2_sm)
    br_dR_lep_lep_sm[0]    = math.sqrt(pow(dphi_lep1_lep2_sm,2)+pow(deta_lep1_lep2_sm,2))
    br_ptsum_lep1_lep2_met_sm[0] = lep1_pt_sm + lep2_pt_sm + met_et_sm
    br_ptsum_lep1_lep2_sm[0]     = lep1_pt_sm + lep2_pt_sm
    br_pttot_lep1_lep2_met_sm[0] = math.sqrt(pow(lep1_pt_sm*math.cos(lep1_phi_sm)+lep2_pt_sm*math.cos(lep2_phi_sm)+met_et_sm*math.cos(met_phi_sm),2)+pow(lep1_pt_sm*math.sin(lep1_phi_sm)+lep2_pt_sm*math.sin(lep2_phi_sm)+met_et_sm*math.sin(met_phi_sm),2)) / ( lep1_pt_sm + lep2_pt_sm + met_et_sm )
    br_pttot_lep1_lep2_sm[0]     = math.sqrt(pow(lep1_pt_sm*math.cos(lep1_phi_sm)+lep2_pt_sm*math.cos(lep2_phi_sm),2)+pow(lep1_pt_sm*math.sin(lep1_phi_sm)+lep2_pt_sm*math.sin(lep2_phi_sm),2)) / ( lep1_pt_sm + lep2_pt_sm )
    br_ptdiff_lep1_lep2_sm[0]    = math.sqrt(pow(lep1_pt_sm*math.cos(lep1_phi_sm)-lep2_pt_sm*math.cos(lep2_phi_sm),2)+pow(lep1_pt_sm*math.sin(lep1_phi_sm)-lep2_pt_sm*math.sin(lep2_phi_sm),2)) / ( lep1_pt_sm + lep2_pt_sm )

    #for imc in xrange( itree.mc_n ):        
    #    if itp<0 or itm<0:        
    #        if abs(itree.mc_pdgId[imc]) in [15] and 120 < itree.mc_status[imc] < 130:
    #            decayed = False
    #            jmc = imc
    #            while not decayed:
    #                if itree.mc_pdgId[jmc] > 0: ilm = imc
    #                else:                       ilp = imc
    #                decayed = True
    #                for kmc in itree.mc_child_index[jmc]:
    #                    decayed = False
    #                    if abs(itree.mc_pdgId[jmc]) in [11,13,15]:
    #                        if itree.mc_pdgId[jmc] > 0: ilm = jmc
    #                        else:                       ilp = jmc
    #                        jmc = kmc
    #                    
    #    if itp >= 0 and itm >= 0:
    #        break
    
    otree.Fill()

ofile.cd()
otree.Write()
th1_boson_pt.Write()
th1_boson_pt_w.Write()
th1_ditau_pt.Write()
th1_ditau_pt_w.Write()
ofile.Close()

#------------------------------------
