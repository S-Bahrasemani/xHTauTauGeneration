


vbf = """
#--------------------------------------------------------------
# Powheg VBF_H setup
#--------------------------------------------------------------
include('PowhegControl/PowhegControl_VBF_H_Common.py')

# Set Powheg variables, overriding defaults
PowhegConfig.mass_H  = {0}
PowhegConfig.width_H = {1}

# Increase number of events requested to compensate for filter efficiency
PowhegConfig.nEvents *= 15

# Generate Powheg events
PowhegConfig.generateRunCard()
PowhegConfig.generateEvents()

#--------------------------------------------------------------
# Pythia8 showering
#--------------------------------------------------------------
include('MC12JobOptions/PowhegPythia8_AU2_CT10_Common.py')
include('MC12JobOptions/Pythia8_Photos.py')
topAlg.Pythia8.Commands += [
                            '25:onMode = off',#decay of Higgs
                            '25:onIfMatch = 15 15'
                            ]

#-----------------------------------------------------------------------
# Filter the tau kinematics: 2 taus with p_T(l/h)>10/20 GeV and |eta(l/h)| <2.5 
# l: e, mu; h: hadronic tau
#-----------------------------------------------------------------------
include('MC12JobOptions/TauFilter.py')
topAlg.TauFilter.Ntaus = 2
topAlg.TauFilter.Ptcute = 12000.0 #MeV
topAlg.TauFilter.Ptcutmu = 12000.0 #MeV
topAlg.TauFilter.Ptcuthad = 20000.0 #MeV

#--------------------------------------------------------------
# EVGEN configuration
#--------------------------------------------------------------
evgenConfig.description = 'POWHEG+PYTHIA8 VBF H->tautau with AU2 CT10'
evgenConfig.keywords    = [ 'SMhiggs', 'VBF', 'tau']
evgenConfig.contact     = [ 'Junichi.Tanaka@cern.ch' ]
evgenConfig.generators += [ 'Powheg', 'Pythia8' ]
"""
########

ggh = """
#--------------------------------------------------------------
# Powheg ggH_quark_mass_effects setup
#--------------------------------------------------------------
include('PowhegControl/PowhegControl_ggH_quark_mass_effects_Common.py')
#--------------------------------------------------------------

# Set Powheg variables, overriding defaults
PowhegConfig.mass_H  = {0}
PowhegConfig.width_H = {1}

# Turn on the heavy quark effect
PowhegConfig.use_massive_b = True
PowhegConfig.use_massive_c = True

# Increase number of events requested to compensate for filter efficiency
PowhegConfig.nEvents *= 20

# Generate Powheg events
PowhegConfig.generateRunCard()
PowhegConfig.generateEvents()

#--------------------------------------------------------------
# Pythia8 showering
#--------------------------------------------------------------
include('MC12JobOptions/PowhegPythia8_AU2_CT10_Common.py')
include('MC12JobOptions/Pythia8_Photos.py')
topAlg.Pythia8.Commands += [
                            '25:onMode = off',#decay of Higgs
                            '25:onIfMatch = 15 15'
                            ]

#-----------------------------------------------------------------------
# Filter the tau kinematics: 2 taus with p_T(l/h)>20 GeV and |eta(l/h)| <2.5 
# l: e, mu; h: hadronic tau
#-----------------------------------------------------------------------
include('MC12JobOptions/TauFilter.py')
topAlg.TauFilter.Ntaus = 2
topAlg.TauFilter.Ptcute = 12000.0 #MeV
topAlg.TauFilter.Ptcutmu = 12000.0 #MeV
topAlg.TauFilter.Ptcuthad = 20000.0 #MeV

#--------------------------------------------------------------
# EVGEN configuration
#--------------------------------------------------------------
evgenConfig.description = 'POWHEG+PYTHIA8 ggH->tautau with AU2 CT10'
evgenConfig.keywords    = [ 'SMhiggs', 'ggH', 'tau']
evgenConfig.contact     = [ 'Quentin.Buat@cern.ch' ]
evgenConfig.generators += [ 'Powheg', 'Pythia8' ]
"""



from yellowhiggs import br
from yellowhiggs.xsbr import __BR as BRs
masses_yh = BRs['2boson'].keys()

def get_jo(mass=125, mode='VBF'):
    """
    Determine width from yellowhiggs report
    and output a jo file for a given m_H
    """
    best_mass = sorted(masses_yh, key=lambda m: abs(m - float(mass)))[0]
    width = br(best_mass, '2boson')[0]
    print '%s: Mass = %s; Width = %s (calculated with m_H = %s)' % (mode, mass, width, best_mass)
    if mode == 'VBF':
        return vbf.format(mass, width)
    elif mode == 'gg':
        return ggh.format(mass, width)
    else:
        raise RuntimeError('Wrong mode !')
        
