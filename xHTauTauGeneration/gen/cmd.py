import re
import os
JO_PATTERN = re.compile('(MC12.(?P<id>\d+).PowhegPythia8_AU2CT10_(?P<mode>gg|VBF)H(?P<mass>\d+)_tautau.py)?$')
DIR_PATTERN = re.compile('(MC12.(?P<id>\d+).PowhegPythia8_AU2CT10_(?P<mode>gg|VBF)H(?P<mass>\d+)_tautau_run(?P<seed>\d+))?$')

def get_dir_mode(d):
    match = re.match(DIR_PATTERN, d)
    if match:
        return match.group('mode')
    else:
        return None

def get_dir_mass(d):
    match = re.match(DIR_PATTERN, d)
    if match:
        return int(match.group('mass'))
    else:
        return None

def get_dir_seed(d):
    match = re.match(DIR_PATTERN, d)
    if match:
        return int(match.group('seed'))
    else:
        return None


def generate_cmd(jo, seed=1234, nevents=5000):
    """
    EVNT GENERATION COMMAND LINE
    """
    match = re.match(JO_PATTERN, jo)
    if match:
        datasetid = match.group('id')
        run_dir = jo.replace('.py', '_run{0}'.format(seed))
        
        cmd_args = [
            "Generate_trf.py",
            "ecmEnergy=13000",
            "runNumber={0}".format(datasetid),
            "firstEvent=0",
            "maxEvents={0}".format(nevents), 
            "randomSeed={0}".format(seed), 
            "jobConfig={0}".format(jo),
            "outputEVNTFile=EVNT_{0}_seed{1}.root".format(datasetid, seed),
            "evgenJobOpts=MC12JobOpts-00-14-70_v8.tar.gz"
            ]
        cmd = ' '.join(cmd_args)
        cmd = 'mkdir -p prod/{0} && mv joboptions/{1} prod/{0}/ && cd prod/{0} && '.format(run_dir, jo) + cmd
        return cmd
    else:
        return None


def d3pd_cmd(run_dir, input_root):
    """
    EVNT -> NTUP_TRUTH COMMAND LINE
    """
    cmd_args = [
        'Reco_trf.py',
        'preExec="rec.doApplyAODFix.set_Value_and_Lock(False);from D3PDMakerConfig.D3PDMakerFlags import D3PDMakerFlags;D3PDMakerFlags.TruthWriteEverything=True"',
        'inputEVNTFile={0}'.format(input_root),
        'outputNTUP_TRUTHFile=ntup.{0}'.format(input_root),
        ]
    cmd = ' '.join(cmd_args)
    cmd = 'cd {0} && '.format(run_dir) + cmd
    return cmd


from . import FLAT_DRIVER
def flat_cmd(run_dir, input_root):
    """
    TRUTH -> FLAT TREE COMMAND LINE
    """
    output_root = 'flat.' + input_root
    cmd = 'python {0} {1} {2}'.format(
        FLAT_DRIVER, input_root, output_root)
    cmd = 'cp {0} {1} && cd {1} && {2}'.format(FLAT_DRIVER, run_dir, cmd)
    if os.path.exists(os.path.join(run_dir, output_root)):
        return 'echo "Output already exists !"'
    return cmd

from . import FLAT_DRIVER_2, PROD_DIR, NTUPLE_DIR
def flat2_cmd(d, input_root):
    """
    TRUTH -> FLAT TREE COMMAND LINE
    """
    mode = get_dir_mode(d)
    mass = get_dir_mass(d)
    seed = get_dir_seed(d)
    output_root = 'flat_{0}{1}_seed{2}.root'.format(mode, mass, seed)
    input_abs = os.path.join(PROD_DIR, d, input_root)
    output_abs = os.path.join(NTUPLE_DIR, 'running', output_root)
    cmd = 'python {0} {1} {2}'.format(
        FLAT_DRIVER_2, input_abs, output_abs)
    if os.path.exists(output_abs):
        return 'echo "Output already exists !"'
    return cmd
