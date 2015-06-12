"""
"""
import os
from argparse import ArgumentParser

from rootpy.io import root_open
from rootpy.plotting import Hist, Canvas
from rootpy.plotting.style import set_style
from rootpy.tree import Cut
from rootpy import ROOT
from rootpy import asrootpy
from ROOT import TLatex

from plotting.variables import VARIABLES, get_label
from gen import MASSES
from gen.cmd import get_dir_mode, get_dir_mass

parser = ArgumentParser(os.path.basename(__file__))
parser.add_argument('--var', type=str, default='')
parser.add_argument('--tree-name', type=str, default='Tree')
parser.add_argument('--mode', nargs='*', default=['VBF', 'gg'], choices=['VBF', 'gg'])
parser.add_argument('--mass', type=int, default=-1)
parser.add_argument('--cut', type=str, default='')
args = parser.parse_args()

if args.mass in MASSES:
    masses = [args.mass]
elif args.mass == -1:
    masses = MASSES
else:
    raise RuntimeError('Wrong mass argument!')

if args.var in VARIABLES.keys():
    variables = {args.var: VARIABLES[args.var]}
elif args.var == '':
    variables = VARIABLES
else:
    raise RuntimeError('Wrong input variable')



set_style('ATLAS', shape='rect')
ROOT.gROOT.SetBatch(True)

hist_array = {}
cut = Cut(args.cut)
DIR_PATH = './prod/'
# DIR_PATH = '/cluster/data03/qbuat/HtautauGeneration/prod/'
for dirpath, dirnames, files in os.walk(DIR_PATH):
    if dirpath != DIR_PATH:
        break
        
    dirnames = filter(lambda d: get_dir_mode(d) in args.mode, dirnames)
    if args.mass != -1:
        dirnames = filter(lambda d: get_dir_mass(d) in masses, dirnames)
    for d in dirnames:
        for f in os.listdir(os.path.join(DIR_PATH, d)):
            if f[0:5] == 'flat2' and '.root' in f:
                with root_open(os.path.join(
                        DIR_PATH, d, f)) as froot:
                    tree = froot[args.tree_name]
                    for key, var_info in variables.items():
                        h = Hist(var_info['bins'], var_info['range'][0], var_info['range'][1])
                        expr = '{0}>>{1}{2}'
                        if 'scale' in var_info.keys():
                            expr = expr.format(
                                var_info['name'] + '*' + str(var_info['scale']),
                                h.name,
                                (var_info['bins'], var_info['range'][0], var_info['range'][1]))
                        else:
                            expr = expr.format(
                                var_info['name'], h.name,
                                (var_info['bins'], var_info['range'][0], var_info['range'][1]))
                        print expr
                        tree.Draw(expr, cut, var_info['style'])
                        h = asrootpy(ROOT.gPad.GetPrimitive(h.name))
                        h.xaxis.title = get_label(var_info)
                        h.yaxis.title = 'Number of Events'
                        if h.integral() != 0:
                            h /= h.integral() 
                        h.SetLineWidth(2)
                        h.yaxis.SetRangeUser(0, 0.5)
                        if 'binlabels' in var_info.keys():
                            for ib, lab in enumerate(var_info['binlabels']):
                                h.xaxis.SetBinLabel(ib + 1, lab)

                        c = Canvas()
                        # c.SetLogy()
                        c.SetTopMargin(0.08)
                        h.Draw('HIST')
                        label = ROOT.TLatex(
                            c.GetLeftMargin() + 0.05, 1 - c.GetTopMargin() + 0.02,
                            '{0}H #rightarrow #tau#tau (m_{{H}} = {1} GeV) - Hist Integral = {2}'.format(
                                get_dir_mode(d), get_dir_mass(d), h.Integral()))
                        label.SetNDC()
                        label.SetTextSize(20)
                        label.Draw('same')
                        c.SaveAs('./plots/{0}_mass{1}_mode{2}.png'.format(
                                key, get_dir_mass(d), get_dir_mode(d)))
