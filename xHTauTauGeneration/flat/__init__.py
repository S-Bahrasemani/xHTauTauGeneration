import rootpy
import logging
import os
import ROOT

ROOT.gROOT.SetBatch(True)

rootpy.log.basic_config_colorized()

log = logging.getLogger('flat')
if not os.environ.get("DEBUG", False):
    log.setLevel(logging.INFO)

if hasattr(logging, 'captureWarnings'):
    logging.captureWarnings(True)
