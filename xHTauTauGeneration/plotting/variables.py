
def get_label(variable):
    label = variable['root']
    if 'units' in variable.keys():
        label += ' [{0}]'.format(variable['units'])
    return label


VARIABLES = {
    'tau1_pt': {
        'name': 'tau1_pt',
        'root': 'Leading p_{T}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 40),
        'style': 'HIST',
        },
    'tau2_pt': {
        'name': 'tau2_pt',
        'root': 'Subleading p_{T}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 40),
        'style': 'HIST',
        },
    'tau1_vis_pt': {
        'name': 'tau1_vis_pt',
        'root': 'Leading p_{T, visible}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 40),
        'style': 'HIST',
        },
    'tau2_vis_pt': {
        'name': 'tau2_vis_pt',
        'root': 'Subleading p_{T, visible}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 40),
        'style': 'HIST',
        },

    'taum_pt': {
        'name': 'taum_pt',
        'root': '#tau^{-} p_{T}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 120),
        'style': 'HIST',
        },
    'taup_eta': {
        'name': 'taup_eta',
        'root': '#tau^{+} #eta',
        'bins': 12,
        'range': (-3, 3),
        'style': 'HIST',
        },
    'taum_eta': {
        'name': 'taum_eta',
        'root': '#tau^{-} #eta',
        'bins': 12,
        'range': (-3, 3),
        'style': 'HIST',
        },
    'lep1_pt': {
        'name': 'lep1_pt',
        'root': 'Leading p_{T}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 120),
        'style': 'HIST',
        },
    'lep2_pt': {
        'name': 'lep2_pt',
        'root': 'Subleading p_{T}',
        'units': 'GeV',
        'scale': 0.001,
        'bins': 10,
        'range': (0, 120),
        'style': 'HIST',
        },
    'lep1_eta': {
        'name': 'lep1_eta',
        'root': 'Leading #eta',
        'units': 'GeV',
        'bins': 12,
        'range': (-3, 3),
        'style': 'HIST',
        },
    'lep2_eta': {
        'name': 'lep2_eta',
        'root': 'Subleading #eta',
        'units': 'GeV',
        'bins': 12,
        'range': (-3, 3),
        'style': 'HIST',
        },
    'channel': {
        'name': 'hadhad + 2 * lephad + 3 * leplep',
        'root': 'Channel',
        'bins': 3,
        'range': (0.5, 3.5),
        'binlabels': ['hadhad', 'lephad', 'leplep'],
        'style': 'HIST',
        },
}
