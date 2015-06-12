from pbs import qsub

def get_setup(filename):
    with open(filename) as f:
        return ' && '.join([line.strip() for line in f.readlines()])
