from os import mkdir
from os.path import join, exists
import argparse
from subprocess import call
from shutil import copy
from os import listdir
from os import remove


def parse_args():
    """Parsing variables for plotting."""
    parser = argparse.ArgumentParser(
        description="Data grabber for the ChiBio.")
    parser.add_argument("fname", help="fname consisting of date and time")
    parser.add_argument("experiment", help='experiment name')
    return parser.parse_args()


def get_files(fname):
    cmd = ['scp', join('root@192.168.7.2:','root','chibio',fname+'*'), join('data', 'tmp')]
    call(' '.join(cmd), shell=True)


def copy_files(experiment):
    for f in listdir(join('data', 'tmp')):
        m = 'M'+f.split('M')[-1][0]
        dir = join('data', experiment, m)
        if not exists(dir):
            mkdir(dir)
        copy(join('data', 'tmp',f), join('data', experiment, m,f))
        remove(join('data', 'tmp',f))


fname = parse_args().fname
experiment = parse_args().experiment
dir = join('data', experiment)
if not exists(dir):
    mkdir(dir)
get_files(fname)
copy_files(experiment)
