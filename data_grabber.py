from glob import glob
from os import mkdir
from os.path import join, exists
import argparse
from subprocess import call
from shutil import move


def parse_args():
    """Parsing variables for plotting."""
    parser = argparse.ArgumentParser(
        description="Data grabber for the ChiBio.")
    parser.add_argument("fname", help="fname consisting of date and time")
    parser.add_argument("experiment", help='experiment name')
    return parser.parse_args()


def get_files(fname):
    cmd = ['scp','root@192.168.7.2:chibio:'+fname+'*',join('data','tmp')]
    call(' '.join(cmd),shell=True)
    

def copy_files(experiment):
    for f in glob(join('data','tmp','*')):
        print(f)
        m = 'M'+f.split('M')[-1][0]
        print(m)
        dir = join('data',experiment,m)
        if not exists(dir):
            mkdir(dir)
        move(f,join('data',experiment,m))
        


fname = parse_args().fname
experiment = parse_args().experiment
dir = join('data',experiment)
if not exists(dir):
    mkdir(dir)
get_files(fname)
copy_files(experiment)
