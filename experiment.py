from subprocess import call
from os.path import exists
from os.path import join
from os import mkdir
import pandas as pd


class Experiment():
    def __init__(self, experiment):
        self.experiment = experiment

    def get_data(self,reactor):
        path = join('data', self.experiment, reactor)
        src = 'eulrich@curnagl.dcsr.unil.ch:'+join(path, '*')
        cmd = ['scp', src, path]
        call(' '.join(cmd), shell=True)

    def write_data(self,reactor):
        path = join('data', self.experiment, reactor, 'injections')
        if not exists(path):
            mkdir(path)