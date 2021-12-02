from subprocess import call
from os.path import exists
from os.path import join
from os import mkdir
from .board import Board


class Experiment():
    def __init__(self, experiment, reactor):
        self.experiment = experiment
        self.reactor = reactor
        pumps = {
            'pump1': 'clock',
            'pump2': 'counter_clock',
            'pump3': 'clock',
            'pump4': 'clock'
        }
        self.board = Board(pumps)

    def get_data(self):
        path = join('data', self.experiment, self.reactor)
        if not exists(path):
            mkdir(path)

        src = 'eulrich@curnagl.dcsr.unil.ch:'+join(path, '*')
        cmd = ['scp', src, path]
        call(' '.join(cmd), shell=True)
