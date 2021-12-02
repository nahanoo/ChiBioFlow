from pump import Pump
from os.path import join
import pandas as pd
from .board import Board


def create_template():
    pump_names = ['pump1', 'pump2', 'pump3', 'pump4']
    df = pd.DataFrame(columns=['name', 'before', 'after', 'fit'],
                      index=pump_names)
    df['name'] = pump_names
    df.to_csv('calibrations.csv', index=False)


def run_pumps(t):
    pumps = {
        'pump1': 'clock',
        'pump2': 'counter_clock',
        'pump3': 'clock',
        'pump4': 'clock'
    }
    b = Board(pumps)
    for pump in b.pumps:
        pump.run_time(t)


def get_fit(t):
    df = pd.read_csv('calibrations.csv')
    df['fit'] = (df['after']-df['before'])/t
