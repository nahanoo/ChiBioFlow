from pump import Pump
from os.path import join
from os.path import exists
from os import mkdir
import pandas as pd
import time

pump_names = ['pump1','pump2','pump3','pump4']

def create_template(reactor_name):
    if not exists('calibrations'):
        mkdir('calibrations')
    df = pd.DataFrame(columns=['before','after'],\
        index=pump_names)
    df.index.name = 'pump_names'
    df.to_csv(join('calibrations',reactor_name+'.csv'),index_label='pump_names')

def run_pumps(reactor_name,pump_interval):
    for pump_name in pump_names:
        p = Pump(reactor_name,pump_name)
        p.turnOn()
        time.sleep(pump_interval)
        p.turnOff()

def get_fit(reactor_name,pump_name,pump_interval):
    df = pd.read_csv(join('calibrations',reactor_name+'.csv'),index_col='pump_names')
    x = pump_interval
    y = df.loc[pump_name]['after']
    b = df.loc[pump_name]['before']
    m = (y-b)/x
    return m