from pump import Pump
from os.path import join
from os.path import exists
from os import mkdir
import pandas as pd
import time

def create_template(reactor_name):
    pump_names = ['pump1','pump2','pump3','pump4']
    if not exists('calibrations'):
        mkdir('calibrations')
    df = pd.DataFrame(columns=['before','after'],\
        index=pump_names)
    df.index.name = 'pump_names'
    df.to_csv(join('calibrations',reactor_name+'.csv'),index_label='pump_names')

def run_pumps(reactor_name,pump_interval):
    pump_names = ['pump1','pump2','pump3','pump4']
    for pump_name in pump_names:
        p = Pump(reactor_name,pump_name)
        p.turnOn()
        time.sleep(pump_interval)
        p.turnOff()