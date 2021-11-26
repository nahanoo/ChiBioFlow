from app import I2CCom
import time
import pandas as pd
from os.path import join

class Reactor():
    def __init__(self,reactor_name,pump_names,directions):
        self.reactor_name = reactor_name
        self.pumps = {pump_name:Pump(reactor_name,pump_name,direction) \
            for pump_name,direction in zip(pump_names,directions)}

class Pump():
    def __init__(self,reactor_name,pump_name,direction):
        registers = {'clock':
            {
                'pump1':0x06,
                'pump2':0x0E,
                'pump3':0x16,
                'pump4':0x1E
            },
            'counter_clock':
            {
                'pump1':0x0A,
                'pump2':0x12,
                'pump3':0x1A,
                'pump4':0x22
            }
        }
        self.reactor_name = reactor_name
        self.pump_name = pump_name
        self.running = False
        self.fit = self.get_fit()
        self.direction = direction
        if self.direction == 'clock':
            self.register = registers['clock'][pump_name]
        if self.direction == 'counter_clock':
            self.register = registers['counter_clock'][pump_name]

    def get_fit(self):
        df = pd.read_csv(join('calibrations',self.reactor_name+'.csv'),index_col='pump_names')
        x = 10
        y = df.loc[self.pump_name]['after']
        b = df.loc[self.pump_name]['before']
        m = (y-b)/x
        return m

    def turnOn(self):
        I2CCom(self.reactor_name,'Pumps',0,8,self.register,1,0)
        self.running = True

    def turnOff(self):
        I2CCom(self.reactor_name,'Pumps',0,8,self.register,0,0)
        self.running = False

    def run_time(self,t):
        self.turnOn()
        time.sleep(t)
        self.turnOff()

    def inject_volume(self,ml):
        runtime = ml/self.fit
        self.turnOn()
        time.sleep(runtime)
        self.turnOff()

def create_reactors():
    pump_mappings = {
        'M0':['pump1','pump2','pump3','pump4']
    }
    directions = {
        'M0':{
            'pump1':'counter_clock',
            'pump2':'clock',
            'pump3':'counter_clock',
            'pump4':'counter_clock'
        }
    }
    reactors = dict()
    for reactor_name, pump_names in pump_mappings.items():
        reactors[reactor_name] = \
            Reactor(reactor_name,pump_names,directions[reactor_name].values())
    return reactors
    
def test_pumps():
    pump_names = ['pump1','pump2','pump3','pump4']
    for pump_name in pump_names:
        p = Pump('M0',pump_name)
        if not p.running:
            p.turnOn()
            time.sleep(0.1)
            p.turnOff()
            time.sleep(1)

r = create_reactors()
m = r['M0']
