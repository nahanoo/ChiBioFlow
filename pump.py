from app import I2CCom
import time
from calibration import get_fit

class Reactor():
    def __init__(self,reactor_name,pump_names):
        self.reactor_name = reactor_name
        self.pumps = {pump_name:Pump(reactor_name,pump_name) for pump_name in pump_names}
        for pump_name,pump in self.pumps.items():
            pump.fit = get_fit(reactor_name,pump_name,5)

class Pump():
    def __init__(self,reactor_name,pump_name):
        registers = {
            'pump1':0x06,
            'pump2':0x0E,
            'pump3':0x16,
            'pump4':0x1E
            }
        self.reactor_name = reactor_name
        self.pump_name = pump_name
        self.register = registers[pump_name]
        self.running = False
        self.fit = None

    def turnOn(self):
        I2CCom(self.reactor,'Pumps',0,8,self.register,1,0)
        self.running = True

    def turnOff(self):
        I2CCom(self.reactor,'Pumps',0,8,self.register,0,0)
        self.running = False

    def inject_volume(self,ml):
        runtime = ml/self.fit
        self.turnOn()
        time.sleep(runtime)
        self.turnOff

def create_reactors():
    pump_mappings = {
        'M0':['pump1','pump2','pump3','pump4']
    }
    reactors = dict()
    for reactor_name, pump_names in pump_mappings.items():
        reactors[reactor_name] = Reactor(reactor_name,pump_names)
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