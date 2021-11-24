from app import I2CCom
import time

registers = {
    'pump1':0x06,
    'pump2':0x0E,
    'pump3':0x16,
    'pump4':0x1E
    }

class Pump():
    def __init__(self,reactor,name):
        self.reactor = reactor
        self.name = name
        self.register = registers[name]
        self.running = False

    def turnOn(self):
        I2CCom(self.reactor,'Pumps',0,8,self.register,1,0)
        self.running = True

    def turnOff(self):
        I2CCom(self.reactor,'Pumps',0,8,self.register,0,0)
        self.running = False


def test_pumps():
    pumps = ['pump1','pump2','pump3','pump4']
    for pump in pumps:
        p = Pump('M0',pump)
        if not p.running:
            p.turnOn()
            time.sleep(0.1)
            p.turnOff()
            time.sleep(1)
