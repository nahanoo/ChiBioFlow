#from app import I2CCom

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
        self.register = registers

p = Pump()
