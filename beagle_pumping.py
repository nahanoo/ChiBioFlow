from app import I2CCom
import time


class Reactor():
    def __init__(self, reactor_name, pumps):
        self.reactor_name = reactor_name
        self.pumps = {pump_name: Pump(reactor_name, pump_name, direction)
                      for pump_name, direction in pumps.items()}


class Pump():
    def __init__(self, reactor_name, pump_name, direction):
        registers = {'clock':
                     {
                         'pump1': 0x06,
                         'pump2': 0x0E,
                         'pump3': 0x16,
                         'pump4': 0x1E
                     },
                     'counter_clock':
                     {
                         'pump1': 0x0A,
                         'pump2': 0x12,
                         'pump3': 0x1A,
                         'pump4': 0x22
                     }
                     }
        self.reactor_name = reactor_name
        self.running = False
        self.direction = direction
        if self.direction == 'clock':
            self.register = registers['clock'][pump_name]
        if self.direction == 'counter_clock':
            self.register = registers['counter_clock'][pump_name]

    def turnOn(self):
        I2CCom(self.reactor_name, 'Pumps', 0, 8, self.register, 1, 0)
        self.running = True

    def turnOff(self):
        I2CCom(self.reactor_name, 'Pumps', 0, 8, self.register, 0, 0)
        self.running = False

    def run_time(self, t):
        self.turnOn()
        time.sleep(t)
        self.turnOff()


def create_reactors():
    pump_mappings = {
        'M0': {
            'pump1': 'counter_clock',
            'pump2': 'counter_clock',
            'pump3': 'counter_clock',
            'pump4': 'counter_clock'
        }
    }
    reactors = dict()
    for reactor_name, pumps in pump_mappings.items():
        reactors[reactor_name] = Reactor(reactor_name, pumps)
    return reactors
