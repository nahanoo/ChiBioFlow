import board
import busio
from adafruit_bus_device.i2c_device import I2CDevice
from pump import Pump
import time


class Reactor():
    def __init__(self, name, pumps):
        self.name = name
        self.registers = {'clock':
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

        i2c = busio.I2C(board.SCL, board.SDA)
        self.device = I2CDevice(i2c, 0x61)

        self.boot()
        self.init_pumps(pumps)

    def boot(self):
        with self.device:
            self.device.write(bytes([0, 0]))
            self.device.write(bytes([253, 0]))

    def init_pumps(self, pumps):
        self.pumps = {name: None for name in pumps.keys()}
        for name, direction in pumps.items():
            register = self.registers[direction][name]
            self.pumps[name] = Pump(name, self.device, register)

    def cycle(self, cycle_time):
        self.pumps['pump1'].inject_volume(1)
        time.sleep(cycle_time)