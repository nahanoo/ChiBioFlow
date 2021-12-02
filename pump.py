from datetime import time
import pandas as pd


class Pump():
    def __init__(self, name, device, register):
        self.name = name
        self.device = device
        self.register = register
        self.running = False
        self.fit = self.get_fit()

    def get_fit(self):
        df = pd.read_csv('calibrations.csv', index_col='name')
        return df.loc[self.name]['fit']

    def turn_on(self):
        if self.running == False:
            with self.device:
                self.device.write(bytes([self.register, 1]))
        self.running = True

    def turn_off(self):
        with self.device:
            self.device.write(bytes([self.register, 0]))
        self.running = False

    def run_time(self, t):
        self.turn_on()
        time.sleep(t)
        self.turn_off()

    def inject_volume(self, ml):
        runtime = ml/self.fit
        self.run_time(runtime)
