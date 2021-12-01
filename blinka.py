import board
import busio
from adafruit_bus_device.i2c_device import I2CDevice
import time

i2c = busio.I2C(board.SCL,board.SDA)
device = I2CDevice(i2c,0x61)

with device:
    device.write(bytes([0,0]))
    device.write(bytes([253,0]))
    device.write(bytes([6,1]))
    time.sleep(1)
    device.write(bytes([6,0]))