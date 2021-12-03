from experiment import Experiment
from reactor import Reactor

e = Experiment('Test')
pumps = {
    'pump1': 'clock',
    'pump2': 'counter_clock',
    'pump3': 'clock',
    'pump4': 'clock'
}
r = Reactor('M0',pumps)
e.get_data()
r.cycle(10)