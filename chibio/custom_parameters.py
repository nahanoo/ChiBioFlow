##### PARAMETERS #####

# Global clock time         Time in seconds for each cycle. Keep more than maximum of what one dilution cyle will take. Default is 120.
# Default stir time         Time in seconds for pausing after stirring when waiting for settling. Default is 10.
# Reactor Acvitity          1 if on and 0 if off. Default is 0
# reactor running cycles     Multiples of global clock time at which custom program of reactor should run. Default is 5
# number of dilutions        Number of dilution cycles to be done in the reactor. Default is 2
# Pump order                 Order of pumps. Default is ['Pump1','Pump2','Pump3','Pump4']
# Pump direction             1 or -1. Entered as dictionary for the four pumps. Default is [1,1,1,1]
# Pump flow time             Time for the pump to flow in seconds. Entered as list for the four pumps in the pump order. Default is [2,2,2,2] 
                             ### NOTE: Keep the flow time of outflow pump higher than inflow pump to avoid overflow
# Pump pause time            Time for the pump to pause in seconds when stirring of reactor occurs. Entered as list for the four pumps in the pump order. Default is [2,2,2,2]

global_clock_time = 1200
default_stir_time = 10

# Each reactor will have all the pump parameters. Keep activity of the reactor 0 if not in use.

### Reactor 1 : M1 ###

m1_activity = 1
m1_running_cycles = 1
m1_dilutions = 2
m1_pump_order = ['Pump4','Pump1','Pump2','Pump3']
m1_pump_direction = {'Pump1':1,'Pump2':-1,'Pump3':1,'Pump4':1}
m1_pump_flow_time = {'Pump1':10,'Pump2':15,'Pump3':15,'Pump4':15}
m1_pump_pause_time = {'Pump1':0,'Pump2':0,'Pump3':0,'Pump4':0}

### Reactor 2 : M2 ###

m2_activity = 1
m2_running_cycles = 1
m2_dilutions = 1
m2_pump_order = ['Pump4','Pump1','Pump2','Pump3']
m2_pump_direction = {'Pump1':1,'Pump2':-1,'Pump3':1,'Pump4':1}
m2_pump_flow_time = {'Pump1':10,'Pump2':15,'Pump3':15,'Pump4':15}
m2_pump_pause_time = {'Pump1':0,'Pump2':0,'Pump3':0,'Pump4':0}

### Reactor 3 : M3 ###
m3_activity = 1
m3_running_cycles = 1
m3_dilutions = 1
m3_pump_order = ['Pump4','Pump1','Pump2','Pump3']
m3_pump_direction = {'Pump1':1,'Pump2':-1,'Pump3':1,'Pump4':1}
m3_pump_flow_time = {'Pump1':10,'Pump2':15,'Pump3':15,'Pump4':15}
m3_pump_pause_time = {'Pump1':0,'Pump2':0,'Pump3':0,'Pump4':0}

### Reactor 4 : M4 ###
m4_activity = 1
m4_running_cycles = 1
m4_dilutions = 2
m4_pump_order = ['Pump1','Pump2','Pump3','Pump4']
m4_pump_direction = {'Pump1':1,'Pump2':-1,'Pump3':1,'Pump4':1}
m4_pump_flow_time = {'Pump1':3,'Pump2':5,'Pump3':5,'Pump4':5}
m4_pump_pause_time = {'Pump1':0,'Pump2':0,'Pump3':10,'Pump4':0}

### Reactor 5 : M5 ###
m5_activity = 0
m5_running_cycles = 5
m5_dilutions = 2
m5_pump_order = ['Pump1','Pump2','Pump3','Pump4']
m5_pump_direction = {'Pump1':1,'Pump2':1,'Pump3':1,'Pump4':1}
m5_pump_flow_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}
m5_pump_pause_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}

### Reactor 6 : M6 ###
m6_activity = 0
m6_running_cycles = 5
m6_dilutions = 2
m6_pump_order = ['Pump1','Pump2','Pump3','Pump4']
m6_pump_direction = {'Pump1':1,'Pump2':1,'Pump3':1,'Pump4':1}
m6_pump_flow_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}
m6_pump_pause_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}

### Reactor 7 : M7 ###
m7_activity = 0
m7_running_cycles = 5
m7_dilutions = 2
m7_pump_order = ['Pump1','Pump2','Pump3','Pump4']
m7_pump_direction = {'Pump1':1,'Pump2':1,'Pump3':1,'Pump4':1}
m7_pump_flow_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}
m7_pump_pause_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}

### Reactor 8 : M8 ###
m8_activity = 0
m8_running_cycles = 5
m8_dilutions = 2
m8_pump_order = ['Pump1','Pump2','Pump3','Pump4']
m8_pump_direction = {'Pump1':1,'Pump2':1,'Pump3':1,'Pump4':1}
m8_pump_flow_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}
m8_pump_pause_time = {'Pump1':2,'Pump2':2,'Pump3':2,'Pump4':2}




