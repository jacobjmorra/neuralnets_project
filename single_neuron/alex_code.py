#midterm code 

from brian2 import *
import numpy as np

start_scope()

num_neurons = 3
duration = 1*second

# Parameters
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2 * area
gl = 5e-5*siemens*cm**-2 * area
El = -65*mV
EK = -90*mV
ENa = 50*mV
g_na = 100*msiemens*cm**-2 * area
g_kd = 30*msiemens*cm**-2 * area
VT = -63*mV

# Constant current, determines neuron excitability
Iext_S = 28e-11 *amp
Iext_I = 28e-11 *amp
Iext_R = 28e-11 *amp # variable
Iext_values = [Iext_S, Iext_I, Iext_R]

w_e = 0.1*msiemens
w_i = 0.4*msiemens
Ee = 0*mV
Ei = -80*mV
taue = 5*ms
taui = 10*ms

# The model
eqs = Equations('''
dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + Iext + g_e*(Ee-v) + g_i*(Ei-v))/Cm : volt
dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
dg_e/dt = -g_e/taue : siemens
dg_i/dt = -g_i/taui : siemens
Iext : amp
''')

# Threshold and refractoriness are only used for spike counting
group = NeuronGroup(num_neurons, eqs,
                    threshold='v > -40*mV',
                    refractory='v > -40*mV',
                    method='exponential_euler')
group.v = El
group.Iext = Iext_values
#group.I = '0.7*nA * i / num_neurons'

A = Synapses(group,group, 'w: siemens', on_pre='g_e += w_e')
              

# Excitatory synapses from Sender to Reveiver (0 to 1) 
# and Receiver to Interneuron (1 to 2)          
A.connect(i = [0, 1], j = [1, 2])

# Inhibitory synapse
G = Synapses(group,group, 'w: siemens', on_pre='g_i += w_i')
                          

# # Inhibitory synapses from Interneuron to Reveiver (2 to 1) 
G.connect(i = [2], j = [1]) 

# Monitor variables
monitorV = StateMonitor(group, 'v', record=True)
monitorS = SpikeMonitor(group)

run(duration)

# plot(group.I/nA, monitor.count / duration)
# xlabel('I (nA)')
# ylabel('Firing rate (sp/s)')
# show()

plot(monitorV.v[0])
plot(monitorV.v[1])
plot(monitorV.v[2])