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

# Synapses parameters
gA = 10 * msiemens
gG = 20 * msiemens # variable 
EA = 0*mV #shifted down 60mV from paper
EG = -80*mV #shifted down 60mV from paper
alphaA = 1.1 * mM**-1*ms**-1 #FIGURE OUT THESE UNITS
betaA = 0.19 * ms**-1 #UNITS
alphaG = 5.0 * mM**-1*ms**-1 #UNITS 
betaG = 0.30 * ms**-1 #UNITS

# The model
eqs = Equations('''
dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + Iext + Isyn)/Cm : volt
dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
Iext : amp
Isyn: amp
''')

# Threshold and refractoriness are only used for spike counting
group = NeuronGroup(num_neurons, eqs,
                    threshold='v > -40*mV',
                    refractory='v > -40*mV',
                    method='exponential_euler')
group.v = El
group.Iext = Iext_values
#group.I = '0.7*nA * i / num_neurons'

A = Synapses(group,group, 
             model = ''' 
               dr/dt = alphaA*T*(1-r)-betaA*r : 1
               T = 1 * mM/(1+exprel(-((0.001*v_pre/mV)-62)/5)) : mM
               IA = gA*r*(v - EA) : amp''', 
               on_pre = 'Isyn += IA', method = 'euler'
               )
              

# Excitatory synapses from Sender to Reveiver (0 to 1) 
# and Receiver to Interneuron (1 to 2)          
A.connect(i = [0, 1], j = [1, 2])

# Inhibitory synapse
G = Synapses(group,group, model = ''' 
             dr/dt = alphaG*T*(1-r)-betaG*r : 1
             T = 1 * mM/(1+exprel(-((0.001*v_pre/mV)-62)/5)) : mM
             IG = gG*r*(v - EG) : amp''', 
             on_pre = 'Isyn += IG', method = 'euler'
             )
                          

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

plot(monitorV.v[1])
