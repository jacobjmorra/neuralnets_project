s#code for midterm project

from brian2 import *
import numpy as np

start_scope()

#PARAMETERS

# Morphology parameters???
#area = 30 * 30 * pi * umeter**2

# Reversal potentials
El = 10.6*mV #leak reversal potential
ENa = 115*mV #sodium reversal potential
EK = -12*mV #potassium reversal potential 
EA = 60*mV 
EG = -20*mV 

# Conductances
gl = 2.7*np.pi*msiemens #leak current conductance 
gNa = 1080*np.pi*msiemens #sodium current max conductance 
gK = 324*np.pi*msiemens #potassium current max conductance
gA = 10 * msiemens
gG = 20 * msiemens # variable 

# Constant current, determines neuron excitability
Iext_S = 28e-11 *amp
Iext_I = 28e-11 *amp
Iext_R = 28e-11 *amp # variable
Iext = [Iext_S, Iext_I, Iext_R]

# Membrane capacitance
Cm = pi*9e-6 *farad

alphaA = 1.1 * mM**-1*ms**-1 #FIGURE OUT THESE UNITS
betaA = 0.19 * ms**-1 #UNITS
alphaG = 5.0 * mM**-1*ms**-1 #UNITS 
betaG = 0.30 * ms**-1 #UNITS

# Typical equations
# the membrane potential is determined by: leaky current, Na and K currents, 
# a constant current Iext0 that has specific values for each neuron
# and the synaptic currents.
eqs = '''
dv/dt = (gl * (El-v) + gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) + Iext0 + Isyn)/Cm : volt
dm/dt = alpham * (1-m) - betam * m : 1
dn/dt = alphan * (1-n) - betan * n : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = (0.1/mV) * (10*mV-v)/exprel((-v+25*mV)/(10*mV))/ms : Hz
betam = 4 * exp(-v/(18*mV))/ms : Hz
alphah = 0.07 * exp(-v/(20*mV))/ms : Hz
betah = 1/(exp((-v+30*mV) / (10*mV)) + 1)/ms : Hz
alphan = (0.01/mV) * (10*mV-v)/exprel((-v+10*mV)/(10*mV))/ms : Hz
betan = 0.125*exp(-v/(80*mV))/ms : Hz
Iext0 : amp
Isyn : amp
'''

# All the neurons are created in the same group. their role is determined by their connections
# Neuron with index 0 is Sender, index 1 is Receiver, index 2 is Interneuron
neurons = NeuronGroup(3, model = eqs, method = 'exponential_euler', threshold="v > -30*mV") 
neurons.Iext0 = Iext

# Excitatory synapse
A = Synapses(neurons,neurons, 
             model = ''' 
             dr/dt = alphaA*T*(1-r)-betaA*r : 1
             T = 1 * mM/(1+exprel(-((0.001*v_pre/mV)-62)/5)) : mM
             IA = gA*r*(v - EA) : amp''', 
             on_pre = 'Isyn_post += IA'
             )

# Excitatory synapses from Sender to Reveiver (0 to 1) 
# and Receiver to Interneuron (1 to 2)          
A.connect(i = [0, 1], j = [1, 2])

# Inhibitory synapse
G = Synapses(neurons,neurons, 
             model = ''' 
             dr/dt = alphaG*T*(1-r)-betaG*r : 1
             T = 1 * mM/(1+exprel(-((0.001*v_pre/mV)-62)/5)) : mM
             IG = gG*r*(v - EG) : amp''', 
             on_pre = 'Isyn_post += IG'
             )

# Inhibitory synapses from Interneuron to Reveiver (2 to 1) 
G.connect(i = [2], j = [1]) 

# Monitor variables
M = StateMonitor(neurons, 'v', record=True)
S = SpikeMonitor(neurons, record=True)

# Run simulation
run( 1 * second, report='text' )

plot(M.v[0])
plot(M.v[1])
plot(M.v[2])
