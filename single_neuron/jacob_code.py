"""Alex code"""
#code for midterm project

from brian2 import *
import numpy as np

#PARAMETERS

#Iext = constant input current that sets the excitability for the neuron
#Isyn = synaptic current - these depend on synapses and neuron - define later

El = 10.6*mV #leak reversal potential
ENa = 115*mV #sodium reversal potential
EK = -12*mV #potassium reversal potential 

gl = 2.7*np.pi*msiemens #leak current conductance 
gNamax = 1080*np.pi*msiemens #sodium current max conductance 
gK = 324*np.pi*msiemens #potassium current conductance

alphaA = 1.1 * nM-1ms-1 () #FIGURE OUT THESE UNITS
betaA = 0.19 * ms-1 #UNITS
alphaG = 5.0 * mM-1ms-1 #UNITS 
betaG = 0.30 * ms-1 #UNITS
Ea = 60*mV
EG = -20*mV 

#NEED TO FIGURE OUT MORPHOLOGY 
morpho = Section(diameter = )
morpho.area = 30*np.pi*um**2

# Typical equations
eqs = '''
Cm*dv/dt = gl * (El-v) + gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) + Iext + Isyn: amp/meter**2

dm/dt = alpham * (1-m) - betam * m : 1
dn/dt = alphan * (1-n) - betan * n : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = (0.1/mV) * (10*mV-v)/exprel((-v+25*mV)/(10*mV))/ms : Hz
betam = 4 * exp(-v/(18*mV))/ms : Hz
alphah = 0.07 * exp(-v/(20*mV))/ms : Hz
betah = 1/(exp((-v+30*mV) / (10*mV)) + 1)/ms : Hz
alphan = (0.01/mV) * (10*mV-v)/exprel((-v+10*mV)/(10*mV))/ms : Hz
betan = 0.125*exp(-v/(80*mV))/ms : Hz
gNa : siemens/meter**2

'''

neuronM = NeuronGroup(1, model = eqs, method = 'euler', ) 
#not sure if we should be using a spatial neuron? or how to set up the spiking thresholds etc in the HH model
#want one excitatory synapse from M to S 
#FIX mM-1 UNIT, how to express Vpre as presynaptic potential ??
MS = Synapses(neuronM,neuronS, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1*mM-1 /(1+exprel(-(Vpre-62*mV)/5*mV) : 1 #T is neurotransmitter  on entration in the synaptic flect 
              # --> Maybe mM-1 or 1/mM is actually 1 / concentration as this would cancel out neurotransmitter concentration units
              Isyn = g*r*(v-Ea)'''
              )

neuronS = NeuronGroup(1, model=eqs, method ='euler')
#want 1 excitatory synapse from S to I 
SI = Synapses(neuronM,neuronS, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1*mM-1 /(1+exprel(-(Vpre-62*mV)/5*mV) : 1 #T is neurotransmitter  on entration in the synaptic flect 
              # --> Maybe mM-1 or 1/mM is actually 1 / concentration as this would cancel out neurotransmitter concentration units
              Isyn = g*r*(v-Ea)'''
              )

neuronI = NeuronGroup(1, model=eqs, method ='euler')
#want 1 inhibitory synapse from I to S 
IS = Synapses(neuronI,neuronS, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1*mM-1 /(1+exprel(-(Vpre-62*mV)/5*mV) : 1 #T is neurotransmitter  on entration in the synaptic flect 
              # --> Maybe mM-1 or 1/mM is actually 1 / concentration as this would cancel out neurotransmitter concentration units
              Isyn = g*r*(v-Ea)'''
              )

indices = array([0, 2, 1])
times = array([1, 2, 3])*ms
G = SpikeGeneratorGroup(3, indices, times)

neuron.v = 0*mV
neuron.h = 1
neuron.m = 0
neuron.n = .5
neuron.gNa = gNamax

M = StateMonitor(neuron, 'v', record=True)




"""Code from Feb 11"""
#FROM PAPER

r = ? #resistance
i(t) = ? #input current - modelled as step function
a = 10*mV #oscillation amplitude
f = ? #oscillation frequency (determined by entrainment function)
b = ? #phase shift (determined by entrainment function)
taum = 20*ms #membrane time constant
vrest = -70*mV #membrane resting potential
vth = -50*mV #membrane threshold
#El = -80*mV #leak reversal potential
#Ee = 0 * mV

Fmin = 500*Hz; Fmax = 32000*Hz
Fsd = 2500*Hz

pf = ? #preferred frequency
#stim = single attended stimululs
#stim.rhythm = rhythm of stimululs in Hz
#31.5*Hz = normalization factor
f(stim) = (|pf-stim.freq|)/31.5*kHz*(1/(2*stim.rhythm))

eqs = '''
dv/dt = -((v-vrest)+r*i(t)+a*sin(2*np.pi*f*t+b))/taum : volt
a : 1
'''

input = a * sin(2 * np.pi) #step current which, when turned on is a sine wave

neuron = NeuronGroup(1, model=eqs, threshold='v > vth', reset='v = vrest',
                      method='euler')



#code for midterm project

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt

start_scope()
A = 2.5
f = 10*Hz
tau = 5*ms
eqs = '''
dv/dt = (I-v)/tau : 1
I = A*sin(2*pi*f*t) : 1
'''
G = NeuronGroup(1, eqs, threshold='v>1', reset='v=0', method='euler')
M = StateMonitor(G, variables=True, record=True)
run(200*ms)
plot(M.t/ms, M.v[0], label='v')
plot(M.t/ms, M.I[0], label='I')
xlabel('Time (ms)')
ylabel('v')
legend(loc='best')
plt.show()
