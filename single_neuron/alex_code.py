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
MS = Synapses(M,S, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1*mM-1 /(1+exprel(-(Vpre-62*mV)/5*mV) : 1 #T is neurotransmitter  on entration in the synaptic flect
              Isyn = g*r*(v-Ea)'''
              )

neuronS = 
#want 1 excitatory synapse from S to I 

neuronI = 
#want 1 inhibitory synapse from I to S 


neuron.v = 0*mV
neuron.h = 1
neuron.m = 0
neuron.n = .5
neuron.gNa = gNamax

M = StateMonitor(neuron, 'v', record=True)

