#code for midterm project

from brian2 import *
import numpy as np

#PARAMETERS

defaultclock.dt = 0.01*ms

morpho = Cylinder(length=10*cm, diameter=2*238*um, n=1000, type='axon')

El = 10.6*mV #leak reversal potential
ENa = 115*mV #sodium reversal potential
EK = -12*mV #potassium reversal potential 

gl = 0.3*msiemens/cm**2
gNa0 = 120*msiemens/cm**2
gK = 36*msiemens/cm**2

alphaA = 1.1 * nM**-1 * ms**-1 
betaA = 0.19 * ms**-1 
alphaG = 5.0 * mM**-1 * ms**-1 
betaG = 0.30 * ms** -1 
EA = 60*mV
EG = -20*mV 

# MODEL EQUATIONS 
eqs = '''
Im = gl * (El-v) + gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) : amp/meter**2
I : amp (point current) # applied current
dm/dt = alpham * (1-m) - betam * m : 1
dn/dt = alphan * (1-n) - betan * n : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = (0.1/mV) * 10*mV/exprel((-v+25*mV)/(10*mV))/ms : Hz
betam = 4 * exp(-v/(18*mV))/ms : Hz
alphah = 0.07 * exp(-v/(20*mV))/ms : Hz
betah = 1/(exp((-v+30*mV) / (10*mV)) + 1)/ms : Hz
alphan = (0.01/mV) * 10*mV/exprel((-v+10*mV)/(10*mV))/ms : Hz
betan = 0.125*exp(-v/(80*mV))/ms : Hz
gNa : siemens/meter**2
'''

#CREATE NEURONS

neuronM = SpatialNeuron(morphology=morpho, model=eqs, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)
neuronM.v = 0*mV
neuronM.h = 1
neuronM.m = 0
neuronM.n = .5
neuronM.I = 0*amp
neuronM.gNa = gNa0
monM = StateMonitor(neuronM, 'v', record=True)
spikesM = SpikeMonitor(neuronM)


neuronS = SpatialNeuron(morphology=morpho, model=eqs, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)
neuronS.v = 0*mV
neuronS.h = 1
neuronS.m = 0
neuronS.n = .5
neuronS.I = 0*amp
neuronS.gNa = gNa0
monS = StateMonitor(neuronS, 'v', record=True)
spikesS = SpikeMonitor(neuronS)

neuronI = SpatialNeuron(morphology=morpho, model=eqs, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)
neuronI.v = 0*mV
neuronI.h = 1
neuronI.m = 0
neuronI.n = .5
neuronI.I = 0*amp
neuronI.gNa = gNa0
monI = StateMonitor(neuronS, 'v', record=True)
spikesI = SpikeMonitor(neuronS)

#CREATE SYNAPSES

#one excitatory synapse from M to S
MS = Synapses(neuronM,neuronS, '''
              dr/dt = (alphaA*T*(1-r)-betaA*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              
              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')
MS.connect()

#one excitatory synapse from S to I 
SI = Synapses(neuronS,neuronI, '''
              dr/dt = (alphaA*T*(1-r)-betaA*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              
              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')
SI.connect()

#one inhibitory synapse from I to S 
IS = Synapses(neuronI,neuronS, '''
              dr/dt = (alphaG*T*(1-r)-betaG*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              
              on_post = '''
              I = gK*r*(v_post - EG)*(meter**2) ''')
IS.connect()

run(100*ms)


