# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 13:03:20 2021

@author: Alex
"""

from brian2 import *
from scipy import stats

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

# Typical equations
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

MS = Synapses(neuronM,neuronS, '''
              dr/dt = (alphaA*T*(1-r)-betaA*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              
              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')
MS.connect()


run(5*ms, report='text')
neuronM.I[0] = 1*uA # current injection at one end
neuronS.I[0] = 1*uA # current injection at one end
'''
run(50*ms, report='text')
neuronM.I[0] = 1*uA # current injection at one end
neuronS.I[0] = 1*uA # current injection at one end
run(3*ms)
neuronM.I = 0*amp
neuronS.I = 0*amp
run(50*ms, report='text')
'''

plot()
for i in range(10):
    plot(monM.t/ms, monM.v.T[:, i*100]/mV)
ylabel('v')

plot()
for i in range(10):
    plot(monS.t/ms, monS.v.T[:, i*100]/mV)
ylabel('v')
