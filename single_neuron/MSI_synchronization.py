#code for midterm project

from brian2 import *
import numpy as np

#PARAMETERS

#Iext = constant input current that sets the excitability for the neuron
#Isyn = synaptic current - these depend on synapses and neuron - define in synapse equations?

El = 10.6*mV #leak reversal potential
ENa = 115*mV #sodium reversal potential
EK = -12*mV #potassium reversal potential

#note: in original HH model from Brian, these units were msiemens/cm**2 and had different values (what does this change?)
gl = 0.3*msiemens/cm**2
gNa0 = 120*msiemens/cm**2
gK = 36*msiemens/cm**2

#note: these are given in paper but not sure what these units are or how to call them in brian
#Jacob note: this is "molar" in brian2 --- https://brian2.readthedocs.io/en/stable/user/equations.html
alphaA = 1.1 * nmolar**-1 * ms**-1 #FIGURE OUT THESE UNITS --> DONE
betaA = 0.19 * ms**-1 #UNITS? --> DONE
alphaG = 5.0 * mmolar**-1 * ms**-1 #UNITS? --> DONE
betaG = 0.30 * ms**-1 #UNITS? -->DONE
Ea = 60*mV
EG = -20*mV

#NEED TO FIGURE OUT MORPHOLOGY (do we need this?? they don't mention anything other than the area listed below)
#it's not needed to my knowledge, but a motor neuron is around 100 microns or 0.1 mm in diameter
#morpho = Section(diameter = 0.1 * mm * np.array([1,1]))
morpho = Soma(50*um)  # chosen for a target Rm -> taken from https://brian2.readthedocs.io/en/2.0rc/examples/frompapers.Brette_2012.Fig1.html
#morpho.area = 30*np.pi*um**2

# MODEL EQUATIONS - not sure if the units line up here!!!! most of this is from the HH example on brian, but i've filled things
#in / made small changes so they match the paper

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

neuronS = SpatialNeuron(morphology=morpho, model=eqs, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)

neuronI = SpatialNeuron(morphology=morpho, model=eqs, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)

#want one excitatory synapse from M to S
#FIX mM-1 UNIT, how to express Vpre as presynaptic potential ??
MS = Synapses(neuronM,neuronS, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',

              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')

#want 1 excitatory synapse from S to I
SI = Synapses(neuronS,neuronI, model = ''' 
              dr/dt = alphaA*T*(1-2)-betaA*r : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',

              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')

#want 1 inhibitory synapse from I to S
IS = Synapses(neuronI,neuronS, model = ''' 
              dr/dt = alphaG*T*(1-2)-betaG*r : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',

              on_post = '''
              I = gK*r*(v_post - EG)*(meter**2) ''')

MS.connect()
SI.connect()
IS.connect()

"""
#these initial values come from the example on brian
neuronM.v = 0*mV
neuronM.h = 1
neuronM.m = 0
neuronM.n = .5
neuronM.gNa = gNa0

neuronS.v = 0*mV
neuronS.h = 1
neuronS.m = 0
neuronS.n = .5
neuronS.gNa = gNa0

neuronI.v = 0*mV
neuronI.h = 1
neuronI.m = 0
neuronI.n = .5
neuronI.gNa = gNa0
"""
run(10 * msecond)

M = StateMonitor(neuronI, 'v', record=True)
S = SpikeMonitor(neuronI, variables='v')

print(S.t[:])
#import matplotlib.pyplot as plt
#plt.figure(figsize=(15, 5))
#plt.plot(S.t, S.v[0])
#plt.show()
