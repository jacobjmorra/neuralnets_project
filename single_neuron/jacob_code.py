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
