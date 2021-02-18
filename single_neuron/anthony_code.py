# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html
from brian2 import *

def hh(duration = 2*second, plots = True):
    start_scope()
    Cm = 9*pi*ufarad
    gl = 2.7*pi*msiemens
    g_na, g_kd = 1080*pi*msiemens, 324*pi*msiemens # Conductances
    El, EK, ENa = 10.6*mV, -12*mV, 115*mV # Potentials, El is resting
    VT = -63*mV

    eqs = Equations('''
    dv/dt = (g_na*(m**3)*h*(ENa-v) + g_kd*(n**4)*(EK-v) + gl*(El-v) + I)/Cm : volt
    dm/dt = (Am*v*(1-m) - Bm*v*m)/volt: 1
    dn/dt = (An*v*(1-m) - Bn*v*m)/volt : 1
    dh/dt = (Ah*v*(1-m) - Bh*v*m)/volt : 1

    An = ((10*mV - v)/volt) / ((100*exp( ((10*mV - v)/10)/volt )-1)*ms) : Hz
    Bn = .125 * exp(-v/(80*volt)) / ms : Hz
    Am = ((25*mV - v)/volt) / ((10*exp( ((25*mV - v)/10)/volt )-1)*ms) : Hz
    Bm = 4 * exp(-v/(18*volt)) / ms : Hz
    Ah = .07 * exp(-v/(20*volt)) / ms : Hz
    Bh = 1 / ((exp( ((30*mV - v)/10)/volt )+1)*ms) : Hz
    I : amp
    ''')

    group = NeuronGroup(3, eqs, threshold='v > -40*mV', refractory='v > -40*mV', method='exponential_euler')
    group.v = El; group.I = '.7*nA * i/3'
    spikes = SpikeMonitor(group); M = StateMonitor(group, 'v', record = 0)
    S = StateMonitor(group, 'v', record = 1); I_m = StateMonitor(group, 'v', record = 2)
    run(duration)

    if plots:
        subplot(1,2,1); plot(group.I/nA, spikes.count / duration, "k.-", linewidth = .5) # Current vs. FR, 1 datum per neuron
        xlabel('Input Current (nA)'); ylabel('Firing rate (Hz)')
        subplot(1,2,2); plot(M.t/ms, M.v[0], "g--", linewidth = .5)
        plot(S.t/ms, S.v[0], "r--", linewidth = .5); plot(I_m.t/ms, I_m.v[0], "b--", linewidth = .5)
        xlabel('Time (ms)'); ylabel('Voltage (mV)')
        subplots_adjust(wspace=.5); show()
hh(10*second)