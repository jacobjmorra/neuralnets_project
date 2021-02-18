from brian2 import *

# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html

def hh_fire(duration = 2*second, plots = True):
    start_scope()
    area = 20000*umetre**2
    Cm = 9*pi*ufarad*cm**-2 * area
    gl = 2.7*pi*msiemens*cm**-2 * area
    g_na, g_kd = 1080*pi*msiemens*cm**-2 * area, 324*pi*msiemens*cm**-2 * area # Conductances
    El, EK, ENa = 10.6*mV, -12*mV, 115*mV # Potentials, El is resting
    VT = -63*mV

    eqs = Equations('''
    dv/dt = (g_na*(m**3)*h*(ENa-v) + g_kd*(n**4)*(EK-v) + gl*(El-v) + I)/Cm : volt
    dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
    dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
    dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
    I : amp
    ''')

    # Threshold and refractoriness are only used for spike counting
    group = NeuronGroup(3, eqs, threshold='v > -40*mV', refractory='v > -40*mV', method='exponential_euler')
    group.v = El; group.I = '.7*nA * i / 3'
    S = SpikeMonitor(group); M = StateMonitor(group, 'v', record = 0)
    run(duration)

    if plots:
        subplot(1,2,1); plot(group.I/nA, S.count / duration, "k.-", linewidth = .5) # Current vs. FR, 1 datum per neuron
        xlabel('Input Current (nA)'); ylabel('Firing rate (Hz)')
        subplot(1,2,2); plot(M.t/ms, M.v[0], "k-", linewidth = .5)
        xlabel('Time (ms)'); ylabel('Voltage (mV)')
        subplots_adjust(wspace=.5); show()

hh_fire(2*second)