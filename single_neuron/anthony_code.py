from brian2 import *

# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html

def hh_fire(duration = 2*second):
    # Parameters
    area = 20000*umetre**2
    Cm = 1*ufarad*cm**-2 * area
    gl = 5e-5*siemens*cm**-2 * area
    El, EK, ENa = -65*mV, -90*mV, 50*mV # Potentials
    g_na, g_kd = 100*msiemens*cm**-2 * area, 30*msiemens*cm**-2 * area # Conductances
    VT = -63*mV

    # The model
    eqs = Equations('''
    dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + I)/Cm : volt
    dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
    dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
    dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
    I : amp
    ''')

    # Threshold and refractoriness are only used for spike counting
    group = NeuronGroup(3, eqs,
                        threshold='v > -40*mV',
                        refractory='v > -40*mV',
                        method='exponential_euler')
    group.v = El
    group.I = '.7*nA * i / 3'

    monitor = SpikeMonitor(group)

    run(duration)
    plot(group.I/nA, monitor.count / duration)
    xlabel('I (nA)'); ylabel('Firing rate (Hz)'); show()
hh_fire()