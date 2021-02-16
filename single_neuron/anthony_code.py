# Auditory stimuli are modelled as step currents.
# Frequencies in range [.5,32]kHz.
# 100ms in length, 100ms offset, and 200ms rhythm.
from brian2 import *

def pl_lif(pf = 5000*Hz, f = 10*Hz, runtime = 1000*ms, plot_spikes = True):
    start_scope()
    tau = 20*ms
    Vr, Vth, Vrest = -70*mV, -54*mV, -60*mV # Voltages reset/threshold/rest
    I, R = 10*nA, 1/(1000*nS) # Input current and resistance, not specified in paper. 0 rn.
    eqs = '''
    dv/dt = (-(v-Vrest) + I*R + A*sin(2*pi*f*t + b*ms)) / tau : volt
    b = abs(pf - f)/(31500*Hz*2*200*ms) : Hz
    A = 10*mV : volt # Oscillation Amplitude, kept constant
    '''
    neuron = NeuronGroup(1, model=eqs, threshold='v>Vth', reset='v=Vr', refractory=5*ms, method='euler')
    neuron.v = Vrest
    S = SpikeMonitor(neuron)
    trace = StateMonitor(neuron, 'v', record = True)

    run(runtime)
    if plot_spikes:
        plot(trace.t/ms, trace.v.T); plot(S.t/ms, [Vth]*S.num_spikes, "ro"); xlabel('Time (ms)'); ylabel('Voltage (V)')
        tight_layout(); show()
    return(S.num_spikes / runtime)

pl_lif(pf = 5000*Hz, f = 6000*Hz, runtime = 1000*ms)