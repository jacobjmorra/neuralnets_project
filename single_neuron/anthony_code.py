# Auditory stimuli are modelled as step currents.
# Frequencies in range [.5,32]kHz.
# 100ms in length, 100ms offset --> 200ms rhythm.

# Next step: Need to configure the step input properly.
from brian2 import *

def pl_lif(pf = 5000*Hz, f = 10*Hz, runtime = 1000*ms, plot_spikes = True):
    start_scope()
    tau = 20*ms
    Vr, Vth, Vrest = -70*mV, -54*mV, -60*mV # Voltages reset/threshold/rest
    I, R = 10*nA, 1/(1000*nS) # Input current and resistance, not specified in paper, except that input is a step current.
    eqs = '''
    dv/dt = (-(v-Vrest) + I*R + A*sin(2*pi*f*t + b)) / tau : volt
    db/dt = abs(pf - f)/(31500*Hz*2*200*ms) : 1
    A = 10*mV : volt # Oscillation Amplitude, kept constant
    '''
    neuron = NeuronGroup(1, model=eqs, threshold='v>Vth', reset='v=Vr', refractory=5*ms, method='euler')
    neuron.v = Vrest; neuron.b = 0
    S = SpikeMonitor(neuron)
    trace = StateMonitor(neuron, 'v', record = True)

    run(runtime - (runtime - 500*ms))
    I = 0*nA; run(runtime - 500*ms)
    if plot_spikes:
        plot(trace.t/ms, trace.v.T); plot(S.t/ms, [Vth]*S.num_spikes, "ro"); xlabel('Time (ms)'); ylabel('Voltage (V)')
        tight_layout(); show()
    return(S.num_spikes / runtime)

pl_lif(pf = 5000*Hz, f = 10*Hz, runtime = 2000*ms)