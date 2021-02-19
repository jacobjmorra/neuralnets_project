# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html
from brian2 import *

def hh(duration = 2*second, gA = 10*nsiemens, gG = 20*nsiemens, plots = True):
    start_scope(); area = 20000*umetre**2 

    # Parameters
    Cm = 1*ufarad*cm**-2 * area
    gl, g_na, g_kd = 5e-5*siemens*cm**-2 * area, 100*msiemens*cm**-2 * area, 30*msiemens*cm**-2 * area
    El, EK, ENa = -65*mV, -90*mV, 50*mV
    VT = -63*mV
    w_e, w_i = gA, gG # Syanptic weights
    Ee, Ei = 0*mV, -80*mV # Synaptic potentials
    taue, taui = 5*ms, 10*ms # Synapse time constants

    eqs = Equations('''
    dv/dt = (g_na*(m**3)*h*(ENa-v) + g_kd*(n**4)*(EK-v) + gl*(El-v) + ge*(Ee-v) + gi*(Ei-v) + Iext)/Cm : volt
    dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
    dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
    dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
    dge/dt = -ge/taue : siemens
    dgi/dt = -gi/taui : siemens
    Iext : amp
    ''')

    group = NeuronGroup(3, eqs, threshold='v > -40*mV', refractory='v > -40*mV', method='exponential_euler')
    group.v = El; group.Iext = ".32*nA" # Initializing voltage and external current

    # Creating spike monitor, voltage monitors, and synapses
    spikes = SpikeMonitor(group)
    M = StateMonitor(group, 'v', record = 0) # Master
    S = StateMonitor(group, 'v', record = 1) # Slave
    I = StateMonitor(group, 'v', record = 2) # Inter
    MSSI = Synapses(group, group, 'w: siemens', on_pre='ge += w_e'); MSSI.connect(i = [0,1], j = [1, 2])
    IS = Synapses(group, group, 'w: siemens', on_pre='gi += w_i'); IS.connect(i = [2], j = [1])
    run(duration)
    # print("\nFiring Rates (M/S/I):\t\t%s" % spikes.count / duration)

    if plots:
        plot(M.t/ms, M.v[0], "k-", linewidth = .5)
        plot(S.t/ms, S.v[0], "r-", linewidth = .5, alpha = .8)
        plot(I.t/ms, I.v[0], "g--", linewidth = .5, alpha = .8)
        xlabel('Time (ms)'); ylabel('Voltage (mV)'); show()

    if spikes.count[0] == spikes.count[1]:
        tau_MS = spikes.spike_trains()[0] - spikes.spike_trains()[1]
        return(mean(tau_MS), std(tau_MS))
    else:
        return (float("NaN")*ms, float("NaN")*ms)

gA_max, gG_max = 15, 15
heat = zeros((gA_max, gG_max))
for i in range(gA_max):
    for j in range(gG_max):
        heat[i,j] = hh(gA = i*nsiemens, gG = j*nsiemens, plots = False)[0]/ms

heat_ma = ma.array(heat, mask = isnan(heat))
cmap = matplotlib.cm.seismic; cmap.set_bad('green', 1.)
imshow(heat_ma, interpolation='nearest', cmap=cmap, vmax = nanmax(abs(heat)), vmin = -nanmax(abs(heat)))
gca().invert_yaxis()
xlabel("AMPA synapse conductance (nS)"); ylabel("GABA synapsa conductance (nS)"); title("Tau")
colorbar(); show()