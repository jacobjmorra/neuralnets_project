# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html
from brian2 import *

def tau_calc(M_spiketrain, S_spiketrain):
    if (len(M_spiketrain) == 0) or (len(S_spiketrain) == 0): return float("NaN")

    output = zeros(len(M_spiketrain))
    for i in range(len(M_spiketrain)):
        closest_S_spike = argmin(abs(array(S_spiketrain) - M_spiketrain[i]))
        output[i] = M_spiketrain[i] - S_spiketrain[closest_S_spike]
    return(output)

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
    M_end_spikes = array(spikes.spike_trains()[0]/ms)[array(spikes.spike_trains()[0]/ms) >= ((duration/ms) - 1000)]
    S_end_spikes = array(spikes.spike_trains()[1]/ms)[array(spikes.spike_trains()[1]/ms) >= ((duration/ms) - 1000)]

    
    return tau_calc(M_end_spikes, S_end_spikes)

def heatmap(gA_max, gG_max, step_size):
    if gA_max % step_size != 0 or gG_max % step_size != 0: return(float("nan"))

    heat = zeros((int(gG_max/step_size)+1, int(gA_max/step_size)+1))
    for j in range(int(gA_max/step_size)+1):
        for i in range(int(gG_max/step_size)+1):
            heat[i,j] = mean(hh(gA = step_size*j*nsiemens, gG = step_size*i*nsiemens, plots = False))

    heat_ma = ma.array(heat, mask = isnan(heat))
    cmap = matplotlib.cm.seismic; cmap.set_bad('green', 1.)
    imshow(heat_ma, interpolation='nearest', cmap=cmap, vmax = nanmax(abs(heat)), vmin = -nanmax(abs(heat)))
    gca().invert_yaxis()
    xticks(list(range(int(gA_max/step_size)+1)), [step_size*i for i in (range(int(gA_max/step_size)+1))])
    yticks(list(range(int(gG_max/step_size)+1)), [step_size*i for i in (range(int(gG_max/step_size)+1))])
    ylabel("GABA synapse conductance (nS)"); xlabel("AMPA synapse conductance (nS)"); title("\u03C4 (ms)")
    colorbar(); show()

gA_tau_means = zeros(100)
for i in range(100):
    print(i)
    gA_tau_results = hh(gA = i*nsiemens, plots=False)
    if not isnan(any(gA_tau_results)):
        gA_tau_means[i] = mean(gA_tau_results)
        if not isinstance(gA_tau_results, float):
            for j in gA_tau_results: plot(i, j, "b.")
        else:
            plot(i, gA_tau_results, "b.")
    else:
        gA_tau_means[i] = float("NaN")

mask = isfinite(gA_tau_means)
x_values, y_values = [], []
for i in range(100):
    if mask[i] == True:
        x_values.append(i)
        y_values.append(gA_tau_means[i])
plot(x_values, y_values, "k-", linewidth = 3)
xlabel("AMPA Synapse Conductance (nS)"); ylabel("\u03C4 (ms)")
show()

# heatmap(10, 10, 2)