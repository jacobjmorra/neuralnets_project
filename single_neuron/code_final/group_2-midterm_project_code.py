"""
Group 2 members: Aymee Alvarez-Rivero (S/N: ), Alex Busch (S/N: ), Anthony Cruz (S/N: ), Jacob Morra (S/N: 251117805) 
"""
# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html

from brian2 import *

def plot_diff(num_incr=19):
    """Find differences between signals for varying gA from 1 ns to 20 ns, gG = 18"""
    import matplotlib.pyplot as plt

    gA_list = [1 + i for i in range(num_incr)] * nsiemens
    gG=18 * nsiemens

    #print(gA_list)

    diff_listMS = []
    diff_listSI = []
    diff_listMI = []
    total_diff_list = []

    for i in range(num_incr):
        #compare the 3 spike_train times for each gA value
        M_spiketrain,S_spiketrain,I_spiketrain = MSI(gA_list[i], gG)

        #check dims
        #print(M_spiketrain.size, S_spiketrain.size, I_spiketrain.size)

        max_len_to_consider = min(M_spiketrain.size, S_spiketrain.size, I_spiketrain.size)

        #if not all same, take smallest dimension from all
        diff_listMS.append(np.absolute(S_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))
        diff_listSI.append(np.absolute(I_spiketrain[0:max_len_to_consider]/ms - S_spiketrain[0:max_len_to_consider]/ms))
        diff_listMI.append(np.absolute(I_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))

        total_diff_list.append(np.absolute(S_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms)
                               + np.absolute(I_spiketrain[0:max_len_to_consider]/ms - S_spiketrain[0:max_len_to_consider]/ms)
                               + np.absolute(I_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))


    #MS diffs plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    for i in range(num_incr):
        x = [i for i in range(len(diff_listMS[i]))]
        ax1.plot(x, diff_listMS[i], label="gA="+str(gA_list[i]))
    ax1.legend(prop={'size': 6})
    ax1.set_title("MS")

    #SI diffs plot
    for i in range(num_incr):
        x = [i for i in range(len(diff_listSI[i]))]
        ax2.plot(x, diff_listSI[i], label="gA="+str(gA_list[i]))
    ax2.legend(prop={'size': 6})
    ax2.set_title("SI")

    #MI diffs plot
    for i in range(num_incr):
        x = [i for i in range(len(diff_listMI[i]))]
        ax3.plot(x, diff_listMI[i], label="gA="+str(gA_list[i]))
    ax3.legend(prop={'size': 6})
    ax3.set_title("MI")

    #Total diffs plot
    for i in range(num_incr):
        x = [i for i in range(len(total_diff_list[i]))]
        ax4.plot(x, total_diff_list[i], label="gA="+str(gA_list[i]))
    ax4.legend(prop={'size': 6})
    ax4.set_title("Differences (total)")

    plt.show()

def MSI(gA, gG):
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

    M_spikes = SpikeMonitor(group, record=True)
    S_spikes = SpikeMonitor(group, record=True)
    I_spikes = SpikeMonitor(group, record=True)

    run(500 * ms)
    # print("\nFiring Rates (M/S/I):\t\t%s" % spikes.count / duration)


    #plot(M.t/ms, M.v[0], "k-", linewidth = .5)
    #plot(S.t/ms, S.v[0], "r-", linewidth = .5, alpha = .8)
    #plot(I.t/ms, I.v[0], "g--", linewidth = .5, alpha = .8)
    #xlabel('Time (ms)'); ylabel('Voltage (mV)'); show()
    #plt.show()

    #print(get_diff(M.t, S.t))
    #print(M.t)
    #print(S.t)
    #print(I.t)

    #print("SPIKES,", M_spikes.spike_trains())

    #print(get_diff(M.t, I.t))

    return M_spikes.spike_trains()[0], M_spikes.spike_trains()[1], M_spikes.spike_trains()[2]
    #print(group.spikes)
    #return(M_spikes.spike_trains()[0], M.v[0], S.t, S.v[0], I.t, I.v[0])

#M_ts, S_ts, I_ts = MSI(gA = 10*nsiemens, gG = 20*nsiemens)

#print(M_ts)
#print(S_ts)
#print(I_ts)

plot_diff()
