"""
Midterm Project Code for APPLMATH 9624B/4264B (Introduction to Neural Networks); Instructors: Dr. Lyle Muller, Dr. Marieke Mur
Group 2 members: Aymee Alvarez-Rivero (S/N: ), Alex Busch (S/N: ), Anthony Cruz (S/N: ), Jacob Morra (S/N: 251117805)

Brief Description:
The following code illustrates our implementation of [1]; specifically, we highlight the phenomenon that "two identical
autonomous systems in a master-slave configuration can have Anticipated Synchronization; that is, wherein the slave's
spiking switches from following the master's spikes (i.e. "Delayed Synchronization [DS]) to preceding them. [1]

Selected Paper:
[1] F. S. Matias, P. V. Carelli, C. R. Mirasso, and M. Copelli, “Anticipated synchronization in a biologically plausible model of
neuronal motifs,” Phys. Rev. E Stat. Nonlin. Soft Matter Phys., vol. 84, no. 2 Pt 1, p. 021922, Aug. 2011.
"""

from brian2 import *

def MSI(gA, gG):
    """
    This function initializes and runs a brian2 simulation which is based off of brian2's Hodgkin-Huxley model (please see
    https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html); from here 3 neurons and their synapses
    are created (1 master, 1 slave, 1 interneuron) --- with excitatory conductance gA and inhibitory conductance gG. A
    constant current of 320 pA is initialized, the simulation is run for 500 ms, and all spike trains are returned.
    :param gA: excitatory conductance (allocated as the excitatory weight)
    :param gG: inhibitory conductance (allocated as the inhibitory weight)
    :return: spiketrains for M, S, and I neurons in the MSI configuration
    """

    start_scope()

    area = 20000*umetre**2

    # Parameters
    Cm = 1*ufarad*cm**-2 * area
    gl, g_na, g_kd = 5e-5*siemens*cm**-2 * area, 100*msiemens*cm**-2 * area, 30*msiemens*cm**-2 * area
    El, EK, ENa = -65*mV, -90*mV, 50*mV
    VT = -63*mV
    w_e, w_i = gA, gG # Syanptic weights
    Ee, Ei = 0*mV, -80*mV # Synaptic potentials
    taue, taui = 5*ms, 10*ms # Synapse time constants

    #from the brian2 HH example: # https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html
    eqs = Equations('''
    dv/dt = (g_na*(m**3)*h*(ENa-v) + g_kd*(n**4)*(EK-v) + gl*(El-v) + ge*(Ee-v) + gi*(Ei-v) + Iext)/Cm : volt
    dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
    dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
    dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
    dge/dt = -ge/taue : siemens
    dgi/dt = -gi/taui : siemens
    Iext : amp
    ''')

    #Create the MIS neurons and initialize voltage and external current parameters
    group = NeuronGroup(3, eqs, threshold='v > -40*mV', refractory='v > -40*mV', method='exponential_euler')
    group.v = El; group.Iext = ".32*nA"

    # Create state monitors, spike monitors, synapses
    M = StateMonitor(group, 'v', record = 0) # Master
    S = StateMonitor(group, 'v', record = 1) # Slave
    I = StateMonitor(group, 'v', record = 2) # Inter

    MSSI = Synapses(group, group, 'w: siemens', on_pre='ge += w_e'); MSSI.connect(i = [0,1], j = [1, 2])
    IS = Synapses(group, group, 'w: siemens', on_pre='gi += w_i'); IS.connect(i = [2], j = [1])

    MSI_spikes = SpikeMonitor(group, record=True)

    #run simulation for .5 seconds
    run(500 * ms)

    #print("\nFiring Rates (M/S/I):\t\t%s" % spikes.count / duration)


    #Uncomment below to plot M,I,S potentials over run time
    """
    #import matplotlib.pyplot as plt
    #plot(M.t/ms, M.v[0], "k-", linewidth = .5)
    #plot(S.t/ms, S.v[0], "r-", linewidth = .5, alpha = .8)
    #plot(I.t/ms, I.v[0], "g--", linewidth = .5, alpha = .8)
    #xlabel('Time (ms)'); ylabel('Voltage (mV)'); show()
    #plt.show()
    """

    return MSI_spikes.spike_trains()[0], MSI_spikes.spike_trains()[1], MSI_spikes.spike_trains()[2]

def heatmap(gA_max, gG_max, step_size):
    """
    This function creates a heat map ... ( *** will finish this today)
    :param gA_max:
    :param gG_max:
    :param step_size:
    :return:
    """
    if gA_max % step_size != 0 or gG_max % step_size != 0: return(float("nan"))

    heat = zeros((int(gG_max/step_size)+1, int(gA_max/step_size)+1))
    for j in range(int(gA_max/step_size)+1):
        for i in range(int(gG_max/step_size)+1):
            heat[i,j] = mean(MSI(gA = step_size*j*nsiemens, gG = step_size*i*nsiemens))

    heat_ma = ma.array(heat, mask = isnan(heat))
    cmap = matplotlib.cm.seismic; cmap.set_bad('green', 1.)
    imshow(heat_ma, interpolation='nearest', cmap=cmap, vmax = nanmax(abs(heat)), vmin = -nanmax(abs(heat)))
    gca().invert_yaxis()
    xticks(list(range(int(gA_max/step_size)+1)), [step_size*i for i in (range(int(gA_max/step_size)+1))])
    yticks(list(range(int(gG_max/step_size)+1)), [step_size*i for i in (range(int(gG_max/step_size)+1))])
    ylabel("GABA synapse conductance (nS)"); xlabel("AMPA synapse conductance (nS)"); title("\u03C4 (ms)")
    colorbar(); show()

def plot_diff(num_incr=19):
    """
    The function calculates the differences in ISI spike times between signals for varying gA from 1 ns to 20 ns, gG = 18.
    The differences between each of M, S, and I spiketrains are plotted, as well as the overall differences.
    :param num_incr: the number of desired steps of 1 * ns, starting from gA = 1 * ns. i.e. num_incr = 19 would indicate
    19 steps of 1, and therefore gA would range from 1 * ns to 20 * ns.
    :return: None
    """
    import matplotlib.pyplot as plt

    #store all of the gA values, gG is fixed at 18 ns
    gA_list = [1 + i for i in range(num_incr)] * nsiemens
    gG=18 * nsiemens

    #store all of the ISI differences between M, S, I, and all neurons
    diff_listMS = []
    diff_listSI = []
    diff_listMI = []
    total_diff_list = []

    #for each gA value, compare the M, S, I, and all spiketrains
    for i in range(num_incr):
        #compare the 3 spike_train times for each gA value
        M_spiketrain,S_spiketrain,I_spiketrain = MSI(gA_list[i], gG)

        #check dims
        #print(M_spiketrain.size, S_spiketrain.size, I_spiketrain.size)

        #for dimensions to work, truncate all vectors to the lowest common size
        max_len_to_consider = min(M_spiketrain.size, S_spiketrain.size, I_spiketrain.size)

        #add differences for particular gA value to respective lists
        diff_listMS.append(np.absolute(S_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))
        diff_listSI.append(np.absolute(I_spiketrain[0:max_len_to_consider]/ms - S_spiketrain[0:max_len_to_consider]/ms))
        diff_listMI.append(np.absolute(I_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))
        total_diff_list.append(np.absolute(S_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms)
                               + np.absolute(I_spiketrain[0:max_len_to_consider]/ms - S_spiketrain[0:max_len_to_consider]/ms)
                               + np.absolute(I_spiketrain[0:max_len_to_consider]/ms - M_spiketrain[0:max_len_to_consider]/ms))

    #plot the differences in ISIs

    #MS differences plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    for i in range(num_incr):
        x = [i for i in range(len(diff_listMS[i]))]
        ax1.plot(x, diff_listMS[i], label="gA="+str(gA_list[i]))
    ax1.legend(prop={'size': 6})
    ax1.set_title("MS")

    #SI differences plot
    for i in range(num_incr):
        x = [i for i in range(len(diff_listSI[i]))]
        ax2.plot(x, diff_listSI[i], label="gA="+str(gA_list[i]))
    ax2.legend(prop={'size': 6})
    ax2.set_title("SI")

    #MI differences plot
    for i in range(num_incr):
        x = [i for i in range(len(diff_listMI[i]))]
        ax3.plot(x, diff_listMI[i], label="gA="+str(gA_list[i]))
    ax3.legend(prop={'size': 6})
    ax3.set_title("MI")

    #Total differences plot
    for i in range(num_incr):
        x = [i for i in range(len(total_diff_list[i]))]
        ax4.plot(x, total_diff_list[i], label="gA="+str(gA_list[i]))
    ax4.legend(prop={'size': 6})
    ax4.set_title("Differences (total)")

    plt.show()

def plot_S_I_diffs(gA = 20 * nsiemens):
    """
    Compare the spiketrain times for current and next spiketrain times on the M neuron from the S neuron at time t
    :return: None
    """
    import matplotlib.pyplot as plt

    #initialize parameters
    gG = 18 * nsiemens

    # store all of the ISI differences between S neuron and M at time t and time t+1
    diff_S_to_M_t = []
    diff_S_to_M_t_next = []

    #get the spiketrains
    M, S, _ = MSI(gA, gG)

    #is S closer to t, or t+1?
    #print(M)
    #print(S)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # for gA = 10
    print("gA = 10 ********************************************************************************")
    gA = 10
    M, S, _ = MSI(10 * nsiemens, gG)
    print(M)
    print(S)
    for i in range(min(len(M), len(S))-1):
        #print(i)
        diff_S_to_M_t.append(abs(M[i]-S[i]))
        diff_S_to_M_t_next.append(abs((S[i+1]-M[i])))
    print(diff_S_to_M_t)
    print(diff_S_to_M_t_next)
    ax1.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="S[t] - M[t]")
    ax1.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="S[t+1] - M[t]")
    ax1.legend()
    ax1.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 15
    print("gA = 15 ********************************************************************************")
    gA = 15
    diff_S_to_M_t = []; diff_S_to_M_t_next = []
    M, S, _ = MSI(15 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(abs(M[i]-S[i]))
        diff_S_to_M_t_next.append(abs((S[i+1]-M[i])))
    ax2.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="S[t] - M[t]")
    ax2.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="S[t+1] - M[t]")
    ax2.legend()
    ax2.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 20
    print("gA = 20 ********************************************************************************")
    gA = 20
    M, S, _ = MSI(20 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(abs(M[i]-S[i]))
        diff_S_to_M_t_next.append(abs((S[i+1]-M[i])))
    ax3.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="S[t] - M[t]")
    ax3.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="S[t+1] - M[t]")
    ax3.legend()
    ax3.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 25
    print("gA = 25 ********************************************************************************")
    gA = 25
    M, S, _ = MSI(25 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(abs(M[i]-S[i]))
        diff_S_to_M_t_next.append(abs((S[i+1]-M[i])))
    ax4.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="S[t] - M[t]")
    ax4.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="S[t+1] - M[t]")
    ax4.legend()
    ax4.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))
    plt.show()

plot_S_I_diffs()