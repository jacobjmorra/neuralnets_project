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

#Sender = M
#Receiver = S
#I = Interneuron

from brian2 import *

def check_for_candidate(tau):
    for i in range(len(tau) - 10):
        #want to check for 5 steps below 0, 5 steps above 0 (arbitrary)
        if(sum(np.sign(tau[i:i+5])) == -5):
            #check the rest of tau for 5 steps above 0
            for j in range(i + 5, len(tau) - 5):
                if(sum(np.sign(tau[j:j+5])) == 5):
                    return True
    return False

def tau_calc2(M_spiketrain, S_spiketrain):
    """
    Given M and S spiketrains, return tau (i.e. all time differences between the M spiketrain and
    all closest S spiketrain times).
    :param M_spiketrain: an array representing the spiketrain of the Master neuron
    :param S_spiketrain: an array representing the spiketrain of the Slave neuron
    :return: tau
    """
    tau = zeros(len(M_spiketrain))

    #print(M_spiketrain)
    #print(S_spiketrain)
    for i in range(len(M_spiketrain)):
        closest_index = argmin(abs(S_spiketrain - M_spiketrain[i]))
        #print("closest index is ", closest_index)
        tau[i] = M_spiketrain[i] - S_spiketrain[closest_index]
    return tau


def tau_calc(M_spiketrain, S_spiketrain):
    """
    Given a spiketrain for the Master neuron (M_spiketrain) and slave neuron (S_spiketrain), determine the difference
    in spike time between each slave spike and its closest master spike. The resulting array is returned.
    :param M_spiketrain: an array of spike times; the spiketrain for the master neuron
    :param S_spiketrain: an array of spike times; the spiketrain for the slave neuron
    :return: an array of time differences between each master spike time and its closest slave spike time
    """
    # if no spikes, return NaN
    if (len(M_spiketrain) == 0) or (len(S_spiketrain) == 0): return float("NaN")

    #otherwise, create array of 0's, fill with time differences between each master spike and its closest slave spike
    output = zeros(len(M_spiketrain))
    for i in range(len(M_spiketrain)):
        closest_S_spike = argmin(abs(array(S_spiketrain) - M_spiketrain[i]))
        output[i] = M_spiketrain[i] - S_spiketrain[closest_S_spike]
    return(output)


def MSI2(gA, gG):
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

    #specify area of region (necessary for initialization of parameters Cm, gl, g_na, and g_k)
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
    run(4000 * ms)

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

    M_end_spikes = array(MSI_spikes.spike_trains()[0] / ms)[array(MSI_spikes.spike_trains()[0] / ms) >= ((500*ms / ms) - 250)]
    S_end_spikes = array(MSI_spikes.spike_trains()[1] / ms)[array(MSI_spikes.spike_trains()[1] / ms) >= ((500*ms / ms) - 250)]

    return MSI_spikes.spike_trains()[0], MSI_spikes.spike_trains()[1], MSI_spikes.spike_trains()[2], tau_calc(M_end_spikes, S_end_spikes)

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

    #specify area of region (necessary for initialization of parameters Cm, gl, g_na, and g_k)
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

    M_end_spikes = array(MSI_spikes.spike_trains()[0] / ms)[array(MSI_spikes.spike_trains()[0] / ms) >= ((500*ms / ms) - 250)]
    S_end_spikes = array(MSI_spikes.spike_trains()[1] / ms)[array(MSI_spikes.spike_trains()[1] / ms) >= ((500*ms / ms) - 250)]

    return MSI_spikes.spike_trains()[0], MSI_spikes.spike_trains()[1], MSI_spikes.spike_trains()[2], tau_calc(M_end_spikes, S_end_spikes)

def SRI():
    """
    This function runs an alternate simulation of the MSI model. The SRI or "Sender-Receiver-Interneuron" model runs
    for 1000 * ms and uses a tuned combination of parameters from the brian2 Hodgkin-Huxley example and [1]. The
    spiketrains for the S, R, and I neurons are plotted. No values are returned.
    :return: None
    """
    start_scope()

    num_neurons = 3
    duration = 1 * second

    # Parameters
    area = 20000 * umetre ** 2
    Cm = 1 * ufarad * cm ** -2 * area
    gl = 5e-5 * siemens * cm ** -2 * area
    El = -65 * mV
    EK = -90 * mV
    ENa = 50 * mV
    g_na = 100 * msiemens * cm ** -2 * area
    g_kd = 30 * msiemens * cm ** -2 * area
    VT = -63 * mV

    # Constant current, determines neuron excitability
    Iext_S = 28e-11 * amp
    Iext_I = 28e-11 * amp
    Iext_R = 28e-11 * amp  # variable
    Iext_values = [Iext_S, Iext_I, Iext_R]

    w_e = 1 * nsiemens
    w_i = 10 * nsiemens  # 60*nsiemens 40*nsiemens #20*nsiemens
    Ee = 0 * mV
    Ei = -80 * mV
    taue = 5 * ms
    taui = 10 * ms

    # The model
    eqs = Equations('''
    dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + Iext + g_e*(Ee-v) + g_i*(Ei-v))/Cm : volt
    dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
    dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
    dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
    dg_e/dt = -g_e/taue : siemens
    dg_i/dt = -g_i/taui : siemens
    Iext : amp
    ''')

    # Threshold and refractoriness are only used for spike counting
    group = NeuronGroup(num_neurons, eqs,
                        threshold='v > -40*mV',
                        refractory='v > -40*mV',
                        method='exponential_euler')
    group.v = El
    group.Iext = Iext_values
    # group.I = '0.7*nA * i / num_neurons'

    A = Synapses(group, group, 'w: siemens', on_pre='g_e += w_e')

    # Excitatory synapses from Sender to Reveiver (0 to 1)
    # and Receiver to Interneuron (1 to 2)
    A.connect(i=[0, 1], j=[1, 2])

    # Inhibitory synapse
    G = Synapses(group, group, 'w: siemens', on_pre='g_i += w_i')

    # # Inhibitory synapses from Interneuron to Reveiver (2 to 1)
    G.connect(i=[2], j=[1])

    # Monitor variables
    monitorV = StateMonitor(group, 'v', record=True)
    monitorS = SpikeMonitor(group)

    run(duration)

    plt.plot()
    suptitle('SRI Neurons', size=20)
    title('Excitatory Weight = 1nS, Inhibitory Weight = 10nS', size=30)
    plt.plot(monitorV.v[0], c='blue', label='Neuron S')
    plt.plot(monitorV.v[1], c='red', label='Neuron R')
    plt.plot(monitorV.v[2], c='green', label='Neuron I')
    xlabel('Time', size=20), ylabel('Membrane Potential', size=20)
    legend()
    xticks(size=20), yticks(size=20)
    plt.show()


def heatmap(gA_max, gG_max, step_size):
    """
    This function creates and plots a heat map of mean AMPA vs GABA synapse conductances based on
    input gA and gG values in a given range based on the specified step_size parameter.
    :param gA_max: the max AMPA synaptic conductance to test
    :param gG_max: the max GABA synaptic conductance to test
    :param step_size: specifies the number of gA and gG values to test (i.e. 100 if step_size = 10)
    :return: None
    """
    if gA_max % step_size != 0 or gG_max % step_size != 0: return(float("nan"))

    heat = zeros((int(gG_max/step_size)+1, int(gA_max/step_size)+1))
    for j in range(int(gA_max/step_size)+1):
        print(str(100* (j / (int(gA_max/step_size)+1))) + " % done ...")
        for i in range(int(gG_max/step_size)+1):
            #store the mean value into a variable
            _,_,_,mean_value = MSI(gA = step_size*j*nsiemens, gG = step_size*i*nsiemens)
            #set each mean value as a heat map entry
            heat[i,j] = mean(mean_value)

    heat_ma = ma.array(heat, mask = isnan(heat))
    cmap = matplotlib.cm.seismic; cmap.set_bad('green', 1.)
    imshow(heat_ma, interpolation='nearest', cmap=cmap, vmax = nanmax(abs(heat)), vmin = -nanmax(abs(heat)))
    gca().invert_yaxis()
    xticks(list(range(int(gA_max/step_size)+1)), [step_size*i for i in (range(int(gA_max/step_size)+1))])
    yticks(list(range(int(gG_max/step_size)+1)), [step_size*i for i in (range(int(gG_max/step_size)+1))])
    ylabel("GABA synapse conductance (nS)"); xlabel("AMPA synapse conductance (nS)"); title("\u03C4 (ms)")
    colorbar(); show()

def plot_synapse_conductances():
    """
    The function plots average AMPA synaptic conductances vs. time in ms. It tests for a range
    of gA values (specified within the for loop range) and runs a simulation for each, returning
    (and ultimately plotting) the mean spike time differences (i.e. tau).
    :return: None
    """

    gA = 10 * nsiemens
    gA_tau_means = zeros(100)
    for i in range(100):
        print(str(100* (i / (len(gA_tau_means)))) + " % done...")
        _,_,_,MSImean = MSI(gA + i * nsiemens, 18 * nsiemens)
        gA_tau_results = MSImean
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
    for i in range(len(gA_tau_means)):
        if mask[i] == True:
            x_values.append(i)
            y_values.append(gA_tau_means[i])
    plot(x_values, y_values, "k-", linewidth=3)
    xlabel("AMPA Synapse Conductance (nS)");
    ylabel("\u03C4 (ms)")
    show()

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
        M_spiketrain,S_spiketrain,I_spiketrain,_ = MSI(gA_list[i], gG)

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
    M, S, _, _ = MSI(gA, gG)

    #is S closer to t, or t+1?
    #print(M)
    #print(S)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # for gA = 10
    print("gA = 10 ********************************************************************************")
    gA = 10
    M, S, _, _= MSI(10 * nsiemens, gG)
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
    M, S, _, _ = MSI(15 * nsiemens, gG)
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
    M, S, _, _ = MSI(20 * nsiemens, gG)
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
    M, S, _, _ = MSI(25 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(abs(M[i]-S[i]))
        diff_S_to_M_t_next.append(abs((S[i+1]-M[i])))
    ax4.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="S[t] - M[t]")
    ax4.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="S[t+1] - M[t]")
    ax4.legend()
    ax4.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))
    plt.show()

def plot_S_I_diffs2(gA = 20 * nsiemens):
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
    #diff_S_to_M_t_best = []

    #get the spiketrains
    M, S, _, _ = MSI(gA, gG)

    #is S closer to t, or t+1?
    #print(M)
    #print(S)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # for gA = 10
    print("gA = 10 ********************************************************************************")
    gA = 10
    M, S, _, _= MSI(10 * nsiemens, gG)
    print(M)
    print(S)
    for i in range(min(len(M), len(S))-1):
        #print(i)
        diff_S_to_M_t.append(M[i]-S[i])
        diff_S_to_M_t_next.append(M[i]- S[i+1])
    diff_S_to_M_t_best = tau_calc2(M,S)

    print(diff_S_to_M_t)
    print(diff_S_to_M_t_next)
    ax1.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="M[t] - S[t]")
    ax1.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="M[t] - S[t+1]")
    ax1.plot([i for i in range(len(diff_S_to_M_t_best))], diff_S_to_M_t_best, label="M[t] - S[t_closest]")
    ax1.legend()
    ax1.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 15
    print("gA = 15 ********************************************************************************")
    gA = 15
    diff_S_to_M_t = []; diff_S_to_M_t_next = []
    M, S, _, _ = MSI(15 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(M[i]-S[i])
        diff_S_to_M_t_next.append(M[i]-S[i+1])
    diff_S_to_M_t_best = tau_calc2(M, S)
    ax2.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="M[t] - S[t]")
    ax2.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="M[t] - S[t+1]")
    ax2.plot([i for i in range(len(diff_S_to_M_t_best))], diff_S_to_M_t_best, label="M[t] - S[t_closest]")
    ax2.legend()
    ax2.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 20
    print("gA = 20 ********************************************************************************")
    gA = 20
    M, S, _, _ = MSI(20 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(M[i]-S[i])
        diff_S_to_M_t_next.append(M[i]-S[i+1])
    diff_S_to_M_t_best = tau_calc2(M, S)
    ax3.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="M[t] - S[t]")
    ax3.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="M[t] - S[t+1]")
    ax3.plot([i for i in range(len(diff_S_to_M_t_best))], diff_S_to_M_t_best, label="M[t] - S[t_closest]")
    ax3.legend()
    ax3.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))

    # for gA = 25
    print("gA = 25 ********************************************************************************")
    gA = 25
    M, S, _, _ = MSI(25 * nsiemens, gG)
    for i in range(min(len(M), len(S))-1):
        diff_S_to_M_t.append(M[i]-S[i])
        diff_S_to_M_t_next.append(M[i]-S[i+1])
    diff_S_to_M_t_best = tau_calc2(M, S)
    ax4.plot([i for i in range(len(diff_S_to_M_t))], diff_S_to_M_t, label="M[t] - S[t]")
    ax4.plot([i for i in range(len(diff_S_to_M_t_next))], diff_S_to_M_t_next, label="M[t] - S[t+1]")
    ax4.plot([i for i in range(len(diff_S_to_M_t_best))], diff_S_to_M_t_best, label="M[t] - S[t_closest]")
    ax4.legend()
    ax4.set_title("Differences between S[t+1] and M[t], S[t] and M[t], gA = " + str(gA))
    plt.show()


def plot_tau_diff():
    #import matplotlib
    #matplotlib.rcParams['axes.facecolor'] = "white"
    #matplotlib.rcParams['figure.facecolor'] = "white"
    #matplotlib.rcParams['lines.linewidth'] = 1.1
    #import matplotlib.pyplot as plt
    #plt.style.use('fivethirtyeight')


    gA = 10 * nsiemens
    gA_list = [1 * nsiemens + i * nsiemens for i in range(20)]
    gG = 10 * nsiemens
    gG_list = [1 * nsiemens + i * nsiemens for i in range(20)]
    gA_rand = [gA_list[np.random.randint(0,len(gA_list)-1)] for i in range(len(gA_list))]
    gG_rand = [gG_list[np.random.randint(0,len(gG_list)-1)] for i in range(len(gG_list))]

    good_taus_gA = []
    good_taus_gG = []
    good_taus_rnd = []
    #M, S, _, _ = MSI(gA, gG)

    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2)

    for i in range(len(gA_list)):
        print("i1 = ", i)
        M, S, _, _ = MSI(gA_list[i], gG)
        #print(tau_calc2(M,S))
        if(tau_calc2(M,S)[1] < 0 and tau_calc2(M,S)[len(tau_calc2(M,S))-1]>0):
            #we may have a candidate for DS -> AS
            print("HAHA")
            good_taus_gA.append([gA_list[i], gG, tau_calc2(M,S)])
        ax1.plot([i for i in range(len(M))], tau_calc2(M,S), label= "gA = " + str(gA_list[i]) + ", gG = " + str(gG))
    #ax1.legend(prop={'size': 4})
    ax1.set_title(r'$\tau$' + " for gA = " + str(gA_list[0]) + " to " + str(gA_list[len(gA_list)-1]), size=13)
    #ax1.set_xlabel(r"$t_i$")
    ax1.set_ylabel(r"$\tau$")

    for i in range(len(gG_list)):
        print("i2 = ", i)
        M, S, _, _ = MSI(gA, gG_list[i])
        if(tau_calc2(M,S)[1]<0 and tau_calc2(M,S)[len(tau_calc2(M,S))-1]>0):
            # we may have a candidate for DS -> AS
            good_taus_gG.append([gA, gG_list[i], tau_calc2(M,S)])
            print("YO!!")
        ax2.plot([i for i in range(len(M))], tau_calc2(M,S), label= "gA = " + str(gA) + ", gG = " + str(gG_list[i]))
    #ax2.legend(prop={'size': 4})
    ax2.set_title(r'$\tau$' + " for gG = " + str(gG_list[0]) + " to " + str(gG_list[len(gG_list) - 1]), size=13)
    #ax2.set_xlabel(r"$t_i$")
    #ax2.set_ylabel(r"$\tau$")

    for i in range(len(gA_rand)):
        print("i3 = ", i)
        #pick random gA, gG pair
        M, S, _, _ = MSI(gA_rand[i], gG_rand[i])
        if(tau_calc2(M,S)[1] < 0 and tau_calc2(M,S)[len(tau_calc2(M,S))-1]>0):
            #we may have a candidate for DS -> AS
            print("HE!!!")
            good_taus_rnd.append([gA_rand[i], gG_rand[i], tau_calc2(M,S)])
        ax3.plot([i for i in range(len(M))], tau_calc2(M,S))
    ax3.set_title(r'$\tau$' + " for randomly-selected gA, gG", size=13)
    ax3.set_xlabel(r"$t_i$")
    ax3.set_ylabel(r"$\tau$")

    for i in range(len(good_taus_rnd)):
        if(len(good_taus_rnd)>0):
            #print([i for i in range(len(good_taus_rnd[i]))])
            #print(len([i for i in range(len(good_taus_rnd[i]))]))
            #print(good_taus_rnd[i])
            #print(len(good_taus_rnd[i]))
            ax4.plot([i for i in range(len(good_taus_rnd[i][2]))], good_taus_rnd[i][2], label= "gA=" + str(good_taus_rnd[i][0]) + ", gG=" + str(good_taus_rnd[i][1]))
    for i in range(len(good_taus_gA)):
        ax4.plot([i for i in range(len(good_taus_gA[i][2]))], good_taus_gA[i][2], label= "gA=" + str(good_taus_gA[i][0]) + ", gG=" + str(good_taus_gA[i][1]))
    for i in range(len(good_taus_gG)):
        ax4.plot([i for i in range(len(good_taus_gG[i][2]))], good_taus_gG[i][2], label= "gA=" + str(good_taus_gG[i][0]) + ", gG=" + str(good_taus_gG[i][1]))
    ax4.legend(prop={'size': 6})
    ax4.set_title(r'$\tau$' + " for gA and gG candidate values (DS to AS)", size=13)

    ax4.set_xlabel(r"$t_i$")
    #ax4.set_ylabel(r"$\tau$")
    plt.show()

    print("!!!!!!!!!!!!", good_taus_rnd)
    print("!!!!!!!!!!!!", good_taus_gA)
    print("!!!!!!!!!!!!", good_taus_gG)


def plot_tau_diff2():
    #import matplotlib
    #matplotlib.rcParams['axes.facecolor'] = "white"
    #matplotlib.rcParams['figure.facecolor'] = "white"
    #matplotlib.rcParams['lines.linewidth'] = 1.1
    #import matplotlib.pyplot as plt
    #plt.style.use('fivethirtyeight')


    gA_list = [0 * nsiemens + 2*i * nsiemens for i in range(20)]
    gG_list = [0 * nsiemens + 2*i * nsiemens for i in range(40)]
    gA_rand = [gA_list[np.random.randint(0,len(gA_list)-1)] for i in range(len(gA_list))]
    gG_rand = [gG_list[np.random.randint(0,len(gG_list)-1)] for i in range(len(gG_list))]

    good_taus = []

    for i in range(len(gA_list)):
        for j in range(len(gG_list)):
            print(i, ", ", j)
            M, S, _, _ = MSI(gA_list[i], gG_list[j])

            if(check_for_candidate(tau_calc2(M,S))):
                good_taus.append([gA_list[i], gG_list[j], tau_calc2(M,S)])
                print("gA = ", gA_list[i], "; gG = ", gG_list[j])

    fig, ax = plt.subplots()

    if(len(good_taus) > 1):
        for i in range(len(good_taus)):
            ax.plot([i for i in range(len(good_taus[i][2]))], good_taus[i][2], label="gA=" + str(good_taus[i][0]) + ", gG=" + str(good_taus[i][1]))
            ax.legend(prop={'size': 6})
            ax.set_title(r'$\tau$' + " for gA and gG candidate values (DS to AS)", size=15)
            ax.set_xlabel(r"$t_i$")
            ax.set_ylabel(r"$\tau$")
        plt.show()

    return good_taus
    """
    for i in range(len(gA_rand)):
        print(i)
        M, S, _, _ = MSI(gA_rand[i], gG_rand[i])

        if(check_for_candidate(tau_calc2(M,S))):
            good_taus.append([gA_rand[i], gG_rand[i], tau_calc2(M,S)])
            print("oooo")
    """

def plot_periodicity():
    import pickle
    with open('candidate_params.txt', 'rb') as f:
        good_taus = pickle.load(f)

    print(good_taus)

    fig, ax = plt.subplots()

    for i in range(len(good_taus)):
        ax.plot([i for i in range(len(good_taus[i][2]))], good_taus[i][2],
                    label="gA=" + str(good_taus[i][0]) + ", gG=" + str(good_taus[i][1]), color='black', linewidth=0.5)
    ax.plot([i for i in range(len(good_taus[0][2]))], good_taus[len(good_taus)-1][2],
            label="gA=" + str(good_taus[len(good_taus)-1][0]) + ", gG=" + str(good_taus[len(good_taus)-1][1]), color='red', linewidth=3)
    ax.legend(prop={'size': 6})
    ax.set_title(r'$\tau$' + " for gA and gG candidate values (DS to AS)", size=15)
    ax.set_xlabel(r"$t_i$")
    ax.set_ylabel(r"$\tau$")
    plt.show()
"""
########################################################################################################
MAIN PROGRAM
########################################################################################################
"""
# ~~~~~~ Uncomment to test functionality ~~~~~~~~

#M, S, _, _ = MSI(10 * nsiemens, 18 * nsiemens)
#tau_calc2(M,S)
#SRI() #plot the spiketrains
#plot_diff() #plot differences in spike trains
#plot_S_I_diffs() #plot prev + current spike train differences
#plot_S_I_diffs2() #plot prev + current spike train differences
#plot_tau_diff()
#a = plot_tau_diff2()
#import pickle
#with open('candidate_params.txt', 'wb') as f: pickle.dump(a, f)
#plot_synapse_conductances() #plot synaptic conductances vs time
#heatmap(10,10,2) #plot heat map of mean synaptic conductances
#print(check_for_candidate([-2,-2,4,3,2,4,5,-2,-3,-5,-3,-1,3,2,34,5,2,3,4]))

plot_periodicity()