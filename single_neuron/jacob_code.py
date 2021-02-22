#Update Feb 20
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


"""
#Start with HH, adapt to our paper

from brian2 import *

#HH PARAMS
num_neurons = 2
duration = 2*second

# Parameters
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2 * area
gl = 5e-5*siemens*cm**-2 * area
El = -65*mV
EK = -90*mV
ENa = 50*mV
g_na = 100*msiemens*cm**-2 * area
g_kd = 30*msiemens*cm**-2 * area
VT = -63*mV

# The model
eqs = Equations('''
dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + I)/Cm : volt
dm/dt = 0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m : 1
dn/dt = 0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
I : amp
''')


#INITIALIZE 3 NEURONS
# Threshold and refractoriness are only used for spike counting
M = NeuronGroup(num_neurons, eqs,
                    threshold='v > -40*mV',
                    refractory='v > -40*mV',
                    method='exponential_euler')
M.v = El
M.I = '0.7*nA * i / num_neurons'

I = NeuronGroup(num_neurons, eqs,
                    threshold='v > -40*mV',
                    refractory='v > -40*mV',
                    method='exponential_euler')
I.v = El
I.I = '0.7*nA * i / num_neurons'

S = NeuronGroup(num_neurons, eqs,
                    threshold='v > -40*mV',
                    refractory='v > -40*mV',
                    method='exponential_euler')
S.v = El
S.I = '0.7*nA * i / num_neurons'




#PAPER PARAMS
# Reversal potentials
El = 10.6*mV #leak reversal potential
ENa = 115*mV #sodium reversal potential
EK = -12*mV #potassium reversal potential
EA = 60*mV
EG = -20*mV
# Conductances
gl = 2.7*np.pi*msiemens #leak current conductance
gNa = 1080*np.pi*msiemens #sodium current max conductance
gK = 324*np.pi*msiemens #potassium current max conductance
gA = 10 * msiemens
gG = 20 * msiemens # variable
# Constant current, determines neuron excitability
Iext_S = 28e-11 *amp
Iext_I = 28e-11 *amp
Iext_R = 28e-11 *amp # variable
Iext = [Iext_S, Iext_I, Iext_R]
# Membrane capacitance
Cm = pi*9e-6 *farad
alphaA = 1.1 * mM**-1*ms**-1 #FIGURE OUT THESE UNITS
betaA = 0.19 * ms**-1 #UNITS
alphaG = 5.0 * mM**-1*ms**-1 #UNITS
betaG = 0.30 * ms**-1 #UNITS
gNa0 = 120*msiemens/cm**2


#INITIALIZE SYNAPSES
# Excitatory synapse
A1 = Synapses(M,S,
             model = ''' 
             dr/dt = (alphaA*T*(1-r)-betaA*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',

              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')

A2 = Synapses(S,I,
             model = ''' 
             dr/dt = (alphaA*T*(1-r)-betaA*r) : 1
              T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              on_post = '''
              I = gNa0*r*(v_post - EA)*(meter**2) ''')

# Excitatory synapses (M to S, S to I)
A1.connect()
A2.connect()

# Inhibitory synapse
B1 = Synapses(I,S,
             model =
             '''dr/dt = (alphaG*T*(1-r)-betaG*r) : 1
             T = 1/(1+exprel(-((0.001*v_pre/mV)-62)/5))*mM : mM''',
              on_post = '''
             I = gK*r*(v_post - EG) ''')

# Inhibitory synapses from Interneuron to Reveiver (2 to 1)
B1.connect()

monitor = SpikeMonitor(M)
SM = StateMonitor(M, 'v', record=1)
SI = StateMonitor(I, 'v', record=1)
SS = StateMonitor(S, 'v', record=1)
run(duration)

#plot(group.I/nA, monitor.count / duration)
#xlabel('I (nA)')
#ylabel('Firing rate (sp/s)')
#plt.show()

plot(SM.v[0], color = 'red')
plot(SI.v[0], color = 'blue')
plot(SS.v[0], color = 'black')
plt.show()
"""
