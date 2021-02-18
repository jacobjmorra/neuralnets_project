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