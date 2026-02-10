import numpy as np
import matplotlib.pyplot as plt

#Input current 0.17 and ginc 0.00057 gives around the ISI changes reported in 

#Try lower voltage reset potentially or lower E_L. 





#ginc_trials = [0,0.06]

Voltage_holder = []
Hz_Holder = []
Hz_Holder2 = []


# for trial2 in range(2):


    # ginc = ginc_trials[trial2]
    # isis_ratio_holder = []

    # for trial in range(20):

thresh = -55
dt = 0.1 #Integration time step

R = 100 #Resistance

tau = 6 #Membrane time constant

smallinjected = 0.000025
largeinjected = 0.000150
current_30hz  = 0.003300  
current_50hz  = 0.022000         
theoretical_2k_current_drive = 5 #nA


#print(input_current)
E_L = -65 #Leak
Voltage = E_L #Resting voltage
V_reset = E_L #Reset voltage after spike
input_current = (thresh-E_L)/R + 0.000025 #0.2 #?
E_k = -80 #Adaptation reversal potential
tau_ad = 100 #Adaptation time constant
grsa = 0
#ginc = 0
ginc = 0.0002/R #Adaptation increment (taken from Dayan and abbott ~0.06 rm*ginc)

#input_current = 0.19 + trial*0.001

spike_times = []
for k in range(10000):
    #No adaptation
    #Voltage_k = ((E_L-Voltage) + R*input_current)/tau

    #With adaptation
    Voltage_k = ((E_L-Voltage) + R*input_current - R*grsa*(Voltage-E_k))/tau

    Voltage = Voltage + dt*Voltage_k
    grsa = grsa - dt*(grsa/tau_ad)
    
    if Voltage > thresh:
        spike_times.append(k)
        Voltage = V_reset
        grsa = grsa + ginc
        Voltage_holder.append(0)
    else:
        Voltage_holder.append(Voltage)
        # if trial2 == 0:
        #     Hz_Holder.append(10000/np.diff(spike_times).mean())
        # else:
        #     isis = np.diff(spike_times)
        #     Hz_Holder2.append(10000/isis.mean())

        #     if len(isis) > 2:
        #         #print(isis)
        #         isis_ratio_holder.append(isis[0]/isis[-1])
        #         #print(isis[0]/isis[-1])
        #     else:
        #         isis_ratio_holder.append(0)

print(len(spike_times))
#Calculate and plot instantaneous firing rate
isi = np.diff(spike_times)
print(isi[0]/isi[-1])

#print(spike_times)

plt.plot(np.arange(0,10000*dt,dt),Voltage_holder)
plt.xlabel('time (ms)')
plt.ylabel('Voltage (mV)')
plt.title(f'Current injection for 50hz (adaptation ratio = {isi[0]/isi[-1]:.2f})')
plt.ylim(-80, 0)
plt.show()




#print(isi)

# plt.scatter(np.arange(20)*0.001 + 0.17, Hz_Holder)
# plt.scatter(np.arange(20)*0.001 + 0.17, Hz_Holder2)
# plt.xlabel('Input Current (nA)')
# plt.ylabel('Firing Rate (Hz)')
# plt.show()


# plt.scatter(np.arange(20)*0.001 + 0.17, isis_ratio_holder)
# plt.xlabel('Input Current (nA)')
# plt.ylabel('Adaption Ratio (1st ISI / Last ISI)')
# plt.show()

#print(spike_times)