import yaml
import numpy as np
import scipy.io
from tqdm import tqdm
from argparse import ArgumentParser
from scipy.io import savemat


class PrepInput(object):
    '''
    Class to prepare input data for IC neuron simulations and generating Poisson spikes based on STRF data.
    Parameters:
        args : argparse.Namespace
            Command-line arguments containing parameters like number of channels, trials, etc.
        config : dict
            Configuration dictionary loaded from a YAML file containing parameters for spike generation.
    '''
    def __init__(self, args, config,batch_size):
        
        self.chans = args.chans
        self.sigma = args.sigma
        self.trials = args.trials
        self.padToTime = args.padToTime
        self.dt = args.dt
        self.strfGain = args.strfGain
        self.std = args.std
        self.fr = args.FR
        self.num_cells = args.num_cells
        
        #self.fr = p[4,:]
        self.batch_size = batch_size
        self.simlen = args.simlen
        
        self.dt = config['dt']
        self.t_ref = config['t_ref']
        self.t_ref_rel = config['t_ref_rel']
        self.rec = config['rec']
        self.refrac = config['refrac']
        self.scale_factor = config['scale_factor']
        self.offset_val = config['offset_val']
        
        self.azi = np.flip(np.linspace(-90,90, self.chans))
        
        self.spatialCurves = self.genSpatiallyTunedChans()

    def genSpatiallyTunedChans(self):
        '''
        Generate spatial tuning curves for each channel based on Gaussian profiles.
        '''
        spatialCurves = np.zeros((self.chans, self.chans))
        for idx in range(self.chans): 
            spatialCurves[idx,:]= np.exp(((-1/2)*(self.azi-self.azi[idx])**2)*(1/self.sigma)) 
            
        return spatialCurves

    def make_grid_target_masker_locs(self):
        '''
        Create grid (coordinates) of masker and target locations based on number of channels.
        '''
        grid_x, grid_y = np.meshgrid(np.arange(len(self.azi)), np.arange(len(self.azi))) # masker, target
        masker_locs = grid_x.flatten() # masker
        target_locs = grid_y.flatten() # target
        return masker_locs, target_locs  

    def gen_IC_spks(self, tmax, locs, fr_targets, newStrfGain, strfGain):
        '''
        Generate adjusted IC based on STRF data for given masker and target locations.
        Parameters:
            tmax : int
                Maximum time for spike train (number of time steps).
            locs : tuple
                (masker_location, target_location)
            fr_targets : np.ndarray
                Firing rates for target stimulus.
            fr_masker : np.ndarray
                Firing rates for masker stimulus.
            newStrfGain : float
                New STRF gain to adjust firing rates.
            strfGain : float    
                Original STRF gain.
        Returns:    
            spks : np.ndarray 
                Generated spike train matrix (time x channels x trials).
        '''

        m_loc, t_loc = locs # masker, target
        singleConfigSpks = np.zeros((self.spatialCurves.shape[0], tmax))  
        

        for ch in range(self.spatialCurves.shape[0]):
            t_wt = self.spatialCurves[ch,self.azi==self.azi[t_loc]] if t_loc is not None else 0.
            m_wt = self.spatialCurves[ch,self.azi==self.azi[m_loc]] if m_loc is not None else 0.

            #if t_wt + m_wt == 0: raise Exception(f'No contribution from target or masker at channel {ch}, trial {t+1}. Both cannot be None.')
                
            singleConfigSpks[ch,:] = t_wt*fr_targets.squeeze()   

            if t_wt + m_wt >= 1: singleConfigSpks[ch,:] = singleConfigSpks[ch,:] / (t_wt + m_wt)
        
        if singleConfigSpks.shape[1] < self.padToTime/self.dt:
            padSize = int(self.padToTime/self.dt)-singleConfigSpks.shape[1]
            singleConfigSpks = np.concatenate([singleConfigSpks, np.zeros((self.spatialCurves.shape[0],padSize))], axis=1)
        
        spks =  singleConfigSpks * newStrfGain / strfGain  # Shape to (time, chans, trials)

        return spks



    

    def spike_generator(self, rate):
        """
        Generate a Poisson spike train with an absolute and relative refractory period.
        parameters:
            rate : np.ndarray
        """
        dt_sec = self.dt / 1000  # ms to seconds

        n = len(rate)
        spike_train = np.zeros(n)
        spike_times = []

        #9/17 refractory preiod seems to low compared to real data. Perhaps extend this?
        n_refab = int(0 / 1000 / dt_sec)  # number of samples for ref. period window
        #n_refab = int(15 / 1000 / dt_sec)  # number of samples for ref. period window
        #n_refab = int(30 / 1000 / dt_sec)  # number of samples for ref. period window
        tw = np.arange(n_refab + 1)

        t_ref_samp = int(self.t_ref / 1000 / dt_sec)
        t_rel_samp = int(self.t_ref_rel / 1000 / dt_sec)


        # Recovery function based on Schaette et al. 2005
        with np.errstate(divide='ignore', invalid='ignore'):
            w = np.power(tw - t_ref_samp, self.rec) / (
                np.power(tw - t_ref_samp, self.rec) + np.power(t_rel_samp, self.rec)
            )
        w[tw < t_ref_samp] = 0
        w = np.nan_to_num(w)

        x = np.random.rand(n)
        

        for i in range(n):
            if spike_times and i - spike_times[-1] < n_refab:
                rate[i] *= w[i - spike_times[-1]]
            if x[i] < dt_sec * rate[i]:
                spike_train[i] = 1
                spike_times.append(i)

        return spike_train
    
    
    def gen_poisson_times(self, chans, FR, std, simlen, trials, num_cells):
        """
        Generate Poisson spike trains with refractory period.

        Parameters:
            chans : int
                Number of neurons in population.
            FR : float
                Firing rate in Hz.
            std : float
                Standard deviation of the firing rate.
            simlen : int
                Number of time steps (default = 35000)

        Returns:
            token : np.array
                Binary spike train matrix of shape (simlen, chans)
        """
        chans = int(chans)
        simlen = int(simlen)
        std = int(std) #set to 0.0 right now
        #FR = int(FR) # set to  8.0 right now

        # Generate Poisson spikes with added noise

        #Convert things from tensors so we can work with them

        #print(np.shape(np.reshape(FR,(100,1))[None][None]))
        #print(np.shape(std * np.random.randn(self.batch_size,simlen, chans, trials)))

        #firing_rate_reshape = np.reshape(FR,(100,1))

        #rand_gauss = firing_rate_reshape[:,:,None,None] + std * np.random.randn(self.batch_size,simlen, chans, trials)
        
        #rand_gauss = FR + std * np.random.randn(simlen, chans, trials)
        #rand_bin = np.random.rand(self.batch_size,simlen, chans, trials) < (rand_gauss * self.dt / 1000)
        firing_rate_reshape = np.reshape(FR,(num_cells,1,1,1,1))
        rand_gauss = firing_rate_reshape + std * np.random.randn(num_cells, self.batch_size,simlen, chans, trials)
        rand_bin = np.random.rand(num_cells,self.batch_size,simlen, chans, trials) < (rand_gauss * self.dt / 1000)

        temp = rand_bin.astype(np.uint8)

        for num in range(num_cells):
            for batch in range(self.batch_size):
                for chan in range(chans):
                    for trial in range(trials):
                        spk_inds = np.where(temp[num,batch,:, chan, trial])[0]
                        if len(spk_inds) > 1:
                            ISIs = np.diff(spk_inds) * self.dt
                            violate_inds = np.where(ISIs < self.refrac)[0] + 1
                            temp[num,batch,spk_inds[violate_inds], chan, trial] = 0

        
        #print(np.shape(np.array(temp)))

        if num_cells == 1:
            return np.squeeze(temp, axis=0)

        return np.array(temp)


    def gen_poisson_inputs(self, spks):
        """
        Generate Poisson spike inputs from spike rates.

        Parameters:
            spks : np.ndarray
                Spike rate matrix (time x neurons x trials)

        Returns:
            s : np.ndarray
                Binary spike train matrix (time x neurons)
        """
        
        trial_rate = spks  # shape: (time, neurons, trials)

        rate = trial_rate[int(self.offset_val):int(self.offset_val+len(trial_rate)*self.scale_factor)] if self.scale_factor != 1 else trial_rate

        #s = np.zeros_like(rate)
        #for chan in range(rate.shape[1]):
        #    for trial_num in range(rate.shape[2]):
        s = np.zeros((rate.shape[0],self.chans,self.trials,self.num_cells,self.batch_size))
        for chan in range(self.chans):
            for trial_num in range(self.trials):
                for batch in range(self.batch_size):
                    for cell in range(self.num_cells):
                        s[:, chan, trial_num, cell, batch] = self.spike_generator(rate[:, cell, batch])
            
        return s,rate
    
    
    def process_input(self, strf_path, list_locs, on_neuron=True, off_neuron=True, target_dict=None):
        '''
        Process input STRF data to generate spike trains for specified masker and target locations.
        Parameters:
            strf_path : str
                Path to the .mat file containing STRF data.
            list_locs : list of tuples
                List of (masker_location, target_location) pairs.
            on_neuron : bool
                Whether to generate spikes for ON neurons.
            off_neuron : bool
                Whether to generate spikes for OFF neurons.
            
            Returns:    
                spks_dict : dict
                    Dictionary containing generated spike trains (adjusted IC and poisson spikes ) for each location and stimulus type.
                
                    '''
        #data = scipy.io.loadmat('C:/Users/ipboy/Documents/GitHub/MouseSpatialGrid/LearningModels/E-prop/strfMat.mat')
        data = target_dict

        fr_target_on = data['fr_on0']
        fr_target_off = data['fr_off0']
        fr_target_on_deriv = data['fr_on_deriv0']
        fr_target_off_deriv = data['fr_off_deriv0']

        #fr_target_on = np.array([data['fr_on0'][0][0],data['fr_on1'][0][0]])
        #fr_target_off = np.array([data['targets']['fr_off0'][0][0],data['targets']['fr_off1'][0][0]])


        #It looks like previously we just used to offset fr for our maskers??
        #fr_masker = []
        #for k in range(10):
            #print(data['maskers'][f'fr_off{k}'][0][0])
        #    fr_masker.append(data['maskers'][f'fr{k}'][0][0])



        #strfGain = self.strfGain
        tmax = fr_target_on.shape[1]
        #newStrfGain = strfGain


        spks_dict = {}
        for locs in list_locs: 


            if on_neuron:
                spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_on'] = {}
                #for stimulus in range(fr_target_on.shape[0]):
                for stimulus in range(1):
                    #on_spks = self.gen_IC_spks(
                    #    tmax=tmax,
                    #    locs=locs,
                    #    fr_targets=fr_target_on[stimulus],
                    #    newStrfGain=newStrfGain,
                    #    strfGain=strfGain
                    #)

                    on_poisson_spks, rate_on = self.gen_poisson_inputs(fr_target_on)
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_on'][f'stimulus_{stimulus}_poisson_spks'] = on_poisson_spks
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_on'][f'stimulus_{stimulus}_rate'] = rate_on
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_on'][f'stimulus_{stimulus}_rate_deriv'] = fr_target_on_deriv
                
            if off_neuron:
                spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'] = {}
                #for stimulus in range(fr_target_off.shape[0]):
                for stimulus in range(1):
                    #off_spks = self.gen_IC_spks(
                    #    tmax=tmax,
                    #    locs=locs,
                    #    fr_targets=fr_target_off[stimulus],
                    #    newStrfGain=newStrfGain,
                    #    strfGain=strfGain
                    #)
                    off_poisson_spks, rate_off = self.gen_poisson_inputs(fr_target_off)
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_poisson_spks'] = off_poisson_spks
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_rate'] = rate_off
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_rate_deriv'] = fr_target_off_deriv
           
            # Generate spiking activity
            spk_noise = self.gen_poisson_times(self.chans, self.fr, self.std, self.simlen, self.trials, self.num_cells)
            spks_dict[f'noise_masker_{locs[0]}_target_{locs[1]}'] = spk_noise
                
        return spks_dict



    
    
    