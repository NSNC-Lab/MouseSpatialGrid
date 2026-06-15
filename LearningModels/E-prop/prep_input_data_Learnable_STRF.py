import yaml
import numpy as np
import scipy.io
from tqdm import tqdm
from argparse import ArgumentParser
from scipy.io import savemat
import matplotlib.pyplot as plt
import os


_PLOTTED_ONSET_POISSON_DIAGNOSTIC = False
_PLOTTED_OFFSET_POISSON_RASTER = False


def _select_first_trace(values):
    values = np.asarray(values)
    if values.ndim == 1:
        return values
    return values[(slice(None),) + (0,) * (values.ndim - 1)]


def _select_first_raster(values):
    values = np.asarray(values)
    if values.ndim == 1:
        return values[:, None]
    if values.ndim >= 5:
        return values[:, 0, :, 0, 0]
    if values.ndim >= 3:
        return values[:, 0, :]
    return values.reshape(values.shape[0], -1)


def _select_low_rate_example(rate, spikes):
    rate = np.asarray(rate)
    spikes = np.asarray(spikes)

    if rate.ndim >= 3 and spikes.ndim >= 5:
        rate_means = np.nanmean(rate, axis=0)
        spike_counts = np.sum(spikes[:, 0, :, :, :], axis=(0, 1))
        if rate_means.shape == spike_counts.shape:
            valid = np.isfinite(rate_means)
            visible = valid & (spike_counts > 0)
            candidates = visible if np.any(visible) else valid
            scores = np.where(candidates, rate_means, np.inf)
            flat_idx = int(np.argmin(scores))
            cell_idx, batch_idx = np.unravel_index(flat_idx, rate_means.shape)
            all_cell_trial_raster = np.transpose(spikes[:, 0, :, :, batch_idx], (2, 1, 0))
            all_cell_trial_raster = all_cell_trial_raster.reshape(-1, spikes.shape[0])
            return {
                'rate_trace': rate[:, cell_idx, batch_idx],
                'spike_raster': spikes[:, 0, :, cell_idx, batch_idx],
                'compact_raster': all_cell_trial_raster,
                'mean_rate': float(rate_means[cell_idx, batch_idx]),
                'spike_count': int(spike_counts[cell_idx, batch_idx]),
                'cell_idx': cell_idx,
                'batch_idx': batch_idx,
            }

    rate_trace = _select_first_trace(rate)
    spike_raster = _select_first_raster(spikes)
    return {
        'rate_trace': rate_trace,
        'spike_raster': spike_raster,
        'compact_raster': spike_raster.T,
        'mean_rate': float(np.nanmean(rate_trace)),
        'spike_count': int(np.sum(spike_raster)),
        'cell_idx': 0,
        'batch_idx': 0,
    }


def _diagnostic_output_dir():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'strf_diagnostic_plots')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _save_and_show(fig, filename):
    save_path = os.path.join(_diagnostic_output_dir(), filename)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show(block=False)
    plt.pause(0.001)
    return save_path


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


    def gen_poisson_inputs(self, spks, response_name=None):
        """
        Generate Poisson spike inputs from spike rates.

        Parameters:
            spks : np.ndarray
                Spike rate matrix (time x neurons x trials)
            response_name : str
                Optional response label for one-time diagnostic plotting.

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

        self.plot_poisson_diagnostic(rate, s, response_name)
            
        return s,rate


    def plot_poisson_diagnostic(self, rate, spikes, response_name):
        global _PLOTTED_ONSET_POISSON_DIAGNOSTIC, _PLOTTED_OFFSET_POISSON_RASTER
        if response_name == 'onset':
            if _PLOTTED_ONSET_POISSON_DIAGNOSTIC:
                return
            _PLOTTED_ONSET_POISSON_DIAGNOSTIC = True
            raster_filename = 'onset_poisson_raster.png'
        elif response_name == 'offset':
            if _PLOTTED_OFFSET_POISSON_RASTER:
                return
            _PLOTTED_OFFSET_POISSON_RASTER = True
            raster_filename = 'offset_poisson_raster.png'
        else:
            return

        example = _select_low_rate_example(rate, spikes)
        rate_trace = example['rate_trace']
        time_ms = np.arange(rate_trace.shape[0]) * self.dt

        if response_name == 'onset':
            fig_rate, ax_rate = plt.subplots(figsize=(10, 4))
            ax_rate.plot(time_ms, rate_trace, linewidth=1.5)
            ax_rate.set_xlabel('Time (ms)')
            ax_rate.set_ylabel('Firing rate (Hz)')
            ax_rate.set_title(
                f"Low-FR onset profile: cell {example['cell_idx']}, batch {example['batch_idx']} "
                f"(mean {example['mean_rate']:.3f} Hz)"
            )
            ax_rate.grid(True, alpha=0.25)
            _save_and_show(fig_rate, 'onset_firing_rate_profile.png')

        self.plot_compact_poisson_raster(example, time_ms, raster_filename)


    def plot_compact_poisson_raster(self, example, time_ms, filename):
        compact_raster = np.asarray(example['compact_raster']) > 0
        row_idx, time_idx = np.nonzero(compact_raster)
        spike_times_ms = time_idx * self.dt
        row_positions = row_idx + 1
        fig_raster, ax_raster = plt.subplots(figsize=(6.0, 1.1))
        if spike_times_ms.size:
            x_points = np.vstack((spike_times_ms, spike_times_ms, np.full_like(spike_times_ms, np.nan))).ravel(order='F')
            y_points = np.vstack((
                row_positions - 0.38,
                row_positions + 0.38,
                np.full_like(row_positions, np.nan, dtype=float),
            )).ravel(order='F')
            ax_raster.plot(
                x_points,
                y_points,
                color='#4ea3d8',
                linewidth=0.28,
                solid_capstyle='butt',
            )
        ax_raster.set_xlim(time_ms[0], time_ms[-1])
        ax_raster.set_ylim(0.5, compact_raster.shape[0] + 0.5)
        ax_raster.invert_yaxis()
        ax_raster.set_xlabel('Time (ms)')
        ax_raster.set_ylabel('')
        ax_raster.set_yticks([])
        ax_raster.grid(False)
        ax_raster.tick_params(axis='both', labelsize=7, length=2)
        for spine in ax_raster.spines.values():
            spine.set_linewidth(0.8)
        _save_and_show(fig_raster, filename)
    
    
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

                    on_poisson_spks, rate_on = self.gen_poisson_inputs(fr_target_on, response_name='onset')
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
                    off_poisson_spks, rate_off = self.gen_poisson_inputs(fr_target_off, response_name='offset')
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_poisson_spks'] = off_poisson_spks
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_rate'] = rate_off
                    spks_dict[f'locs_masker_{locs[0]}_target_{locs[1]}_off'][f'stimulus_{stimulus}_rate_deriv'] = fr_target_off_deriv
           
            # Generate spiking activity
            spk_noise = self.gen_poisson_times(self.chans, self.fr, self.std, self.simlen, self.trials, self.num_cells)
            spks_dict[f'noise_masker_{locs[0]}_target_{locs[1]}'] = spk_noise
                
        return spks_dict



    
    
    
