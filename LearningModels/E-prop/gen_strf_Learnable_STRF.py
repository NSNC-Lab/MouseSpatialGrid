import os
import numpy as np
try:
        import cupy as cp
except ImportError:
        cp = None
import math
import sys
import yaml
import copy

from scipy.io import wavfile
from tqdm import tqdm, trange
from argparse import ArgumentParser
from scipy.signal import resample_poly

import matplotlib.pyplot as plt


_PLOTTED_STIMULUS_SPECTROGRAM = False
_PLOTTED_STRF_DIAGNOSTICS = False


def rms(data):
        """
        Calculates the Root Mean Square (RMS) of a NumPy array.

        Args:
        data (np.ndarray): The input_ array.

        Returns:
        float or np.ndarray: The RMS value(s).
        """
        return np.sqrt(np.mean(data**2))


def get_f(fftLen, samprate, tfParams):
        select = np.linspace(1,int((fftLen+1)/2), int((fftLen+1)/2)) if fftLen%2!=0 else np.linspace(1,int(fftLen/2)+1, int(fftLen/2)+1)
        f0 = (select-1)*samprate/fftLen
        
        maxIndx = np.where(f0 >= tfParams['high_freq'])[0][0]
        minIndx = np.where(f0 < tfParams['low_freq'])[0][-1]+1
        
        fstep = f0[1]
        f = np.linspace(f0[minIndx], f0[maxIndx], 1 + int((f0[maxIndx]-f0[minIndx]+1)/fstep)) #f0(minIndx):fstep:f0(maxIndx)
                
        return f


def _plot_line_trace(ax, x, y, xlabel, ylabel, title):
        ax.plot(x, y, linewidth=1.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.25)


def _select_first_trace(values):
        values = np.asarray(values)
        if values.ndim == 1:
                return values
        return values[(slice(None),) + (0,) * (values.ndim - 1)]


def _select_first_matrix(values):
        values = np.asarray(values)
        if values.ndim <= 2:
                return np.squeeze(values)
        return values[(slice(None), slice(None)) + (0,) * (values.ndim - 2)]


def _flat_index_of_min(values):
        values = np.asarray(values, dtype=np.float64).reshape(-1)
        if values.size == 0:
                return 0
        finite_values = np.where(np.isfinite(values), values, np.inf)
        return int(np.argmin(finite_values))


def _select_trace_at_flat_index(values, flat_index):
        values = np.asarray(values)
        if values.ndim == 1:
                return values
        columns = values.reshape(values.shape[0], -1)
        return columns[:, flat_index % columns.shape[1]]


def _select_matrix_at_flat_index(values, flat_index):
        values = np.asarray(values)
        if values.ndim <= 2:
                return np.squeeze(values)
        matrices = values.reshape(values.shape[0], values.shape[1], -1)
        return matrices[:, :, flat_index % matrices.shape[2]]


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


class GenSTRF(object):
        
        '''
        Class to generate STRFs from stimuli
        Parameters:
                args : argparse.Namespace
                        Command-line arguments containing more dynamic parameters.
                config : dict
                        Configuration dictionary loaded from a YAML file containing static parameters.
                        
                init_stim_path:
                        reference stimulus path
        '''
        
        def __init__(self, args, config, init_stim_path,learnable_params,batch_size,num_cells):
                
                self.args = args
                self.strf_config = config
                self.paramH = self.strf_config['paramH']
                self.paramG = self.strf_config['paramG']
                self.preprocStimParams = self.strf_config['preprocStimParams']
                self.strfGain = learnable_params[0,:]
                self.strfTimeConstant = learnable_params[1,:] #np.squeeze(np.ones((1, batch_size))*0.01)
                self.batch_size = batch_size
                self.num_cells = num_cells
                # generating strfs with a reference stimlus
                #self.strf = self.process_initial_stimulus(init_stim_path, self.strf_config['targetlvl'], self.strf_config['strfGain'])
                #With learnable parameters
                self.strf = self.process_initial_stimulus(init_stim_path, self.strf_config['targetlvl'], self.strfGain)
                self._plotted_stimulus_spectrogram = False

        def process_initial_stimulus(self, stim_path, lvl, strfGain):
                '''
                Computes STRFs based on a reference stimulus
                
                Parameters:
                        stim_path: str
                                path to the stimulus
                        lvl: float
                        strfGain: float

                '''
                fs, data = wavfile.read(stim_path)
                data = data[np.newaxis, :]

                #Append zeros to the beginning (1 second worth)
                #data = np.reshape(np.append(np.zeros((1,200000)), data),(1,846011))

                spec, t, f = self.STRFspectrogram(data/rms(data)*lvl,fs)                 
                strf = self.STRFgen(self.paramH, self.paramG, f, t[1]-t[0]) 
                #Learnable version
                #strf['w1'] = strf['w1'][:,:,None]*np.squeeze(strfGain)[None,None,:]  #The idea here is that for the STRF we have some spectrogram that is 2500 by ~50 and we want to scale it by the strfGain across each batch trial. So we create fake axis and broadcast across them.
                #strf['w2'] = strf['w2'][:,:,None]*np.squeeze(strfGain)[None,None,:]  
                gain = np.asarray(strfGain)
                if self.num_cells > 1:
                        gain = gain.reshape(self.num_cells, self.batch_size)[None, None, :, :]
                else:
                        gain = gain.reshape(self.batch_size)[None, None, :]

                strf['w1'] = strf['w1'] * gain
                strf['w2'] = strf['w2'] * gain

                self.plot_strf_diagnostics(f, strf, stim_path)
                
                return strf
        
        def process_stimulus(self, stim_path, lvl, stimGain):
                '''
                Computes IC values based for onset and offset neurons
                
                Parameters:
                        stim_path: str
                        path to the stimulus
                        lvl: float
                        strfGain: float

                '''
                
                fs, data = wavfile.read(stim_path)


                #plt.plot(data)
                #plt.show()

                #print(fs)
                data = data.astype(np.float64)

                #Add a second to the front of the signal ******FOR LOOKING AT PRE-sim response*******
                #data = np.append(np.zeros((1,200000)), data)
                #print(np.shape(data))

                # #Rsample the data to match the simulation timestep
                
                #For whatever reason in matlab it samples stuff in at 195khz and here it does 200k, so no need to upsample. --- might be because we are pulling from already upsampled stimuli in this version.
                #p, q = 12500, 12207
                #data = resample_poly(data, p, q)

                #print(np.shape(data))

                #stim_pad = 2500 * 20  #This 2500 represents the length of the STRF kernel. \TODO implement changing this

                #pad = np.zeros(stim_pad, dtype=data.dtype)
                #data = np.concatenate([pad, data])

                #print(np.shape(data))

                data = data[np.newaxis, :]

                #print(np.shape(data))

                spec, t, f = self.STRFspectrogram(data/rms(data)*lvl,fs)
                self.plot_stimulus_spectrogram(spec, t, f, stim_path)

                fr_on, fr_off, fr_on_deriv, fr_off_deriv = self.STRFconvolve(self.strf,spec*stimGain, self.strf_config['mean_rate'])

                #print(np.shape(fr_on))
                
                #print(np.shape(fr_on))

                fr_on = fr_on[2500:]
                fr_off = fr_off[2500:]
                fr_on_deriv = fr_on_deriv[2500:]
                fr_off_deriv = fr_off_deriv[2500:]

                #print(np.shape(fr_on))

                return fr_on, fr_off, fr_on_deriv, fr_off_deriv


        def plot_stimulus_spectrogram(self, spec, t, f, stim_path):
                global _PLOTTED_STIMULUS_SPECTROGRAM
                if self._plotted_stimulus_spectrogram or _PLOTTED_STIMULUS_SPECTROGRAM:
                        return

                self._plotted_stimulus_spectrogram = True
                _PLOTTED_STIMULUS_SPECTROGRAM = True
                fig, ax = plt.subplots(figsize=(10, 4))
                extent = [t[0], t[-1], f[0], f[-1]]
                im = ax.imshow(spec.T, aspect='auto', origin='lower', extent=extent)
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('Frequency (Hz)')
                ax.set_title(f'Input stimulus spectrogram: {os.path.basename(stim_path)}')
                fig.colorbar(im, ax=ax, label='Spectrogram power')
                _save_and_show(fig, 'input_stimulus_spectrogram.png')


        def plot_strf_diagnostics(self, f, strf, stim_path):
                global _PLOTTED_STRF_DIAGNOSTICS
                if _PLOTTED_STRF_DIAGNOSTICS:
                        return

                _PLOTTED_STRF_DIAGNOSTICS = True

                time_constants = np.asarray(strf.get('time_constant', []), dtype=np.float64).reshape(-1)
                temporal_example_idx = _flat_index_of_min(time_constants)
                selected_time_constant = time_constants[temporal_example_idx] if time_constants.size else np.nan

                temporal_kernel = _select_trace_at_flat_index(strf['H'], temporal_example_idx)
                temporal_t = np.asarray(strf['t'])[:temporal_kernel.shape[0]]
                frequency_kernel = _select_first_trace(strf['G'])
                w1_matrix = _select_matrix_at_flat_index(strf['w1'], temporal_example_idx)

                fig_temporal, ax_temporal = plt.subplots(figsize=(8, 4))
                _plot_line_trace(
                        ax_temporal,
                        temporal_t,
                        temporal_kernel,
                        'Delay (s)',
                        'Weight',
                        f'STRF temporal kernel (short tau = {selected_time_constant * 1000:.2f} ms)',
                )
                _save_and_show(fig_temporal, 'strf_temporal_kernel.png')

                fig_frequency, ax_frequency = plt.subplots(figsize=(8, 4))
                _plot_line_trace(
                        ax_frequency,
                        f[:frequency_kernel.shape[0]],
                        frequency_kernel,
                        'Frequency (Hz)',
                        'Weight',
                        'STRF frequency kernel',
                )
                _save_and_show(fig_frequency, 'strf_frequency_kernel.png')

                fig_w1, ax_w1 = plt.subplots(figsize=(8, 5))
                w1_extent = [temporal_t[0], temporal_t[-1], f[0], f[-1]]
                im_w1 = ax_w1.imshow(w1_matrix, aspect='auto', origin='lower', extent=w1_extent)
                ax_w1.set_xlabel('Delay (s)')
                ax_w1.set_ylabel('Frequency (Hz)')
                ax_w1.set_title(f'w1 weight matrix / STRF (short tau = {selected_time_constant * 1000:.2f} ms)')
                fig_w1.colorbar(im_w1, ax=ax_w1, label='Weight')
                _save_and_show(fig_w1, 'strf_w1_weight_matrix.png')


        def GaussianSpectrum(self, input_, increment, winLength, samprate, nstd = 6):

                # Enforce even winLength to have a symmetric window
                if winLength%2 == 1: winLength = winLength +1

                # Make input_ it into a row vector if it isn't
                if input_.shape[0] > 1: input_ = input_.transpose()

                # Padd the input_ with zeros
                pinput_ = np.zeros([1,input_.shape[1]+winLength])
                pinput_[:, winLength//2:winLength//2+input_.shape[1]] = input_
                input_Length = pinput_.shape[1]

                # The number of time points in the spectrogram
                frameCount = math.floor((input_Length-winLength)/increment)+1

                # The window of the fft
                fftLen = winLength

                # Guassian window               
                wx2 = (np.linspace(1, winLength, winLength)-((winLength+1)/2))**2
                wvar = (winLength/nstd)**2
                ws = np.exp(-0.5*(wx2/wvar))

                # Initialize output "s" 
                if fftLen%2!=0:  s = np.zeros([int((fftLen+1)/2)+1, frameCount]) # winLength is odd
                else: s = np.zeros([int(fftLen/2)+1, frameCount], dtype=complex) # winLength is even 

                pg = np.zeros([1, frameCount])
                
                for i in trange(frameCount, desc='Computing fft - ', ncols=100):
                        start = i*increment
                        last = start + winLength
                        
                        f = ws*pinput_[:, start:last]
                        pg[:,i] = np.std(f)
                        
                        specslice = np.fft.fft(f)
                        s[:,i] = specslice[:,:(int((fftLen+1)/2)+1)] if fftLen%2!=0 else specslice[:,:(int(fftLen/2)+1)]
                
                # Assign frequency_label
                select = np.linspace(1,int((fftLen+1)/2), int((fftLen+1)/2)) if fftLen%2!=0 else np.linspace(1,int(fftLen/2)+1, int(fftLen/2)+1)
                f0 = (select-1)*samprate/fftLen

                # assign time_label
                t0 = np.linspace(0, s.shape[1]-1, s.shape[1])*(increment/samprate)
                return s, t0, f0, pg


        def timefreq(self, audioWaveform, sampleRate, typeName, params, samplingRatePoints):
                '''
                General purpose time-frequency representation function
                
                parameters:
                        audioWaveform: digitized sound pressure waveform
                
                        sampleRate: rate in Hz of digitized waveform
                
                        typeName: 'stft' for short-time fourier transforms
                        
                        params: depends on typeName, default values used if not given
        
                                'stft' params
                                fband: spacing between frequency bands of spectrogram (125)
                                nstd: # of std deviations that define width of Gaussian window (6)
                                low_freq: lowest frequency in Hz (250)
                                high_freq: highest frequency in Hz (8000)
                
                returns:        
                        tfrep: time-frequency representation
                        t: vector of time points for tf-representation
                        f: vector of frequencies for tf-representation
                        spec: matrix of values for tf-representation at a given time-frequency cell
                        params: the parameters that created the time-frequency representation
                '''
                
                tfrep = {}
                tfrep['params'] = params
                tfrep['params']['rawSampleRate'] = sampleRate

                if typeName == 'stft':
                        # compute raw complex spectrogram
                        twindow = tfrep['params']['nstd']/(tfrep['params']['fband']*2.0*math.pi)   # Window length
                        winLength = int(twindow*sampleRate)  # Window length in number of points
                        winLength = int(winLength/2)*2 # Enforce even window length
                        
                        increment = int(sampleRate/samplingRatePoints) # Sampling rate of spectrogram in number of points
                        
                        s, t0, f0, pg = self.GaussianSpectrum(audioWaveform, increment, winLength, sampleRate, tfrep['params']['nstd']) 
                        
                        # normalize the spectrogram within the specified frequency range
                        maxIndx = np.where(f0 >= tfrep['params']['high_freq'])[0][0]
                        minIndx = np.where(f0 < tfrep['params']['low_freq'])[0][-1]+1
                        
                        normedS = np.abs(s[minIndx:maxIndx+1, :]) # s is output spectrogram
                        
                        # set tfrep values
                        fstep = f0[1]
                        tfrep['t'] = t0
                        tfrep['f'] = np.linspace(f0[minIndx], f0[maxIndx], 1 + int((f0[maxIndx]-f0[minIndx]+1)/fstep)) #f0(minIndx):fstep:f0(maxIndx)
                        tfrep['spec'] = normedS # final output
                else:
                        raise Exception(f'Unknown time-frequency representation type: {typeName}')
                
                return tfrep


        def preprocSound(self, audioWaveform, params):
                '''
                Preprocess audio waveforms into strflab-ready time-frequency representation
                
                Parameters:
                        audioWaveforms: digitized sound waveform
                        
                        params:                        
                                rawSampleRate: sample rates of audio waveforms in Hz, not
                                required if audioWaveforms specifies .wav files.
                        
                                tfType: time-frequency type, either 'stft' for short-time
                                Fourier transform, others are not supported at the moment
                        
                                tfParams: time-frequency parameter structure, specific to the tfType
                                specified. For information on these params, see timefreq
                        
                        returns:
                                wholeStim: an NxM matrix representing the stimuli,
                                        where N is the # of time points and M is the # of features.
                                
                                groupIndex: an 1xN vector representing the group index of each
                                        time point, i.e. the stimulus/response pair it belongs to.
                                
                                stimInfo: information about the preprocessed stimulus:
                                        stimLengths = vector of lengths of each stimulus, in seconds
                                        sampleRate = sample rate of preprocessed stimulus in Hz
                                        numStimFeatures = number of stimulus features (M)
                                        tfType = time frequency representation, same as params['tfType']
                                        tfParams = same as params['tfParams'];
                                        f = vector of frequencies, one for each feature
                                
                                params: same params as passed into function
                '''
                # TODO: implement 'wavelet' and 'lyons' tf representations if needed
                allowedTypes = ['stft'] 
                if params['tfType'] not in allowedTypes:
                        raise Exception(f'Unknown time-frequency representation type: {params["tfType"]}')

                stimStructs = {}
                maxPower = -1
                
                
                
                stim = {}        
                stim['tfrep'] = self.timefreq(audioWaveform, params['rawSampleRate'], params['tfType'], params['tfParams'], params['samplingRatePoints'])  
                stim['stimLength'] = stim['tfrep']['spec'].shape[1] / params['stimSampleRate'] #size(stim.tfrep.spec, 2)
                stim['sampleRate'] = params['stimSampleRate']

                stimStructs[f'waveform']=stim
                
                if params['numStimFeatures'] == -1: params['numStimFeatures'] = stim['tfrep']['spec'].shape[0] #size(stim.tfrep.spec, 1)
                if params['f'] == -1: params['f'] = stim['tfrep']['f']
                
                params['totalStimLength'] = params['totalStimLength'] + stim['tfrep']['spec'].shape[1] #size(stim.tfrep.spec, 2)

                maxPower = np.max([maxPower, np.max(np.max(stim['tfrep']['spec'], axis=0))])        
                
                ## normalize spectrograms and take log if requested
                if params['tfParams']['log']:
                        refpow = maxPower if params['tfParams']['refpow'] == 0 else params['tfParams']['refpow']
                        for key in stimStructs.keys():          
                                stim = stimStructs[key]
                                stim['tfrep']['spec'] = 20*np.log10(stim['tfrep']['spec']/refpow)+params['tfParams']['dbnoise']
                                stim['tfrep']['spec'][stim['tfrep']['spec']<0]=0
                                stimStructs[key] = stim

                ## concatenate stims into big matrix and record info into struct
                stimInfo = {}
                stimInfo['stimLengths'] = np.zeros([1, audioWaveform.shape[1]]) #zeros(1, length(audioWaveforms))
                stimInfo['sampleRate'] = params['stimSampleRate']
                stimInfo['numStimFeatures'] = params['numStimFeatures']
                stimInfo['tfType'] = params['tfType']
                stimInfo['tfParams'] = params['tfParams']
                stimInfo['f'] = params['f']

                groupIndex = np.zeros([1, params['totalStimLength']])
                wholeStim = np.zeros([params['totalStimLength'], params['numStimFeatures']])
                
                cindx = 0
                for idx, key in enumerate(stimStructs.keys()):
                
                        stim = stimStructs[key]
                        slen = stim['tfrep']['spec'].shape[1] #size(stim.tfrep.spec, 2)
                        tend = cindx + slen
                        
                        wholeStim[cindx:tend, :] = stim['tfrep']['spec'].transpose()
                        groupIndex[:, cindx:tend] = idx
                        stimInfo['stimLengths'][idx] = slen / params['stimSampleRate']
                        
                        cindx = tend + 1

                return wholeStim, groupIndex, stimInfo, params

        def STRFspectrogram(self, stim, fs):
                '''
                Parameters
                        stim: array or cell of sound pressure waveform
                        fs:
                returns
                        stim_pec: STRF matrix
                        t: time array for STRF
                        f:  freq array for STRF
                '''
                
                params = copy.deepcopy(self.preprocStimParams)
                params['rawSampleRate']=fs
                params['totalStimLength'] = 0
                params['numStimFeatures'] = -1
                params['f']=-1
                
                ## use preprocSound to generate spectrogram
                stim_spec, groupIndex, stimInfo, preprocStimParams = self.preprocSound(stim, params)
                
                tInc = 1 / stimInfo['sampleRate'] # generate corresponding timeline for spectrogram
                
                t = np.linspace(0, (stim_spec.shape[0]-1)*tInc, 1 + int(((stim_spec.shape[0]-1)*tInc)/tInc))  # 0:tInc:(size(stim_spec, 1)-1)*tInc
                f=stimInfo['f']
                
                return stim_spec, t, f

        def  STRFgen(self, paramH, paramG, f, dt, maxdelay = 2500, nIn=1, outputNL='linear', freqDom=0):
                '''
                Generate STRFs with specified parameters
                Create STRF of [f] frequecy channels, and time delays of 40 dts
                Parameters from Amin et al., 2010, J Neurophysiol
                Temporal parameters from Adelson and Bergen 1985, J Opt Soc Am A
        
                Takes the number of inputs and outputs for a generalized linear model, and
                returns a strf dict. The weights are drawn from a zero mean,
                isotropic Gaussian, with variance scaled by the fan-in of the output
                units.
                
                parameters:
                        nin: number of inputs in one time slice of the strf
                        maxdelay: max latencies for the strf
                        outputNL: string indicating the output function and regression type:
                                                'linear'     (For use with squared error)
                        freqDom: toggles between options for convolution of stimulus with model weight.  For models with long delays
                                                (e.g., auditory data), try freqDom = 1.  For models with shorter delays,
                                                freqDom = 0 will usually produce faster results.
        
                returns:
                        strf: a STRF dict containing fields:
                                type: 'glm'
                                nin: number of inputs in one time slice of the strf (see above)
                                delays: delay vector input (see above)
                                nWts: total number of weights and biases
                                actfn: string describing the output unit activation function: 'linear' for now
                                freqDomain: 1 for convolution in frequency domain.
                                w1: nin x 1 vector of weights
                                b1: scalar bias term                
                '''
                strf = {}
                strf['type'] = 'lin'
                strf['nIn'] = nIn
                strf['t'] = np.linspace(0, maxdelay*dt, int((maxdelay*dt)/dt)) #0:dt:maxdelay*dt
                strf['delays'] = np.linspace(0, maxdelay, maxdelay) #0:maxdelay
                strf['nWts'] = (nIn*len(strf['delays']) + 1)

                # strf.w1=zeros(nIn,length(delays))
                strf['b1']=0

                nlSet=['linear', 'logistic', 'softmax', 'exponential']
                
                if outputNL in nlSet: strf['outputNL'] = outputNL
                else: raise Exception('linInit >< Unknown Output Nonlinearity!')

                strf['freqDomain'] = freqDom

                strf['internal'] = {}
                strf['internal']['compFwd']=1
                strf['internal']['prevResp'] = []
                strf['internal']['prevLinResp'] = []
                strf['internal']['dataHash'] = np.nan

                
                paramH['phase'] = np.pi/2

                ## Generate STRFs with specified parameters
                # Create STRF of [f] frequecy channels, and time delays of 40 dts
                # Parameters from Amin et al., 2010, J Neurophysiol
                # Temporal parameters from Adelson and Bergen 1985, J Opt Soc Am A   
                t = strf['t'][:,None] # time delay
                
                #just a reminder that we are replacing all paramH['alphas'] with the learnable time constant.
                #strf['H'] = np.exp(-t/paramH['alpha'])*(paramH['SC1']*(t/paramH['alpha'])**paramH['N1']/math.factorial(paramH['N1']) - \
                
                self.strfTimeConstant = self.strfTimeConstant[None,:]
                strf['time_constant'] = np.asarray(self.strfTimeConstant).copy()

                if self.num_cells >1:
                        t = t[:,:,None] 

                strf['H'] = np.exp(-t/self.strfTimeConstant)*(paramH['SC1']*(t/self.strfTimeConstant)**paramH['N1']/math.factorial(paramH['N1']) - \
                
                paramH['SC2'] * (t/self.strfTimeConstant)**paramH['N2']/math.factorial(paramH['N2']))
                
                strf['G'] = np.exp(-0.5*((f-paramG['f0'])/paramG['BW'])**2)* np.cos(2*np.pi*paramG['BSM']*(f-paramG['f0']))
                
                if self.num_cells >1:
                        strfG = strf['G'][:, np.newaxis, np.newaxis, np.newaxis]
                else:
                        strfG = strf['G'][:, np.newaxis, np.newaxis]
                
                strf['w1']=strfG*strf['H'][np.newaxis, :]
                strf['f']=f


                #Secondary Temporal kernel for latency derivative
                strf['H2_epsion'] = np.exp(-t/self.strfTimeConstant)
                strf['H2_gamma'] = (paramH['SC1']*(t/self.strfTimeConstant)**paramH['N1']/math.factorial(paramH['N1']) - paramH['SC2'] * (t/self.strfTimeConstant)**paramH['N2']/math.factorial(paramH['N2']))
                strf['H2'] = strf['H2_epsion']*((paramH['SC1']*paramH['N1']*(-t/self.strfTimeConstant**2))*(t/self.strfTimeConstant)**(paramH['N1']-1)/math.factorial(paramH['N1']) - (paramH['SC2']*paramH['N2']*(-t/self.strfTimeConstant**2)) * (t/self.strfTimeConstant)**(paramH['N2']-1)/math.factorial(paramH['N2'])) + strf['H2_epsion']*strf['H2_gamma']*(t/self.strfTimeConstant**2)
                strf['w2'] = strfG*strf['H2'][np.newaxis, :]
                
                return strf



        def linFwd_Junzi(self, strf, stim):
                '''
                
                Forward propagation through generalized linear model.
                
                Takes a generalized linear model data dict strf together with a stimulus matrix stim 
                and forward propagates the inputs through the network to generate a vector 
                resp_strf of output vectors.  Can also return a, the vector of 
                model responses generated before application of the output nonlinearity.
        
                parameters:
                        strf: a strf stucture
                        stim: stim matrix
                
                
                returns:
                        strf: unmodified strf structure
                        resp_strf_dat: Nx1 vector of model outputs elicited by [strf] structure by stimuli identified in datIdx.
                        a: Nx1 vector of response before output nonlinearity.

                '''
                samplesize = stim.shape[0]
                xp = cp if cp is not None else np
                stim_x = xp.asarray(stim)
                w1_x = xp.asarray(strf['w1'])
                w2_x = xp.asarray(strf['w2'])
                b1_x = xp.asarray(strf['b1'])
                a = xp.zeros((samplesize, self.num_cells, self.batch_size))
                
                #if self.num_cells >1:
                #        stim = stim[None,:,:]

                for ti in range(len(strf['delays'])):
                        #print(ti)
                        #at = np.matmul(stim, strf['w1'][:,ti][:, np.newaxis]) 

                        #With learning

                        

                        #at = np.matmul(stim, np.squeeze(strf['w1'][:,ti][:, np.newaxis])) 
                        #if self.num_cells >1:
                        #        at = xp.einsum('tf,fhb->thb', stim_x, xp.squeeze(w1_x[:,ti])) #This is the same as the above line but with broadcasting and should be faster.
                        #else:
                        w1_slice = w1_x[:, ti]
                        if w1_slice.ndim == 2:
                                w1_slice = w1_slice[:, None, :]
                        else:
                                w1_slice = np.squeeze(w1_slice[:, np.newaxis])
                        at = xp.einsum('tf,fhb->thb', stim_x, w1_slice)

                        thisshift = int(strf['delays'][ti])
                        if thisshift>=0: 
                                a[thisshift:] = a[thisshift:] + at[:len(at)-thisshift]
                        else: 
                                offset = thisshift%samplesize
                                a[:offset] = a[:offset] + at[-thisshift:]

                a = a + b1_x
                if strf['outputNL'] == 'linear':   # Linear outputs 
                        resp_strf = a                
                else:
                        raise Exception('Unknown activation function ', strf['outputNL'])

                # mask for nonvalid frames
                nanmask = np.asarray(strf['delays'])%(stim.shape[0]+1)
                nanmask = nanmask[nanmask!=0].astype(np.int64) #no mask for delay 0
                nanmask_x = xp.asarray(nanmask)
                a[nanmask_x] = 0
                resp_strf[nanmask_x] = 0

                #############################################
                #Copy same code as above for the conv deriv.#
                #############################################
                #a_deriv = np.zeros((samplesize, self.batch_size))
                a_deriv = xp.zeros((samplesize, self.num_cells, self.batch_size))
                for ti in range(len(strf['delays'])):
                        #at = np.matmul(stim, strf['w1'][:,ti][:, np.newaxis]) 

                        #With learning
                        #at_deriv = np.matmul(stim, np.squeeze(strf['w2'][:,ti][:, np.newaxis])) 
                        #at_deriv = xp.einsum('tf,fhb->thb', stim_x, xp.squeeze(w2_x[:,ti])) #This is the same as the above line but with broadcasting and should be faster.
                        w2_slice = w2_x[:, ti]
                        if w2_slice.ndim == 2:
                                w2_slice = w2_slice[:, None, :]
                        else:
                                w2_slice = np.squeeze(w2_slice[:, np.newaxis])
                        at_deriv = xp.einsum('tf,fhb->thb', stim_x, w2_slice)

                        thisshift = int(strf['delays'][ti])
                        if thisshift>=0: 
                                a_deriv[thisshift:] = a_deriv[thisshift:] + at_deriv[:len(at_deriv)-thisshift]
                        else: 
                                offset = thisshift%samplesize
                                a_deriv[:offset] = a_deriv[:offset] + at_deriv[-thisshift:]

                a_deriv = a_deriv + b1_x
                if strf['outputNL'] == 'linear':   # Linear outputs 
                        resp_strf_deriv = a_deriv                
                else:
                        raise Exception('Unknown activation function ', strf['outputNL'])

                # mask for nonvalid frames
                nanmask = np.asarray(strf['delays'])%(stim.shape[0]+1)
                nanmask = nanmask[nanmask!=0].astype(np.int64) #no mask for delay 0
                nanmask_x = xp.asarray(nanmask)
                a_deriv[nanmask_x] = 0
                resp_strf_deriv[nanmask_x] = 0

                if xp is cp:
                        return strf, cp.asnumpy(resp_strf), cp.asnumpy(resp_strf_deriv), cp.asnumpy(a)

                return strf, resp_strf, resp_strf_deriv, a


        def STRFconvolve(self, strf, stim_spec, mean_rate):
                '''
                Parameters:
                        strf
                        stim_spec
                        mean_rate
                returns:
                        onset_rate:    firing rate as a function of time for onset
                        offset_rate:    firing rate as a function of time for offset

                '''
                t = strf['t']      
                _, frate, frate_deriv,_ = self.linFwd_Junzi(strf, stim_spec) #strfFwd_Junzi(strf)
                frate = frate*mean_rate

                frate_deriv = frate_deriv*mean_rate
                
                # offset rate
                offset_rate = -frate + np.max(frate,axis=0)*0.75 #-frate + max(frate)*0.6;

                offset_rate_deriv = -frate_deriv + np.max(frate_deriv,axis=0)*0.3 #-frate + max(frate)*0.6

                #print(np.shape(offset_rate))

                #firstneg = np.where(offset_rate <= 0)[0][0] #find(offset_rate <= 0,1,'first')

                #Plot these out to make sure they look right
                
                # offset rate
                #offset_rate[:firstneg] = 0
                for k in range(self.batch_size):
                        for cc in range(self.num_cells):
                                firstneg = np.where(offset_rate[:,cc,k] <= 0)[0][0]
                                offset_rate[:firstneg,cc,k] = 0
                                offset_rate_deriv[:firstneg,cc,k] = 0
                offset_rate[offset_rate < 0] = 0
                offset_rate_deriv[offset_rate_deriv < 0] = 0

                # onset rate
                onset_rate = frate
                onset_rate_deriv = frate_deriv
                onset_rate[onset_rate < 0] = 0
                onset_rate_deriv[onset_rate_deriv < 0] = 0

                return onset_rate, offset_rate, onset_rate_deriv, offset_rate_deriv


