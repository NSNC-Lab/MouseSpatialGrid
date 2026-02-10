import numpy as np
import cupy as cp
import scipy.io
from scipy.io import loadmat
from scipy.signal import lfilter
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.io import loadmat, savemat

def calculate(forwards_output, grads,data):

        # -- Constants
        dt = 0.1 #ms
        dts = dt/1000 #in seconds
        bin_width = 200 # binwidth/10 = ms ex. binwidth 200 = 20 ms

        # -- Load in data

        forwards_out = cp.squeeze(cp.asarray(forwards_output, dtype=cp.float32))
        grads = cp.asarray(grads, dtype=cp.float32)
        data = cp.asarray(data, dtype=forwards_out.dtype)
        data = cp.transpose(data, (0, 3, 1, 2))

        # -- L2 Loss & Deriv Vectorized

        diff = forwards_out - data

        #print(np.shape(diff))

        L2_loss_avg = cp.mean(cp.sum(diff * diff, axis=-1), axis=-1)
        L2_deriv_avg = 2.0 * cp.mean(cp.sum(diff, axis=-1),axis=-1)

        # -- PSTH Average
        num_bins, remainder = divmod(cp.shape(data)[-1], bin_width) 
        forwards_out_r = forwards_out[:,:,:,remainder:]
        data_r = data[:,:,:,remainder:]

        forwards_out_reshaped = forwards_out_r.reshape((cp.shape(forwards_out_r)[0],cp.shape(forwards_out_r)[1],cp.shape(forwards_out_r)[2],num_bins,bin_width))
        data_reshaped = data_r.reshape((cp.shape(data_r)[0],cp.shape(data_r)[1],cp.shape(data_r)[2],num_bins,bin_width))

        #print(np.shape(forwards_out_reshaped))
        forwards_out_hist = cp.sum(cp.sum(forwards_out_reshaped,axis=-1),axis=-2)

        #print(np.shape(forwards_out_hist))

        data_hist = cp.sum(cp.sum(data_reshaped,axis=-1),axis=-2)

        import matplotlib.pyplot as plt

        diff = forwards_out_hist - data_hist

        PSTH_loss_avg = cp.sum(diff * diff, axis=-1)

        PSTH_deriv_avg = 2.0 * cp.sum(diff, axis=-1)

        #print(np.shape(PSTH_deriv_avg))
        #print(np.shape(grads))
        
        #out_grad = cp.squeeze(PSTH_deriv_avg[:,None,:,None] * grads)

        return cp.asnumpy(L2_deriv_avg), [cp.asnumpy(L2_loss_avg), cp.asnumpy(PSTH_loss_avg)]