import numpy as np
import scipy.io
from scipy.io import loadmat
from scipy.signal import lfilter
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.io import loadmat, savemat

def calculate(forwards_output, grads):

        # -- Constants
        dt = 0.1 #ms
        dts = dt/1000 #in seconds
        bin_width = 200 # binwidth/10 = ms ex. binwidth 200 = 20 ms

        # -- Load in data

        forwards_out = np.squeeze(np.asarray(forwards_output, dtype=np.float32))



        #print(np.shape(forwards_out))

        filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{7}contra.mat"
        data = loadmat(filename)['picture'].astype(np.float32)[:,:,None]

        savemat("compare.mat", {"data": data, "forwards_out":forwards_out}, do_compression=True)

        #print(np.shape(data))

        data = np.transpose(data,(2,0,1)) #Transpose things to be Batch,trials,timecouse
        #forwards_out = np.transpose(forwards_out,(2,0,1))

        # -- L2 Loss & Deriv Vectorized

        diff = forwards_out - data
        L2_loss_avg = np.mean(np.sum(diff * diff, axis=-1), axis=-1)

        # -- PSTH Average
        num_bins, remainder = divmod(np.shape(data)[-1], bin_width) 

        forwards_out_r = forwards_out[:,:,remainder:]
        data_r = data[:,:,remainder:]

        forwards_out_reshaped = forwards_out_r.reshape((np.shape(forwards_out_r)[0],np.shape(forwards_out_r)[1],num_bins,bin_width))
        data_reshaped = data_r.reshape((np.shape(data_r)[0],np.shape(data_r)[1],num_bins,bin_width))

        forwards_out_hist = np.sum(np.sum(forwards_out_reshaped,axis=-1),axis=-2)
        data_hist = np.sum(np.sum(data_reshaped,axis=-1),axis=-2)

        import matplotlib.pyplot as plt

        diff = forwards_out_hist - data_hist

        PSTH_loss_avg = np.sum(diff * diff, axis=-1)

        PSTH_deriv_avg = 2.0 * np.sum(diff, axis=-1)

        #print(np.shape(PSTH_deriv_avg))
        #print(np.shape(grads))
        
        out_grad = np.squeeze(PSTH_deriv_avg[None,:,None] * grads)

        return out_grad, [L2_loss_avg, PSTH_loss_avg]