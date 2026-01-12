import set_options
import declarations
from BuildFile import Forwards_Method, Compile_Solve, calculate_loss
import numpy as np
import time
import Update_params
import InitParams
from scipy.io import loadmat, savemat
import yaml
import os
from strf_handler import call_strfs
from input_handler import call_inputs


class runSimulation(object):

    gen_strfs_toggle = 0  #Toggle generating the STRFs

    #Run STRF
    if gen_strfs_toggle == 1:
        call_strfs()

    #Set options
    opts = set_options.options()
    #Declare architecture
    arch = declarations.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method.Euler_Compiler(arch[0],arch[1],arch[2],opts)
    #Compile a solve file (python or c++)
    solve_file = Compile_Solve.solve_file_generator(solve_file_body = file_body_forwards, cpp_gen = 1)
    from BuildFile import generated_solve_file

    #Initilialize parameters
    num_params = 17
    batch_size = opts['N_batch']
    #Learnign rate currently set to 0.05 abs(p). Update once we have multiple parameters
    p,lr = InitParams.pinit(batch_size,num_params)
    
    losses = []
    param_tracker = []

    best_loss = 1e32

    best_output = []

    best_params = []

    start = time.perf_counter()

    #Bring in the data
    filename = f"C:/Users/isaac/OneDrive/Documents/GitHub/MouseSpatialGrid/LearningModels/BridgingMatlabtoPython/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{7}contra.mat"
    data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    data = data[:,:,None]

    spks = call_inputs(batch_size)
    on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
    off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
    noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))


    for epoch in range(1):

        

        #Make it so that you don't have to supply data if you are not running gradients
        output = generated_solve_file.solve_run(on_spks,off_spks,noise,data,p) #Python Verison to build

        #out_grads,loss = calculate_loss.calculate(output,grads)
        #losses.append(loss)
        #param_tracker.append(p)

        #print(f'L2 loss : {np.mean(loss[0]):.2f}  -:-  MSE loss : {np.mean(loss[1]):.2f} ---- Epoch: {epoch}')

        #m, v, p, t = Update_params.adam_update(m, v, p, t, beta1, beta2, lr, eps, out_grads)

        # best_loss_idx = np.argmin(np.array(loss)[1,:])

        # if np.array(loss)[1,best_loss_idx] < best_loss:
        #     best_output = output[best_loss_idx,:,:,:]
        #     best_loss = np.array(loss)[1,best_loss_idx]
        #     best_params = p[:,best_loss_idx]


    elapsed = time.perf_counter() - start
    print(f"{elapsed*1000:.2f} ms")

    #savemat("output_compressed_500.mat", {"output": output, "losses":losses, "params" : param_tracker,  "best_loss" : np.asarray(best_loss, dtype=np.float32),"best_output" : np.asarray(best_output, dtype=np.float32),"best_params" : np.asarray(best_params, dtype=np.float32)}, do_compression=True)

if __name__ == "__main__":


    run_sim = runSimulation()
