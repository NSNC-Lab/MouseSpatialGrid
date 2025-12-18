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

    #Little control pannel for now (Eventually move this to a yaml file or whatever makes the most sense)

    gen_strfs_toggle = 0  #Toggle generating the STRFs
    gradients_toggle = 0  #Toggle generating the graidnets in the forwards process *Also toggles running epochs

    #Run STRF
    if gen_strfs_toggle == 1:
        call_strfs()
    #PreProcessesing  #Note! Recheck everything once you start running multichannel inputs -- check where the gain control for the tuning curves is and make sure we arn't doing extra steps
    #Will also need two worry about how exactly you are going to parse spks once you have multiple data streams
    
    

    #batch,trials,channels,timecourse

    #Set options
    opts = set_options.options()
    #Declare architecture
    arch = declarations.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method.Euler_Compiler(arch[0],arch[1],arch[2],opts)
    #Compile a solve file (python or c++)
    solve_file = Compile_Solve.solve_file_generator(solve_file_body = file_body_forwards, cpp_gen = 1)
    from BuildFile import generated_solve_file


    ############
    #- Move the data loading to a seperate file and make it toggleable

    # -- Load in data
    filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{7}contra.mat"
    data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    data = data[:,:,None]

    num_params = 17
    batch_size = opts['N_batch']
    #Learnign rate currently set to 0.05 abs(p). Update once we have multiple parameters
    p,lr = InitParams.pinit(batch_size,num_params)

    #print('p1')
    #print(p[13:17,1:5])

    #savemat("compare.mat", {"data": data, "forwards_out":on_spks}, do_compression=True)

    #Grad Params
    beta1, beta2 = 0.99, 0.9995 
    eps = 1e-6
    #lr = 0.5

    #Init mvt
    m = np.zeros((num_params,batch_size))
    v = np.zeros((num_params,batch_size))
    t = 0
    
    losses = []
    param_tracker = []

    best_loss = 1e32

    best_output = []

    best_params = []

    start = time.perf_counter()
    for epoch in range(100):

        #spks = call_inputs(p,batch_size)
        #on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
        #off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
        #noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))

        spks = call_inputs(p,batch_size)
        on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
        off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
        noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))

        #print(p)

        #Make it so that you don't have to supply data if you are not running gradients
        output, grads = generated_solve_file.solve_run(on_spks,off_spks,noise,data,p) #Python Verison to build
        
        #print(grads)

        #Calculate loss

        #print(grads)
        #print(np.shape(grads))

        out_grads,loss = calculate_loss.calculate(output,grads)
        losses.append(loss)
        param_tracker.append(p)

        #print(np.max(out_grads))
        #print(np.min(out_grads))

        #print('grads')
        #print(np.shape(out_grads))

        #print(np.shape(loss))
        #print(np.shape(loss[0]))

        print(f'L2 loss : {np.mean(loss[0]):.2f}  -:-  MSE loss : {np.mean(loss[1]):.2f} ---- Epoch: {epoch}')

       # print(out_grads)

        #Use adam optimizer on grads

        #print('p2')
        #print(p[13:17,1:5])

        m, v, p, t = Update_params.adam_update(m, v, p, t, beta1, beta2, lr, eps, out_grads)

        #print('p3')
        #print(p[13:17,1:5])

        #print(p[0][0])

        #print(np.shape(loss))
        #print(np.shape(p))

        best_loss_idx = np.argmin(np.array(loss)[1,:])

        if np.array(loss)[1,best_loss_idx] < best_loss:
            best_output = output[best_loss_idx,:,:,:]
            best_loss = np.array(loss)[1,best_loss_idx]
            best_params = p[:,best_loss_idx]


    elapsed = time.perf_counter() - start
    print(f"{elapsed*1000:.2f} ms")

    savemat("output_compressed_500.mat", {"output": output, "losses":losses, "params" : param_tracker,  "best_loss" : np.asarray(best_loss, dtype=np.float32),"best_output" : np.asarray(best_output, dtype=np.float32),"best_params" : np.asarray(best_params, dtype=np.float32)}, do_compression=True)

if __name__ == "__main__":


    run_sim = runSimulation()