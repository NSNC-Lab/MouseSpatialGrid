import set_options
import declarations
import declarations_output_inhibit
from BuildFile import Forwards_Method_cupy, Compile_Solve, calculate_loss
import numpy as np
import time
from datetime import datetime
import Update_params
import InitParams
from scipy.io import loadmat, savemat
import yaml
import os
from strf_handler import call_strfs
from input_handler import call_inputs
import matplotlib.pyplot as plt

class runSimulation(object):

    #Remove the outer loop later. Just for testing purposes.
    #for holder_thing in range(10):

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
    #arch = declarations.Declare_Architecture(opts)
    arch = declarations_output_inhibit.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method_cupy.Euler_Compiler(arch[0],arch[1],arch[2],opts)
    #Compile a solve file (python or c++)


    #Note!!!! This is commented out. It cannot be left this way.
    solve_file = Compile_Solve.solve_file_generator(solve_file_body = file_body_forwards, cpp_gen = 1)
    from BuildFile import generated_solve_file


    ############
    #- Move the data loading to a seperate file and make it toggleable

    #Calculate the approximate spontaneous firing rate before fitting.
    from pathlib import Path

    # --- paths (MATLAB: cd(userpath); cd('../GitHub/.../OliverDataPlotting')) ---
    userpath = Path(r"C:\Users\ipboy\Documents\MATLAB")  # <-- change this
    plot_dir = userpath / "../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting"
    plot_dir = plot_dir.resolve()
    os.chdir(plot_dir)

    # --- load MAT files ---
    mat = loadmat("all_units_info_with_polished_criteria_modified_perf.mat",variable_names=["all_data"],squeeze_me=True,struct_as_record=False)
    all_data = mat["all_data"]  # numpy array of MATLAB structs

    n = 7  # MATLAB is 1-based
    unit = all_data[n - 1]

    spike_times = unit.ctrl_tar1_timestamps

    # Ensure we have "(:,1)" behavior: take first column if 2D
    if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
        spike_times = spike_times[:, 0]

    #pre_to_aprx_3s_holder = np.zeros((10,39801))


    
    pre_to_aprx_3s_holder = np.zeros((10,29801))
    no_pre_holder = np.zeros((10,29801))

    pre_zeros_list = []
    post_zeros_list = []
    for k in range(10):  # MATLAB: 1:10
        times = np.asarray(spike_times[k]).squeeze()
        pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98) & (times >= 0)]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
        #pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98)]*(1000/opts['dt'])+(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
        #no_pre_holder[k,np.round(times[times < 0]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices

        pre_zeros_list.append(times[times < 0])
        post_zeros_list.append(times[times >= 0])
 
    pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
    FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10

    print(FR)

    # -- Load in data
    #filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{n}contra.mat"
    #data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    
    #data = data[:,:,None]
    data = pre_to_aprx_3s_holder[:,:,None]

    num_params = 5
    batch_size = opts['N_batch']
    #Learnign rate currently set to 0.05 abs(p). Update once we have multiple parameters
    p,lr = InitParams.pinit(batch_size,num_params)


    spks = call_inputs(FR,batch_size)
    on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
    off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
    noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))


    #print('p1')
    #print(p[13:17,1:5])

    savemat("compare3.mat", {"data": data, "forwards_out":on_spks, "noise": noise}, do_compression=True)

    #Grad Params
    beta1, beta2 = 0.5, 0.9995 
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

    

    for epoch in range(1):

        #spks = call_inputs(p,batch_size)
        #on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
        #off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
        #noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))

    

        #print(p)

        #Make it so that you don't have to supply data if you are not running gradients
        output, grads, on_track, off_track, loss_holder = generated_solve_file.solve_run(on_spks,off_spks,noise,data,p) #Python Verison to build
        
        #print(grads)

        #print(np.shape(on_track))

        #print(grads)
        #plt.figure()
        #plt.plot(np.squeeze(np.sum(np.sum(on_track,axis=0),axis=0)))
        #plt.show()

        #plt.figure()
        #plt.plot(np.squeeze(np.sum(np.sum(off_track,axis=0),axis=0)))
        #plt.show()

        #Calculate loss

        #print(grads)
        #print(np.shape(grads))

        out_grads,loss = calculate_loss.calculate(output,grads,data)


        grads = np.squeeze(grads)

        #print(np.shape(out_grads))

        losses.append(loss)
        param_tracker.append(p)

        #print(np.max(out_grads))
        #print(np.min(out_grads))

        #print('grads')
        #print(np.shape(out_grads))

        #print(np.shape(loss))
        #print(np.shape(loss[0]))

        print(f'L2 loss : {np.mean(loss[0]):.2f}  -:-  MSE loss : {np.mean(loss[1]):.2f} ---- Epoch: {epoch}')


        best_loss_idx = np.argmin(np.array(loss)[1,:])

        if np.array(loss)[1,best_loss_idx] < best_loss:
            best_output = output[best_loss_idx,:,:,:]
            best_loss = np.array(loss)[1,best_loss_idx]
            best_params = p[:,best_loss_idx]

        # print(out_grads)

        #Use adam optimizer on grads

        #print('p2')
        #print(p[13:17,1:5])
        #print('grad shape')
        #print(np.shape(grads))
        
        #Using grads here for of e-prop
        m, v, p, t = Update_params.adam_update(m, v, p, t, beta1, beta2, lr, eps, grads)
        #m, v, p, t = Update_params.adam_update(m, v, p, t, beta1, beta2, lr, eps, out_grads)

        #print('p3')
        #print(p[13:17,1:5])

        #print(p[0][0])

        #print(np.shape(loss))
        #print(np.shape(p))

        


    elapsed = time.perf_counter() - start
    print(f"{elapsed:.2f} s")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    savemat(f"output_compressed_Eprop_{timestamp}.mat", {"output": output, "losses":losses, "params" : param_tracker,  "best_loss" : np.asarray(best_loss, dtype=np.float32),"best_output" : np.asarray(best_output, dtype=np.float32),"best_params" : np.asarray(best_params, dtype=np.float32)}, do_compression=True)

    # # ============== PARAMETER EVOLUTION PLOT ==============
    # # param_tracker is a list of arrays, each with shape (1, 100) or (num_params, batch_size)
    # # Stack into array: (num_epochs, num_params, batch_size)
    # param_array = np.array(param_tracker)
    # print(f"param_array shape: {param_array.shape}")
    
    # fig_params, axes = plt.subplots(2, 2, figsize=(14, 10))
    # fig_params.suptitle('Parameter Evolution Across Epochs (All 100 Batch Members)', fontsize=14)
    
    # # Plot 1: Each parameter trajectory over epochs (all 100 batch members)
    # ax = axes[0, 0]
    # epochs_x = np.arange(len(param_tracker))
    # for batch_idx in range(param_array.shape[2]):
    #     ax.plot(epochs_x, param_array[:, 0, batch_idx], alpha=0.3, linewidth=1)
    # ax.set_xlabel('Epoch')
    # ax.set_ylabel('Parameter Value')
    # ax.set_title('All 100 Parameter Trajectories')
    # ax.grid(True, alpha=0.3)
    
    # # Plot 2: Parameter distribution at each epoch (box plot style)
    # # ax = axes[0, 1]
    # # param_at_epochs = [param_array[e, 0, :] for e in range(len(param_tracker))]
    # # bp = ax.boxplot(param_at_epochs, positions=epochs_x)
    # # ax.set_xlabel('Epoch')
    # # ax.set_ylabel('Parameter Value')
    # # ax.set_title('Parameter Distribution per Epoch')
    # # ax.grid(True, alpha=0.3)
    
    # # Plot 3: Parameter changes (delta) per epoch
    # ax = axes[1, 0]
    # if len(param_tracker) > 1:
    #     for batch_idx in range(param_array.shape[2]):
    #         deltas = np.diff(param_array[:, 0, batch_idx])
    #         ax.plot(epochs_x[1:], deltas, alpha=0.3, linewidth=1)
    #     ax.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Zero change')
    #     ax.set_xlabel('Epoch')
    #     ax.set_ylabel('Parameter Change (Δp)')
    #     ax.set_title('Parameter Step Size per Epoch')
    #     ax.legend()
    # ax.grid(True, alpha=0.3)
    
    # # Plot 4: Initial vs Final parameter values (scatter)
    # ax = axes[1, 1]
    # initial_params = param_array[0, 0, :]
    # final_params = param_array[-1, 0, :]
    # ax.scatter(initial_params, final_params, alpha=0.5, s=30)
    # # Add diagonal line for reference
    # min_val = min(initial_params.min(), final_params.min())
    # max_val = max(initial_params.max(), final_params.max())
    # ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='No change line')
    # ax.set_xlabel('Initial Parameter Value')
    # ax.set_ylabel('Final Parameter Value')
    # ax.set_title('Initial vs Final Parameters')
    # ax.legend()
    # ax.grid(True, alpha=0.3)
    
    # plt.tight_layout()
    # plt.savefig('parameter_evolution.png', dpi=150, bbox_inches='tight')
    # plt.show()
    # print("Saved parameter evolution plot to parameter_evolution.png")

if __name__ == "__main__":


    run_sim = runSimulation()