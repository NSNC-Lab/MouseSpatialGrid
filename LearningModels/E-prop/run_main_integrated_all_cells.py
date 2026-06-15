import set_options
import declarations
import declarations_output_inhibit
from BuildFile import Forwards_Method_cupy, Compile_Solve, calculate_loss_all_cells
import numpy as np
import time
from datetime import datetime
import Update_params_all_cells
import InitParams
import InitParamsAllCells
from scipy.io import loadmat, savemat
import yaml
import os
from pathlib import Path
from strf_handler import call_strfs
from input_handler import call_inputs
import matplotlib.pyplot as plt


CHECKPOINT_EVERY_EPOCHS = 1
CHECKPOINT_PREFIX = "checkpoint_Eprop_All_cells"
TOTAL_EPOCHS = 15
#RESUME_CHECKPOINT_PATH = "CW:\\Users\\ipboy\\Documents\\GitHub\\ModelingEffort\\Multi-Channel\\Plotting\\OliverDataPlotting\\checkpoint_Eprop_All_cells_epoch_0230.mat"  # Example: r"checkpoint_Eprop_All_cells_latest.mat"
RESUME_CHECKPOINT_PATH = ""

def save_training_checkpoint(
    epoch_one_indexed,
    p,
    m,
    v,
    t,
    losses,
    best_output_per_cell=None,
    best_loss_per_cell=None,
    best_batch_id_per_cell=None,
    best_epoch_per_cell=None,
    best_params_per_cell=None,
    prefix=CHECKPOINT_PREFIX,
):
    checkpoint_payload = {
        "epoch": np.asarray(epoch_one_indexed, dtype=np.int32),
        "params": np.asarray(p, dtype=np.float32),
        "adam_m": np.asarray(m, dtype=np.float32),
        "adam_v": np.asarray(v, dtype=np.float32),
        "adam_t": np.asarray(t, dtype=np.int32),
        "losses_so_far": np.asarray(losses, dtype=np.float32),
    }

    if best_output_per_cell is not None:
        checkpoint_payload.update(
            {
                "best_output_per_cell": np.asarray(best_output_per_cell, dtype=np.int8),
                "best_loss_per_cell": np.asarray(best_loss_per_cell, dtype=np.float32),
                "best_batch_id_per_cell": np.asarray(best_batch_id_per_cell, dtype=np.int32),
                "best_epoch_per_cell": np.asarray(best_epoch_per_cell, dtype=np.int32),
                "best_params_per_cell": np.asarray(best_params_per_cell, dtype=np.float32),
            }
        )

    checkpoint_path = Path(f"{prefix}_epoch_{epoch_one_indexed:04d}.mat")
    latest_checkpoint_path = Path(f"{prefix}_latest.mat")
    latest_epoch_txt_path = Path(f"{prefix}_latest_epoch.txt")

    savemat(str(checkpoint_path), checkpoint_payload, do_compression=True)
    savemat(str(latest_checkpoint_path), checkpoint_payload, do_compression=True)
    latest_epoch_txt_path.write_text(f"{epoch_one_indexed}\n", encoding="utf-8")

    print(
        f"Checkpoint saved at epoch {epoch_one_indexed}: "
        f"{checkpoint_path} (latest: {latest_checkpoint_path})"
    )


def load_training_checkpoint(checkpoint_path, expected_shape):
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    checkpoint = loadmat(str(checkpoint_path), squeeze_me=True, struct_as_record=False)

    required = ["epoch", "params", "adam_m", "adam_v", "adam_t"]
    missing = [key for key in required if key not in checkpoint]
    if missing:
        raise KeyError(
            f"Checkpoint missing required key(s): {missing}. "
            f"Available keys: {[k for k in checkpoint.keys() if not k.startswith('__')]}"
        )

    p_loaded = np.asarray(checkpoint["params"], dtype=np.float64)
    m_loaded = np.asarray(checkpoint["adam_m"], dtype=np.float64)
    v_loaded = np.asarray(checkpoint["adam_v"], dtype=np.float64)
    t_loaded = int(np.asarray(checkpoint["adam_t"]).squeeze())
    epoch_loaded = int(np.asarray(checkpoint["epoch"]).squeeze())

    if p_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint params shape {p_loaded.shape} does not match expected {expected_shape}."
        )
    if m_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint adam_m shape {m_loaded.shape} does not match expected {expected_shape}."
        )
    if v_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint adam_v shape {v_loaded.shape} does not match expected {expected_shape}."
        )

    best_state = {}
    optional_best_keys = [
        "best_output_per_cell",
        "best_loss_per_cell",
        "best_batch_id_per_cell",
        "best_epoch_per_cell",
        "best_params_per_cell",
    ]
    for key in optional_best_keys:
        if key in checkpoint:
            best_state[key] = checkpoint[key]

    return p_loaded, m_loaded, v_loaded, t_loaded, epoch_loaded, best_state


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

    num_cells = 220

    #num_cells, batch, trials, channels, timecourse
    #batch,trials,channels,timecourse

    #Set options
    opts = set_options.options()
    #Declare architecture
    #arch = declarations.Declare_Architecture(opts)
    #arch = declarations_output_inhibit.Declare_Architecture(opts)
    #Warning! Make sure you know what declarations you are pulling from
    arch = declarations.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method_cupy.Euler_Compiler(arch[0],arch[1],arch[2],opts,num_cells)
    #Compile a solve file (python or c++)


    #Note!!!! This is commented out. It cannot be left this way.
    solve_file = Compile_Solve.solve_file_generator(solve_file_body = file_body_forwards, cpp_gen = 1)
    from BuildFile import generated_solve_file


    ############
    #- Move the data loading to a seperate file and make it toggleable

    #Calculate the approximate spontaneous firing rate before fitting.
    # --- paths (MATLAB: cd(userpath); cd('../GitHub/.../OliverDataPlotting')) ---
    userpath = Path(r"C:\Users\ipboy\Documents\MATLAB")  # <-- change this
    plot_dir = userpath / "../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting"
    plot_dir = plot_dir.resolve()
    os.chdir(plot_dir)

    # --- load MAT files ---
    mat = loadmat("all_units_info_with_polished_criteria_modified_perf.mat",variable_names=["all_data"],squeeze_me=True,struct_as_record=False)
    all_data = mat["all_data"]  # numpy array of MATLAB structs

    #cell_storage = []

    

    no_pre_holder = np.zeros((num_cells,10,29801))
    FRs = []

    #Get peak indicies
    peak_indices = []
    for n in range(num_cells):
        vals = []
        for m in range(4):
            vals2 = 0
            for k in range(10):
                if len(np.shape(all_data[n].ctrl_tar1_timestamps[k, m])) != 0:  #This just skips areas where there are no spikes
                    vals2 += np.shape(all_data[n].ctrl_tar1_timestamps[k, m])[0]
            vals.append(vals2)

        idx = np.where(vals == np.max(np.array(vals)))[0][0]
        vals = []
        peak_indices.append(idx)
       



    for k in range(num_cells):
        spike_times = all_data[k].ctrl_tar1_timestamps
        spike_times = spike_times[:, peak_indices[k]]  #Grab -90 degrees  #4.1.2026 -> Instead of grabbing -90 grab the indicy corresponding to the peak of the tuning curve.
        #cell_storage.append(spike_times)

        pre_zeros_list = []
        post_zeros_list = []

        for m in range(10):  # MATLAB: 1:10
            times = np.asarray(spike_times[m]).squeeze()
            no_pre_holder[k,m,np.round(times[(times < 2.98) & (times >= 0)]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
            #pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98)]*(1000/opts['dt'])+(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
            #no_pre_holder[k,np.round(times[times < 0]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices

            pre_zeros_list.append(times[times < 0])
            #post_zeros_list.append(times[times >= 0])

        pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
        FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10
        FRs.append(FR)

    #n = 7  # MATLAB is 1-based
    #unit = all_data[n - 1]

    #spike_times = unit.ctrl_tar1_timestamps

    # Ensure we have "(:,1)" behavior: take first column if 2D
    #if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
    #    spike_times = spike_times[:, 0]

    #pre_to_aprx_3s_holder = np.zeros((10,39801))


    
    #pre_to_aprx_3s_holder = np.zeros((10,29801))
    

    # pre_zeros_list = []
    # post_zeros_list = []
    # for k in range(10):  # MATLAB: 1:10
    #     times = np.asarray(spike_times[k]).squeeze()
    #     pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98) & (times >= 0)]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
    #     #pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98)]*(1000/opts['dt'])+(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
    #     #no_pre_holder[k,np.round(times[times < 0]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices

    #     pre_zeros_list.append(times[times < 0])
    #     post_zeros_list.append(times[times >= 0])
 
    # pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
    # FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10

    # print(FR)

    # -- Load in data
    #filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{n}contra.mat"
    #data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    
    #data = data[:,:,None]
    data = no_pre_holder[:,:,:,None]

    num_params = 8
    batch_size = opts['N_batch']
    #Learnign rate currently set to 0.05 abs(p). Update once we have multiple parameters
    p,lr = InitParamsAllCells.pinit(batch_size,num_params,num_cells)
    #p,lr = InitParams.pinit(batch_size,num_params)


    #spks = call_inputs(num_cells,FRs,batch_size)
    #on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
    #off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
    #noise = np.transpose(spks['noise_masker_None_target_0'],(0,1,4,2,3))

    #Repeat along the first axis to match the data shape #Note! You might be able to save some space by just using this across all tials... all trials are giong to have the bigger size anyways through so it might not matter much.
    #On second thought this should broadcast correctly.
    #on_spks= np.repeat(on_spks[None, ...], 220, axis=0) 
    #off_spks= np.repeat(off_spks[None, ...], 220, axis=0) 
    #noise= np.repeat(noise[None, ...], 220, axis=0) 

    #print('p1')
    #print(p[13:17,1:5])

    # Saving this full `noise` tensor to a v5 MAT-file will often exceed MATLAB's ~2GB limit.
    # Use the Python raster plot below instead of exporting the full array.

    #if plot_noise_toggle == 1:

    plot_noise_toggle = 1
    plot_noise_cell = 9      # MATLAB-style indexing (cell 1..num_cells)
    plot_noise_batch = 2     # MATLAB-style indexing (batch 1..N_batch)
    plot_noise_channel = 1   # MATLAB-style indexing (channel 1..chans)


    #from plot_noise_raster import plot_noise_raster
    #plot_noise_raster(noise,cell=plot_noise_cell,batch=plot_noise_batch,channel=plot_noise_channel,dt_ms=opts['dt'],one_indexed=True,)
    #Grad Params
    beta1, beta2 = 0.5, 0.9995 
    eps = 1e-6
    #lr = 0.5

    #Init mvt
    m = np.zeros((num_params,num_cells,batch_size))
    v = np.zeros((num_params,num_cells,batch_size))
    t = 0
    

    losses = []
    param_tracker = []
    start_epoch = 0

    best_output_per_cell = np.zeros((num_cells, opts["N_trials"], opts["N_channels"], opts["sim_len"]), dtype=np.int8)
    best_loss_per_cell = np.full((num_cells,), np.inf, dtype=np.float32)
    best_batch_id_per_cell = np.zeros((num_cells,), dtype=np.int32)
    best_epoch_per_cell = np.zeros((num_cells,), dtype=np.int32)
    best_params_per_cell = np.zeros((num_params, num_cells), dtype=np.float32)

    if RESUME_CHECKPOINT_PATH:
        expected_shape = p.shape
        p, m, v, t, loaded_epoch, checkpoint_best_state = load_training_checkpoint(
            checkpoint_path=RESUME_CHECKPOINT_PATH,
            expected_shape=expected_shape,
        )
        start_epoch = loaded_epoch
        if checkpoint_best_state:
            loaded_best_output_per_cell = np.asarray(
                checkpoint_best_state.get("best_output_per_cell", best_output_per_cell),
                dtype=np.int8,
            )
            if loaded_best_output_per_cell.ndim == 3:
                loaded_best_output_per_cell = loaded_best_output_per_cell[:, :, None, :]
            best_output_per_cell = loaded_best_output_per_cell
            best_loss_per_cell = np.asarray(
                checkpoint_best_state.get("best_loss_per_cell", best_loss_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_batch_id_per_cell = np.asarray(
                checkpoint_best_state.get("best_batch_id_per_cell", best_batch_id_per_cell),
                dtype=np.int32,
            ).reshape(num_cells)
            best_epoch_per_cell = np.asarray(
                checkpoint_best_state.get("best_epoch_per_cell", best_epoch_per_cell),
                dtype=np.int32,
            ).reshape(num_cells)
            best_params_per_cell = np.asarray(
                checkpoint_best_state.get("best_params_per_cell", best_params_per_cell),
                dtype=np.float32,
            )
        print(
            f"Resumed from checkpoint '{RESUME_CHECKPOINT_PATH}' "
            f"(saved after epoch {loaded_epoch})."
        )

    best_loss = 1e32

    best_output = []

    best_params = []

    start = time.perf_counter()

    start_epoch = 0

    for epoch in range(start_epoch, TOTAL_EPOCHS):

        print(range(start_epoch, TOTAL_EPOCHS))
        #spks = call_inputs(p,batch_size)
        #on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
        #off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
        #noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))

        target_dict = call_strfs(p,batch_size, num_cells)


        spks = call_inputs(num_cells,FRs,batch_size,target_dict)
        
        #Organize as batch,trial,channel,time 
        
        on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(3,4,2,1,0))
        off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(3,4,2,1,0))
        rate_on = spks[f'locs_masker_None_target_0_on'][f'stimulus_0_rate']
        rate_off = spks[f'locs_masker_None_target_0_off'][f'stimulus_0_rate']
        rate_on_deriv = spks[f'locs_masker_None_target_0_on'][f'stimulus_0_rate_deriv']
        rate_off_deriv = spks[f'locs_masker_None_target_0_off'][f'stimulus_0_rate_deriv']
        noise = np.transpose(spks['noise_masker_None_target_0'],(0,1,4,3,2))
    

        #print(p)

        #Make it so that you don't have to supply data if you are not running gradients
        output, grads, on_track, off_track, loss_holder = generated_solve_file.solve_run(on_spks,off_spks,noise,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p) #Python Verison to build
        
        #generated_solve_file.solve_run(on_spks,off_spks,noise,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p) #Python Verison to build
        
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

        out_grads,loss = calculate_loss_all_cells.calculate(output,grads,data)


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

        psth_loss = np.asarray(loss[1], dtype=np.float32)
        best_batch_idx_this_epoch = np.argmin(psth_loss, axis=1)
        best_loss_this_epoch = psth_loss[np.arange(num_cells), best_batch_idx_this_epoch]
        improved_cells = best_loss_this_epoch < best_loss_per_cell

        if np.any(improved_cells):
            improved_idx = np.where(improved_cells)[0]
            chosen_batch_idx = best_batch_idx_this_epoch[improved_idx]

            best_output_per_cell[improved_idx, :, :, :] = np.asarray(
                output[improved_idx, chosen_batch_idx, :, :, :],
                dtype=np.int8,
            )
            best_loss_per_cell[improved_idx] = best_loss_this_epoch[improved_idx]
            best_batch_id_per_cell[improved_idx] = chosen_batch_idx.astype(np.int32) + 1
            best_epoch_per_cell[improved_idx] = epoch + 1
            best_params_per_cell[:, improved_idx] = p[:, improved_idx, chosen_batch_idx].astype(np.float32)

            print(
                f'Updated best output for {len(improved_idx)} cells '
                f'(median best PSTH loss: {np.nanmedian(best_loss_per_cell):.2f})'
            )

        # print(out_grads)

        #Use adam optimizer on grads

        #print('p2')
        #print(p[13:17,1:5])
        #print('grad shape')
        #print(np.shape(grads))
        
        #Using grads here for of e-prop
        m, v, p, t = Update_params_all_cells.adam_update(m, v, p, t, beta1, beta2, lr, eps, grads)
        #m, v, p, t = Update_params.adam_update(m, v, p, t, beta1, beta2, lr, eps, out_grads)

        epoch_one_indexed = epoch + 1
        if epoch_one_indexed % CHECKPOINT_EVERY_EPOCHS == 0:
            save_training_checkpoint(
                epoch_one_indexed=epoch_one_indexed,
                p=p,
                m=m,
                v=v,
                t=t,
                losses=losses,
                best_output_per_cell=best_output_per_cell,
                best_loss_per_cell=best_loss_per_cell,
                best_batch_id_per_cell=best_batch_id_per_cell,
                best_epoch_per_cell=best_epoch_per_cell,
                best_params_per_cell=best_params_per_cell,
            )

        #print('p3')
        #print(p[13:17,1:5])

        #print(p[0][0])

        #print(np.shape(loss))
        #print(np.shape(p))

        


    elapsed = time.perf_counter() - start
    print(f"{elapsed:.2f} s")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    #savemat(f"output_compressed_Eprop_{timestamp}.mat", {"output": output, "losses":losses, "params" : param_tracker,  "best_loss" : np.asarray(best_loss, dtype=np.float32),"best_output" : np.asarray(best_output, dtype=np.float32),"best_params" : np.asarray(best_params, dtype=np.float32)}, do_compression=True)
    savemat(
        f"output_compressed_Eprop_All_cells_{timestamp}.mat",
        {
            "output": output,
            "losses": losses,
            "params": param_tracker,
            "best_output_per_cell": best_output_per_cell,
            "best_loss_per_cell": best_loss_per_cell,
            "best_batch_id_per_cell": best_batch_id_per_cell,
            "best_epoch_per_cell": best_epoch_per_cell,
            "best_params_per_cell": best_params_per_cell,
        },
        do_compression=True,
    )

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
