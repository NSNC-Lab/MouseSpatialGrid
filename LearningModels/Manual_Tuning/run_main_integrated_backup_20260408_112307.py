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
from strf_handler import call_strfs, export_strf_temporal_kernel_svg, export_onset_offset_rate_svg
from input_handler import call_inputs
import matplotlib.pyplot as plt
from plot_data_raster_svg import make_poster_raster, make_poster_psth_overlay


class runSimulation(object):

    #Remove the outer loop later. Just for testing purposes.
    #for holder_thing in range(10):

    #Little control pannel for now (Eventually move this to a yaml file or whatever makes the most sense)

    gen_strfs_toggle = 1  #Toggle generating the STRFs
    gradients_toggle = 0  #Toggle generating the graidnets in the forwards process *Also toggles running epochs
    export_data_raster_svg_toggle = 1  # Toggle exporting data-only raster as SVG
    export_sim_raster_svg_toggle = 1  # Toggle exporting simulation raster as SVG (same style as data raster)
    export_psth_overlay_svg_toggle = 1  # Toggle exporting minimal sim/data PSTH overlay as SVG
    export_strf_temporal_kernel_svg_toggle = 1  # Toggle exporting STRF temporal kernel as SVG
    export_onset_offset_rate_svg_toggle = 1  # Toggle exporting onset/offset firing-rate SVG
    onset_offset_trial_index = 0  # Trial index used if onset/offset rate has multiple trials
    onset_offset_stimulus_index = 0  # Target stimulus index

    #Run STRF
    if gen_strfs_toggle == 1:
        call_strfs()

    if export_strf_temporal_kernel_svg_toggle == 1:
        export_strf_temporal_kernel_svg("strf_temporal_kernel.svg")

    if export_onset_offset_rate_svg_toggle == 1:
        export_onset_offset_rate_svg(
            output_svg="strf_onset_offset_single_trial.svg",
            stimulus_index=onset_offset_stimulus_index,
            trial_index=onset_offset_trial_index,
        )


    #Set options
    opts = set_options.options()
    #Declare architecture
    arch = declarations.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method.Euler_Compiler(arch[0],arch[1],arch[2],opts)
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

    n = 2  # MATLAB is 1-based
    unit = all_data[n - 1]

    spike_times = unit.ctrl_tar1_timestamps

    # Ensure we have "(:,1)" behavior: take first column if 2D
    if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
        spike_times = spike_times[:, 0]

    pre_zeros_list = []
    for k in range(10):  # MATLAB: 1:10
        times = np.asarray(spike_times[k]).squeeze()
        pre_zeros_list.append(times[times < 0])

    pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
    FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10

    # -- Load in data
    filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{n}contra.mat"
    data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    data = data[:,:,None]

    if export_data_raster_svg_toggle == 1:
        data_trials_time = np.asarray(data).squeeze()
        if data_trials_time.ndim == 1:
            data_trials_time = data_trials_time[None, :]

        if data_trials_time.shape[0] > data_trials_time.shape[1]:
            data_trials_time = data_trials_time.T

        raster_svg_out = Path.cwd() / f"picture_fit{n}_data_raster.svg"
        make_poster_raster(
            data_trials_time=data_trials_time,
            dt_ms=float(opts.get('dt', 1.0)),
            threshold=0.0,
            fig_width_in=3.4,
            fig_height_in=1.5,
            marker_size=6.5,
            line_width=1.0,
            title=f"Unit {n} Data Raster",
            output_svg=raster_svg_out,
        )
        print(f"Saved data raster SVG: {raster_svg_out}")



    num_params = 5
    batch_size = opts['N_batch']
    #Learnign rate currently set to 0.05 abs(p). Update once we have multiple parameters
    p,lr = InitParams.pinit(batch_size,num_params)


    spks = call_inputs(FR,batch_size)
    on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
    off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
    noise = np.transpose(spks['noise_masker_None_target_0'],(0,3,1,2))



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

        output = generated_solve_file.solve_run(on_spks,off_spks,noise,data,p) #Python Verison to build        

    elapsed = time.perf_counter() - start

    print(f"{elapsed:2.2f} seconds")


    print(np.shape(output))


    # Plot model output raster + data raster in the same figure
    # output: batch x trials x neurons x time
    # data: trials x time x 1 (loaded above)
    spikes = np.asarray(output)
    dt_ms = float(opts.get('dt', 1.0))

    if spikes.ndim == 4:
        batch_idx = 0
        spikes_tnt = spikes[batch_idx]  # trials x neurons x time

        if export_sim_raster_svg_toggle == 1 and spikes_tnt.ndim == 3:
            sim_trials_time = (np.sum(spikes_tnt, axis=1) > 0).astype(np.float32)
            sim_raster_svg_out = Path.cwd() / f"picture_fit{n}_sim_raster.svg"
            make_poster_raster(
                data_trials_time=sim_trials_time,
                dt_ms=float(opts.get('dt', 1.0)),
                threshold=0.0,
                fig_width_in=3.4,
                fig_height_in=1.5,
                marker_size=6.5,
                line_width=1.0,
                title=f"Unit {n} Sim Raster",
                output_svg=sim_raster_svg_out,
            )
            print(f"Saved simulation raster SVG: {sim_raster_svg_out}")

        data_trials_time = np.asarray(data).squeeze()
        if data_trials_time.ndim == 1:
            data_trials_time = data_trials_time[None, :]

        fig, (ax_out, ax_data, ax_psth) = plt.subplots(
            3, 1, sharex=True, figsize=(10, 8)
        )

        # --- output raster ---
        if spikes_tnt.ndim == 3:
            n_trials, n_neurons, _ = spikes_tnt.shape
            out_event_times = []
            for trial_idx in range(n_trials):
                for neuron_idx in range(n_neurons):
                    spike_idx = np.flatnonzero(spikes_tnt[trial_idx, neuron_idx, :])
                    out_event_times.append(spike_idx * dt_ms)

            ax_out.eventplot(out_event_times, linewidths=0.8)
            ax_out.set_ylabel('Trial / neuron')
            ax_out.set_title('Model output raster')

        # --- data raster ---
        data_event_times = []
        for trial_idx in range(data_trials_time.shape[0]):
            spike_idx = np.flatnonzero(data_trials_time[trial_idx, :] > 0)
            data_event_times.append(spike_idx * dt_ms)

        ax_data.eventplot(data_event_times, linewidths=0.8)
        ax_data.set_ylabel('Trial')
        ax_data.set_title('Data raster')

        # --- PSTH overlay (trial-averaged), binned at 20ms ---
        bin_width_ms = 20.0
        bin_width_samples = int(round(bin_width_ms / dt_ms))
        if bin_width_samples < 1:
            bin_width_samples = 1
        actual_bin_ms = bin_width_samples * dt_ms

        def bin_1d_sum(x: np.ndarray, width: int) -> np.ndarray:
            x = np.asarray(x)
            if x.ndim != 1:
                x = x.reshape(-1)
            n = (x.shape[0] // width) * width
            if n == 0:
                return np.zeros((0,), dtype=np.float32)
            return x[:n].reshape(-1, width).sum(axis=1)

        # Sim PSTH: sum across neurons per trial, then average across trials
        sim_psth_binned = None
        if spikes_tnt.ndim == 3:
            sim_trial_counts = spikes_tnt.sum(axis=1)  # trials x time
            sim_psth = sim_trial_counts.sum(axis=0)
            sim_psth_binned = bin_1d_sum(sim_psth, bin_width_samples)

        # Data PSTH: average across trials (data may be binary spikes or a continuous “picture”)
        data_psth = data_trials_time.sum(axis=0)
        data_psth_binned = bin_1d_sum(data_psth, bin_width_samples)

        t_ms_bins = (np.arange(data_psth_binned.shape[0]) + 0.5) * actual_bin_ms
        ax_psth.plot(t_ms_bins, data_psth_binned, label=f'Data PSTH ({actual_bin_ms:g} ms bins)')
        if sim_psth_binned is not None and sim_psth_binned.shape[0] > 0:
            t_ms_bins_sim = (np.arange(sim_psth_binned.shape[0]) + 0.5) * actual_bin_ms
            ax_psth.plot(t_ms_bins_sim, sim_psth_binned, label=f'Sim PSTH ({actual_bin_ms:g} ms bins)')

        if export_psth_overlay_svg_toggle == 1:
            psth_svg_out = Path.cwd() / f"picture_fit{n}_psth_overlay.svg"
            make_poster_psth_overlay(
                t_ms=t_ms_bins,
                data_psth=data_psth_binned,
                sim_psth=sim_psth_binned,
                fig_width_in=3.4,
                fig_height_in=1.5,
                title=f"Unit {n} PSTH",
                output_svg=psth_svg_out,
            )
            print(f"Saved PSTH overlay SVG: {psth_svg_out}")

        ax_psth.set_xlabel('Time (ms)')
        ax_psth.set_ylabel('Mean spikes / 20ms bin')
        ax_psth.set_title('PSTH overlay')
        ax_psth.legend(loc='upper right')

        plt.tight_layout()
        plt.show()

    diff = sim_psth_binned - data_psth_binned

    PSTH_loss_avg = np.sum(diff * diff)

    print(f"PSTH loss: {PSTH_loss_avg}")
    


if __name__ == "__main__":


    run_sim = runSimulation()