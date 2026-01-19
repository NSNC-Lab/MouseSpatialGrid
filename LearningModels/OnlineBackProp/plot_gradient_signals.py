"""
Visualize L_t × eligibility trace over time to understand gradient signal dynamics.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.chdir(r'c:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\E-prop')
sys.path.insert(0, r'c:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\E-prop')

from scipy.io import loadmat
from pathlib import Path
from input_handler import call_inputs
import set_options


def load_data():
    """Load the actual data."""
    opts = set_options.options()
    batch_size = opts['N_batch']
    
    # Load target spike data
    filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{1}contra.mat"
    data = loadmat(filename)['picture'].astype(np.float32)
    data = data[:, :, None]
    
    # Calculate FR
    userpath = Path(r"C:\Users\ipboy\Documents\MATLAB")
    plot_dir = userpath / "../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting"
    plot_dir = plot_dir.resolve()
    
    mat = loadmat(str(plot_dir / "all_units_info_with_polished_criteria_modified_perf.mat"),
                  variable_names=["all_data"], squeeze_me=True, struct_as_record=False)
    all_data = mat["all_data"]
    unit = all_data[7 - 1]
    spike_times = unit.ctrl_tar1_timestamps
    
    if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
        spike_times = spike_times[:, 0]
    
    pre_zeros_list = []
    for k in range(10):
        times = np.asarray(spike_times[k]).squeeze()
        pre_zeros_list.append(times[times < 0])
    pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
    FR = pre_zeros.size / 10
    
    # Load inputs
    spks = call_inputs(FR, batch_size)
    on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'], (2, 0, 1))
    off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'], (2, 0, 1))
    noise = np.transpose(spks['noise_masker_None_target_0'], (0, 3, 1, 2))
    
    return on_spks, off_spks, noise, data, batch_size


def run_with_tracking(on_input, off_input, noise_token, data, p):
    """
    Run simulation and track L_t × eligibility at each timestep.
    This is a simplified version of solve_run that captures the gradient signals.
    """
    
    # Key parameters
    ROn_R = 200.0
    ROn_tau = 20.0
    ROn_V_thresh = -47
    ROn_E_L = -65
    SOnOff_R = 100.0
    SOnOff_tau = 10.0
    SOnOff_V_thresh = -47
    SOnOff_E_L = -57
    
    On_ROn_ESYN = 0
    Off_ROn_ESYN = 0
    SOnOff_ROn_ESYN = -80
    Off_SOnOff_ESYN = 0
    On_SOnOff_ESYN = 0
    
    SOnOff_ROn_gSYN = 0.025
    
    # Synapse params
    On_ROn_tauR, On_ROn_tauD = 0.7, 1.5
    Off_ROn_tauR, Off_ROn_tauD = 0.7, 1.5
    On_SOnOff_tauR, On_SOnOff_tauD = 0.1, 1.0
    Off_SOnOff_tauR, Off_SOnOff_tauD = 0.1, 1.0
    SOnOff_ROn_tauR, SOnOff_ROn_tauD = 1.0, 4.5
    
    On_ROn_scale = 1.9481350796278847
    Off_ROn_scale = 1.9481350796278847
    On_SOnOff_scale = 1.2915496650148839
    Off_SOnOff_scale = 1.2915496650148839
    SOnOff_ROn_scale = 1.5368523544529802
    
    Off_SOnOff_gSYN = p[0].reshape(100, 1, 1)
    
    timesteps = 29801
    dt = 0.1
    
    # State variables (simplified)
    ROn_V = np.ones((100, 10, 1)) * ROn_E_L
    SOnOff_V = np.ones((100, 10, 1)) * SOnOff_E_L
    
    # PSC variables
    On_ROn_PSC_s = np.zeros((100, 10, 1))
    On_ROn_PSC_x = np.zeros((100, 10, 1))
    Off_ROn_PSC_s = np.zeros((100, 10, 1))
    Off_ROn_PSC_x = np.zeros((100, 10, 1))
    On_SOnOff_PSC_s = np.zeros((100, 10, 1))
    On_SOnOff_PSC_x = np.zeros((100, 10, 1))
    Off_SOnOff_PSC_s = np.zeros((100, 10, 1))
    Off_SOnOff_PSC_x = np.zeros((100, 10, 1))
    SOnOff_ROn_PSC_s = np.zeros((100, 10, 1))
    SOnOff_ROn_PSC_x = np.zeros((100, 10, 1))
    
    # Spike tracking
    On_last_spike = np.ones((100, 10, 1)) * -100
    Off_last_spike = np.ones((100, 10, 1)) * -100
    SOnOff_last_spike = np.ones((100, 10, 1)) * -100
    ROn_last_spike = np.ones((100, 10, 1)) * -100
    
    # Tracking arrays for plotting
    L_t_track = np.zeros(timesteps)
    eligibility_Off_SOnOff_track = np.zeros(timesteps)
    eligibility_On_ROn_track = np.zeros(timesteps)
    eligibility_SOnOff_ROn_track = np.zeros(timesteps)
    Lt_x_elig_Off_SOnOff_track = np.zeros(timesteps)
    Lt_x_elig_On_ROn_track = np.zeros(timesteps)
    Lt_x_elig_SOnOff_ROn_track = np.zeros(timesteps)
    ROn_V_track = np.zeros(timesteps)
    SOnOff_V_track = np.zeros(timesteps)
    psi_ROn_track = np.zeros(timesteps)
    psi_SOnOff_track = np.zeros(timesteps)
    ROn_spikes_track = np.zeros(timesteps)
    SOnOff_spikes_track = np.zeros(timesteps)
    target_track = np.zeros(timesteps)
    cumulative_grad_Off_SOnOff = np.zeros(timesteps)
    cumulative_grad_On_ROn = np.zeros(timesteps)
    
    grad_Off_SOnOff = 0.0
    grad_On_ROn = 0.0
    
    print("Running simulation with tracking...")
    
    for t_idx in range(timesteps):
        t = t_idx * dt
        
        # Simple PSC dynamics (double exponential)
        On_ROn_PSC_s += dt * (On_ROn_scale * On_ROn_PSC_x - On_ROn_PSC_s) / On_ROn_tauR
        On_ROn_PSC_x += dt * (-On_ROn_PSC_x / On_ROn_tauD)
        Off_ROn_PSC_s += dt * (Off_ROn_scale * Off_ROn_PSC_x - Off_ROn_PSC_s) / Off_ROn_tauR
        Off_ROn_PSC_x += dt * (-Off_ROn_PSC_x / Off_ROn_tauD)
        On_SOnOff_PSC_s += dt * (On_SOnOff_scale * On_SOnOff_PSC_x - On_SOnOff_PSC_s) / On_SOnOff_tauR
        On_SOnOff_PSC_x += dt * (-On_SOnOff_PSC_x / On_SOnOff_tauD)
        Off_SOnOff_PSC_s += dt * (Off_SOnOff_scale * Off_SOnOff_PSC_x - Off_SOnOff_PSC_s) / Off_SOnOff_tauR
        Off_SOnOff_PSC_x += dt * (-Off_SOnOff_PSC_x / Off_SOnOff_tauD)
        SOnOff_ROn_PSC_s += dt * (SOnOff_ROn_scale * SOnOff_ROn_PSC_x - SOnOff_ROn_PSC_s) / SOnOff_ROn_tauR
        SOnOff_ROn_PSC_x += dt * (-SOnOff_ROn_PSC_x / SOnOff_ROn_tauD)
        
        # Get inputs at this timestep
        on_in = on_input[:, t_idx, :] if t_idx < on_input.shape[1] else np.zeros((10, 1))
        off_in = off_input[:, t_idx, :] if t_idx < off_input.shape[1] else np.zeros((10, 1))
        
        # Update PSC_x on presynaptic spikes (simplified - using input directly)
        On_ROn_PSC_x += on_in.reshape(1, 10, 1)
        Off_ROn_PSC_x += off_in.reshape(1, 10, 1)
        On_SOnOff_PSC_x += on_in.reshape(1, 10, 1)
        Off_SOnOff_PSC_x += off_in.reshape(1, 10, 1)
        
        # Simple voltage dynamics for SOnOff
        I_syn_SOnOff = (0.085 * On_SOnOff_PSC_s * (SOnOff_V - On_SOnOff_ESYN) + 
                        Off_SOnOff_gSYN * Off_SOnOff_PSC_s * (SOnOff_V - Off_SOnOff_ESYN))
        dV_SOnOff = ((SOnOff_E_L - SOnOff_V) - SOnOff_R * I_syn_SOnOff) / SOnOff_tau
        SOnOff_V += dt * dV_SOnOff
        
        # SOnOff spikes
        SOnOff_mask = SOnOff_V >= SOnOff_V_thresh
        SOnOff_V = np.where(SOnOff_mask, -52.0, SOnOff_V)
        SOnOff_last_spike = np.where(SOnOff_mask, t, SOnOff_last_spike)
        SOnOff_ROn_PSC_x += SOnOff_mask.astype(np.float64)
        
        # Simple voltage dynamics for ROn
        I_syn_ROn = (0.02 * On_ROn_PSC_s * (ROn_V - On_ROn_ESYN) + 
                     0.04 * Off_ROn_PSC_s * (ROn_V - Off_ROn_ESYN) +
                     SOnOff_ROn_gSYN * SOnOff_ROn_PSC_s * (ROn_V - SOnOff_ROn_ESYN))
        dV_ROn = ((ROn_E_L - ROn_V) - ROn_R * I_syn_ROn) / ROn_tau
        ROn_V += dt * dV_ROn
        
        # ROn spikes
        ROn_mask = ROn_V >= ROn_V_thresh
        ROn_V = np.where(ROn_mask, -54.0, ROn_V)
        ROn_last_spike = np.where(ROn_mask, t, ROn_last_spike)
        
        # === E-PROP COMPUTATIONS ===
        
        # Surrogate gradients
        psi_ROn = 1 - np.tanh(ROn_V - ROn_V_thresh)**2
        psi_SOnOff = 1 - np.tanh(SOnOff_V - SOnOff_V_thresh)**2
        
        # Eligibility traces
        eligibility_On_ROn = dt * psi_ROn * (-ROn_R * On_ROn_PSC_s * (ROn_V - On_ROn_ESYN)) / ROn_tau
        eligibility_Off_SOnOff = dt * psi_SOnOff * (-SOnOff_R * Off_SOnOff_PSC_s * (SOnOff_V - Off_SOnOff_ESYN)) / SOnOff_tau
        eligibility_SOnOff_ROn = dt * psi_ROn * (-ROn_R * SOnOff_ROn_PSC_s * (ROn_V - SOnOff_ROn_ESYN)) / ROn_tau
        
        # Learning signal
        target = data[:, t_idx, 0] if t_idx < data.shape[1] else np.zeros(10)
        L_t = ROn_mask.astype(np.float64).mean(axis=(0, 2)) - target  # average over batch
        
        # For hidden layer
        L_t_S = SOnOff_ROn_gSYN * L_t
        
        # L_t × eligibility
        Lt_x_elig_On_ROn = (L_t.reshape(1, 10, 1) * eligibility_On_ROn).mean()
        Lt_x_elig_Off_SOnOff = (L_t_S.reshape(1, 10, 1) * eligibility_Off_SOnOff).mean()
        Lt_x_elig_SOnOff_ROn = (L_t.reshape(1, 10, 1) * eligibility_SOnOff_ROn).mean()
        
        # Accumulate gradients
        grad_Off_SOnOff += Lt_x_elig_Off_SOnOff
        grad_On_ROn += Lt_x_elig_On_ROn
        
        # Store for plotting (average across batch and trials)
        L_t_track[t_idx] = L_t.mean()
        eligibility_Off_SOnOff_track[t_idx] = eligibility_Off_SOnOff.mean()
        eligibility_On_ROn_track[t_idx] = eligibility_On_ROn.mean()
        eligibility_SOnOff_ROn_track[t_idx] = eligibility_SOnOff_ROn.mean()
        Lt_x_elig_Off_SOnOff_track[t_idx] = Lt_x_elig_Off_SOnOff
        Lt_x_elig_On_ROn_track[t_idx] = Lt_x_elig_On_ROn
        Lt_x_elig_SOnOff_ROn_track[t_idx] = Lt_x_elig_SOnOff_ROn
        ROn_V_track[t_idx] = ROn_V.mean()
        SOnOff_V_track[t_idx] = SOnOff_V.mean()
        psi_ROn_track[t_idx] = psi_ROn.mean()
        psi_SOnOff_track[t_idx] = psi_SOnOff.mean()
        ROn_spikes_track[t_idx] = ROn_mask.mean()
        SOnOff_spikes_track[t_idx] = SOnOff_mask.mean()
        target_track[t_idx] = target.mean()
        cumulative_grad_Off_SOnOff[t_idx] = grad_Off_SOnOff
        cumulative_grad_On_ROn[t_idx] = grad_On_ROn
    
    print(f"Final accumulated gradient (Off→SOnOff): {grad_Off_SOnOff:.6f}")
    print(f"Final accumulated gradient (On→ROn): {grad_On_ROn:.6f}")
    
    return {
        'L_t': L_t_track,
        'eligibility_Off_SOnOff': eligibility_Off_SOnOff_track,
        'eligibility_On_ROn': eligibility_On_ROn_track,
        'eligibility_SOnOff_ROn': eligibility_SOnOff_ROn_track,
        'Lt_x_elig_Off_SOnOff': Lt_x_elig_Off_SOnOff_track,
        'Lt_x_elig_On_ROn': Lt_x_elig_On_ROn_track,
        'Lt_x_elig_SOnOff_ROn': Lt_x_elig_SOnOff_ROn_track,
        'ROn_V': ROn_V_track,
        'SOnOff_V': SOnOff_V_track,
        'psi_ROn': psi_ROn_track,
        'psi_SOnOff': psi_SOnOff_track,
        'ROn_spikes': ROn_spikes_track,
        'SOnOff_spikes': SOnOff_spikes_track,
        'target': target_track,
        'cumulative_grad_Off_SOnOff': cumulative_grad_Off_SOnOff,
        'cumulative_grad_On_ROn': cumulative_grad_On_ROn,
    }


def plot_results(results):
    """Plot the tracked values."""
    
    time_ms = np.arange(len(results['L_t'])) * 0.1
    
    fig, axes = plt.subplots(6, 2, figsize=(16, 18))
    fig.suptitle('E-prop Gradient Signal Analysis', fontsize=14)
    
    # Row 1: L_t and target vs output
    ax = axes[0, 0]
    ax.plot(time_ms, results['target'], 'g-', alpha=0.7, label='Target')
    ax.plot(time_ms, results['ROn_spikes'], 'b-', alpha=0.7, label='ROn spikes (rate)')
    ax.set_ylabel('Spike rate')
    ax.set_title('Target vs Output (ROn)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    ax = axes[0, 1]
    ax.plot(time_ms, results['L_t'], 'r-', alpha=0.7)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('L_t')
    ax.set_title('Learning Signal L_t = output - target')
    ax.set_xlim([0, 3000])
    
    # Row 2: Eligibility traces
    ax = axes[1, 0]
    ax.plot(time_ms, results['eligibility_On_ROn'], 'b-', alpha=0.7, label='On→ROn (exc)')
    ax.plot(time_ms, results['eligibility_SOnOff_ROn'], 'r-', alpha=0.7, label='SOnOff→ROn (inh)')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Eligibility')
    ax.set_title('Output Layer Eligibility Traces')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    ax = axes[1, 1]
    ax.plot(time_ms, results['eligibility_Off_SOnOff'], 'g-', alpha=0.7, label='Off→SOnOff (exc)')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Eligibility')
    ax.set_title('Hidden Layer Eligibility Trace (Off→SOnOff)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    # Row 3: L_t × eligibility
    ax = axes[2, 0]
    ax.plot(time_ms, results['Lt_x_elig_On_ROn'], 'b-', alpha=0.7, label='On→ROn')
    ax.plot(time_ms, results['Lt_x_elig_SOnOff_ROn'], 'r-', alpha=0.7, label='SOnOff→ROn')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('L_t × eligibility')
    ax.set_title('Instantaneous Gradient Signal (Output Layer)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    ax = axes[2, 1]
    ax.plot(time_ms, results['Lt_x_elig_Off_SOnOff'], 'g-', alpha=0.7, label='Off→SOnOff')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('L_t × eligibility')
    ax.set_title('Instantaneous Gradient Signal (Hidden Layer)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    # Row 4: Cumulative gradients
    ax = axes[3, 0]
    ax.plot(time_ms, results['cumulative_grad_On_ROn'], 'b-', alpha=0.7, label='On→ROn')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Cumulative gradient')
    ax.set_title('Cumulative Gradient (On→ROn)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    ax = axes[3, 1]
    ax.plot(time_ms, results['cumulative_grad_Off_SOnOff'], 'g-', alpha=0.7, label='Off→SOnOff')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Cumulative gradient')
    ax.set_title('Cumulative Gradient (Off→SOnOff)')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    # Row 5: Surrogate gradients (psi)
    ax = axes[4, 0]
    ax.plot(time_ms, results['psi_ROn'], 'b-', alpha=0.7, label='ψ_ROn')
    ax.plot(time_ms, results['psi_SOnOff'], 'r-', alpha=0.7, label='ψ_SOnOff')
    ax.set_ylabel('ψ (surrogate grad)')
    ax.set_title('Surrogate Gradients')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    ax = axes[4, 1]
    ax.plot(time_ms, results['ROn_V'], 'b-', alpha=0.7, label='ROn')
    ax.plot(time_ms, results['SOnOff_V'], 'r-', alpha=0.7, label='SOnOff')
    ax.axhline(-47, color='k', linestyle='--', alpha=0.3, label='threshold')
    ax.set_ylabel('Voltage (mV)')
    ax.set_title('Membrane Voltages')
    ax.legend()
    ax.set_xlim([0, 3000])
    
    # Row 6: Full time course of cumulative gradient
    ax = axes[5, 0]
    ax.plot(time_ms, results['cumulative_grad_Off_SOnOff'], 'g-', alpha=0.7)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Cumulative gradient')
    ax.set_title('Full Time Course: Cumulative Gradient (Off→SOnOff)')
    
    ax = axes[5, 1]
    ax.plot(time_ms, results['cumulative_grad_On_ROn'], 'b-', alpha=0.7)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Cumulative gradient')
    ax.set_title('Full Time Course: Cumulative Gradient (On→ROn)')
    
    plt.tight_layout()
    plt.savefig('eprop_gradient_analysis.png', dpi=150)
    plt.show()
    
    # Print summary
    print("\n" + "=" * 60)
    print("GRADIENT ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Final gradient Off→SOnOff: {results['cumulative_grad_Off_SOnOff'][-1]:.6f}")
    print(f"Final gradient On→ROn: {results['cumulative_grad_On_ROn'][-1]:.6f}")
    print()
    print("Sign interpretation:")
    if results['cumulative_grad_Off_SOnOff'][-1] > 0:
        print("  Off→SOnOff gradient > 0 → Adam will DECREASE this weight")
    else:
        print("  Off→SOnOff gradient < 0 → Adam will INCREASE this weight")
    if results['cumulative_grad_On_ROn'][-1] > 0:
        print("  On→ROn gradient > 0 → Adam will DECREASE this weight")
    else:
        print("  On→ROn gradient < 0 → Adam will INCREASE this weight")


if __name__ == "__main__":
    print("Loading data...")
    on_input, off_input, noise, data, batch_size = load_data()
    
    # Initialize parameter
    p = np.full((1, batch_size), 0.045, dtype=np.float64)
    
    print(f"Running with Off_SOnOff_gSYN = {p[0,0]}")
    results = run_with_tracking(on_input, off_input, noise, data, p)
    
    plot_results(results)
