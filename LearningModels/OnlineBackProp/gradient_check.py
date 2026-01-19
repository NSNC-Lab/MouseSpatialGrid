"""
Numerical gradient check for e-prop implementation.
Compares analytical gradient with finite difference approximation.

This script tests whether the analytical gradient computed by e-prop
matches the numerical gradient computed via finite differences:
    numerical_grad = (Loss(p + ε) - Loss(p - ε)) / (2ε)
"""
import numpy as np
import sys
import os

# Change to the E-prop directory to ensure imports work
os.chdir(r'c:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\E-prop')
sys.path.insert(0, r'c:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\E-prop')

from BuildFile.generated_solve_file import solve_run


def compute_loss(spikes, target):
    """
    Compute MSE loss between output spike raster and target.
    
    Args:
        spikes: (100, 10, 1, timesteps) - batch × trials × 1 × time
        target: (10, timesteps, 1) - trials × time × 1
    
    Returns:
        Scalar MSE loss, averaged over all dimensions
    """
    # Average over batch dimension to get expected spike rate
    output = spikes.mean(axis=0)  # (10, 1, timesteps)
    output = output[:, 0, :]  # (10, timesteps)
    
    target_reshaped = target[:, :, 0]  # (10, timesteps)
    
    # Timestep-wise squared error, averaged
    loss = np.mean((output - target_reshaped)**2)
    return loss


def numerical_gradient_check(on_input, off_input, noise_token, data, base_param=0.045, epsilon=1e-4, batch_size=100):
    """
    Compute numerical gradient via central finite differences and compare to analytical e-prop gradient.
    
    The numerical gradient is:
        dL/dp ≈ (L(p + ε) - L(p - ε)) / (2ε)
    
    This is the ground truth we compare the analytical gradient against.
    """
    
    # Create parameter arrays - p is (num_params, batch_size) = (1, 100)
    # This matches InitParams.py format
    p_plus = np.full((1, batch_size), base_param + epsilon, dtype=np.float64)
    p_minus = np.full((1, batch_size), base_param - epsilon, dtype=np.float64)
    p_base = np.full((1, batch_size), base_param, dtype=np.float64)
    
    print("=" * 70)
    print("NUMERICAL GRADIENT CHECK FOR E-PROP")
    print("=" * 70)
    print(f"Testing parameter: Off_SOnOff_gSYN (conductance from Off→SOnOff)")
    print(f"Base value: {base_param}")
    print(f"Epsilon (perturbation): {epsilon}")
    print()
    
    # Run with p + epsilon
    print("Step 1/3: Running simulation with p + ε...")
    spikes_plus, grads_plus, _, _, _ = solve_run(on_input, off_input, noise_token, data, p_plus)
    loss_plus = compute_loss(spikes_plus, data)
    print(f"  Loss(p + ε) = {loss_plus:.8f}")
    
    # Run with p - epsilon
    print("Step 2/3: Running simulation with p - ε...")
    spikes_minus, grads_minus, _, _, _ = solve_run(on_input, off_input, noise_token, data, p_minus)
    loss_minus = compute_loss(spikes_minus, data)
    print(f"  Loss(p - ε) = {loss_minus:.8f}")
    
    # Run with base p to get analytical gradient from e-prop
    print("Step 3/3: Running simulation with base p to get analytical gradient...")
    spikes_base, grads_base, _, _, _ = solve_run(on_input, off_input, noise_token, data, p_base)
    loss_base = compute_loss(spikes_base, data)
    print(f"  Loss(p) = {loss_base:.8f}")
    print()
    
    # Compute numerical gradient via central difference
    numerical_grad = (loss_plus - loss_minus) / (2 * epsilon)
    
    # Get analytical gradient from e-prop (average over batch and trials)
    # grads_base shape is (1, 100, 10) after np.sum(np.stack([grad_Off_SOnOff], axis=0), axis=2)
    analytical_grad = grads_base.mean()
    
    print("=" * 70)
    print("GRADIENT COMPARISON RESULTS")
    print("=" * 70)
    print(f"Numerical gradient  (finite diff):   {numerical_grad:+.10f}")
    print(f"Analytical gradient (e-prop):        {analytical_grad:+.10f}")
    print()
    
    # Sign comparison
    num_sign = "positive" if numerical_grad > 0 else ("negative" if numerical_grad < 0 else "zero")
    ana_sign = "positive" if analytical_grad > 0 else ("negative" if analytical_grad < 0 else "zero")
    
    print(f"Numerical gradient sign:  {num_sign}")
    print(f"Analytical gradient sign: {ana_sign}")
    print()
    
    if np.sign(numerical_grad) == np.sign(analytical_grad):
        print("✓ SIGNS MATCH - Gradient direction is correct!")
    else:
        print("✗ SIGNS DO NOT MATCH - Gradient is INVERTED!")
        print("  This explains why gradient descent moves in the wrong direction.")
    print()
    
    # Magnitude comparison
    if analytical_grad != 0:
        ratio = numerical_grad / analytical_grad
        print(f"Ratio (numerical/analytical): {ratio:.6f}")
        if abs(ratio - 1) < 0.1:
            print("  ✓ Magnitudes are close (ratio ≈ 1)")
        elif abs(ratio + 1) < 0.1:
            print("  ✗ Magnitudes match but OPPOSITE SIGN (ratio ≈ -1)")
    
    if abs(numerical_grad) > 1e-10:
        rel_error = abs(numerical_grad - analytical_grad) / abs(numerical_grad)
        print(f"Relative error: {rel_error:.4f}")
    print()
    
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    if numerical_grad > 0:
        print("Numerical gradient > 0 means:")
        print("  → Increasing weight INCREASES loss")
        print("  → Gradient descent should DECREASE weight to minimize loss")
        print("  → Correct update: p_new = p_old - lr * (positive gradient)")
    elif numerical_grad < 0:
        print("Numerical gradient < 0 means:")
        print("  → Increasing weight DECREASES loss")
        print("  → Gradient descent should INCREASE weight to minimize loss")
        print("  → Correct update: p_new = p_old - lr * (negative gradient) = p_old + lr * |gradient|")
    else:
        print("Numerical gradient ≈ 0: Parameter has minimal effect on loss (near a stationary point)")
    print()
    
    if np.sign(analytical_grad) != np.sign(numerical_grad) and numerical_grad != 0:
        print("DIAGNOSIS: Your analytical gradient has the WRONG SIGN.")
        print("FIX: Either:")
        print("  1. Negate L_t in the e-prop update: L_t = target - output (instead of output - target)")
        print("  2. Or negate the gradient before passing to optimizer")
    
    return numerical_grad, analytical_grad, loss_base


def load_real_data():
    """Load the actual data from the MATLAB files used in run_main_integrated.py"""
    from scipy.io import loadmat
    from pathlib import Path
    from input_handler import call_inputs
    import set_options
    
    print("Loading real experimental data...")
    
    opts = set_options.options()
    batch_size = opts['N_batch']  # Should be 100
    
    # Load target spike data
    filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{1}contra.mat"
    data = loadmat(filename)['picture'].astype(np.float32)  # trials, timecourse
    data = data[:, :, None]  # (10, timesteps, 1)
    
    # Calculate FR for input generation
    userpath = Path(r"C:\Users\ipboy\Documents\MATLAB")
    plot_dir = userpath / "../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting"
    plot_dir = plot_dir.resolve()
    
    mat = loadmat(str(plot_dir / "all_units_info_with_polished_criteria_modified_perf.mat"),
                  variable_names=["all_data"], squeeze_me=True, struct_as_record=False)
    all_data = mat["all_data"]
    unit = all_data[7 - 1]  # n = 7
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
    
    print(f"Raw shapes from data loading:")
    print(f"  on_spks:  {on_spks.shape}")
    print(f"  off_spks: {off_spks.shape}")
    print(f"  noise:    {noise.shape}")
    print(f"  data:     {data.shape}")
    
    # The generated_solve_file expects inputs that can be indexed as on_input[:,timestep,:]
    # and broadcast with (100, 10, 1). 
    # Looking at the code: on_input[:,timestep,:] gives a 2D slice
    # on_spks is (trials=10, time, channels=1) after transpose
    # This seems like it should be (batch, time, ...) to work
    # Let me check if the file expects different input shapes
    
    print()
    
    return on_spks, off_spks, noise, data, batch_size


def generate_synthetic_data():
    """Generate synthetic test data for gradient checking."""
    print("Generating synthetic test data...")
    np.random.seed(42)
    
    timesteps = 29801
    batch = 100
    trials = 10
    
    # Simple synthetic inputs
    on_input = np.zeros((batch, timesteps, 1), dtype=np.float64)
    off_input = np.zeros((batch, timesteps, 1), dtype=np.float64)
    
    # Create sparse input patterns (Poisson-like)
    # On input: bursts at certain times
    on_input[:, 1000:2000, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    on_input[:, 5000:6000, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    on_input[:, 10000:11000, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    
    # Off input: bursts at different times  
    off_input[:, 2500:3500, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    off_input[:, 7000:8000, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    off_input[:, 12000:13000, 0] = (np.random.rand(batch, 1000) < 0.01).astype(np.float64)
    
    # Noise token
    noise_token = np.random.randn(batch, trials, timesteps, 1) * 0.01
    
    # Target data - spike when On fires (with delay for processing)
    data = np.zeros((trials, timesteps, 1), dtype=np.float32)
    data[:, 1500:2500, 0] = 0.1  # Low target rate during On-driven period
    data[:, 5500:6500, 0] = 0.1
    data[:, 10500:11500, 0] = 0.1
    
    return on_input, off_input, noise_token, data


if __name__ == "__main__":
    print()
    print("╔" + "═" * 68 + "╗")
    print("║" + " NUMERICAL GRADIENT VERIFICATION FOR E-PROP ".center(68) + "║")
    print("╚" + "═" * 68 + "╝")
    print()
    
    # Try to load real data, fall back to synthetic if it fails
    try:
        on_input, off_input, noise_token, data, batch_size = load_real_data()
        print("Using REAL experimental data\n")
    except Exception as e:
        print(f"Could not load real data: {e}")
        print("Falling back to synthetic data\n")
        on_input, off_input, noise_token, data = generate_synthetic_data()
        batch_size = 100
    
    print(f"Input shapes:")
    print(f"  on_input:    {on_input.shape}")
    print(f"  off_input:   {off_input.shape}")
    print(f"  noise_token: {noise_token.shape}")
    print(f"  data:        {data.shape}")
    print()
    
    # Run gradient check with the starting parameter value
    numerical_grad, analytical_grad, loss = numerical_gradient_check(
        on_input, off_input, noise_token, data,
        base_param=0.045,  # Starting value for Off_SOnOff_gSYN
        epsilon=1e-4,      # Perturbation size
        batch_size=batch_size
    )
    
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if np.sign(numerical_grad) == np.sign(analytical_grad):
        print("The analytical e-prop gradient matches the numerical gradient.")
        print("If training still moves in wrong direction, check the optimizer.")
    else:
        print("The analytical e-prop gradient is INVERTED relative to numerical gradient.")
        print("This is the root cause of wrong training direction.")
