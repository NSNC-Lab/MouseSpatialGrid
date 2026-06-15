# E-prop vs Backprop Benchmark Report

Cell: 7
LongRunResults source: `C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\JointMDS_Space\LongRunResults.mat`
Sim length: 29801 steps
Bin steps: 200

## Memory Summary

| Method | Batch | Epochs | Wall s | RSS delta MB | Torch peak allocated MB | Torch peak reserved MB | CuPy pool reserved MB | Best loss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| backprop_cpu | 5 | 1 | 224.07 | 6102.4 |  |  |  | 5643.0000 |
| truncated_backprop_cpu | 5 | 1 | 210.39 | 651.3 |  |  |  | 376.0714 |
| backprop_cuda | 5 | 1 | 574.94 | 8682.4 | 997.6 | 1016.0 |  | 5643.0000 |
| truncated_backprop_cuda | 5 | 1 | 347.44 | 581.7 | 88.0 | 106.0 |  | 374.7857 |
| eprop_cpu | 5 | 1 | 12.83 | 4.1 |  |  |  | 4141.0000 |
| eprop_cupy | 5 | 1 | 794.23 | 40.7 |  |  | 25.4 | 4141.0000 |
| backprop_cpu | 20 | 1 | 230.72 | 6921.9 |  |  |  | 5480.0000 |
| truncated_backprop_cpu | 20 | 1 | 214.10 | -1346.7 |  |  |  | 390.3571 |
| backprop_cuda | 20 | 1 | 583.08 | 8369.0 | 1909.6 | 1940.0 |  | 5480.0000 |
| truncated_backprop_cuda | 20 | 1 | 553.54 | -3014.1 | 212.7 | 242.0 |  | 390.7857 |
| eprop_cpu | 20 | 1 | 42.15 | 14.1 |  |  |  | 3471.0000 |
| eprop_cupy | 20 | 1 | 784.90 | 96.6 |  |  | 97.6 | 3471.0000 |

## Speed Summary

| Method | Batch | Epochs | Wall s | s / epoch | Best loss | Best batch | Best epoch |
|---|---:|---:|---:|---:|---:|---:|---:|
| backprop_cpu | 5 | 10 | 2680.75 | 268.08 | 3626.0000 | 5 | 2 |
| truncated_backprop_cpu | 5 | 10 | 1586.54 | 158.65 | 330.7143 | 2 | 2 |
| backprop_cuda | 5 | 10 | 5756.39 | 575.64 | 3626.0000 | 5 | 2 |
| truncated_backprop_cuda | 5 | 10 | 4018.49 | 401.85 | 325.5714 | 2 | 4 |
| eprop_cpu | 5 | 10 | 334.20 | 33.42 | 4179.0000 | 5 | 7 |
| eprop_cupy | 5 | 10 | 7876.87 | 787.69 | 4179.0000 | 5 | 7 |

## Memory Culprits

These culprits are broad, non-overlapping categories. PyTorch peak memory is measured through `torch.cuda.max_memory_allocated`; CuPy/E-prop CUDA memory is pool-measured plus generated-array estimates. CPU memory uses process RSS delta when `psutil` is available.

### memory: backprop_cpu, batch 5
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### memory: truncated_backprop_cpu, batch 5
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### memory: backprop_cuda, batch 5
- autograd saved timestep graph and intermediates: 974.81 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### memory: truncated_backprop_cuda, batch 5
- autograd saved timestep graph and intermediates: 65.24 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### memory: eprop_cpu, batch 5
- output/PSC/loss spike holders: 5.68 MB (estimated generated NumPy/CPU int8 holders)
- input, noise, and target arrays: 18.19 MB (estimated NumPy/CPU float32 inputs/data)
- eligibility and gradient accumulators: 0.19 MB (estimated generated NumPy/CPU float64 accumulators)
- neuron and synapse rolling state: 0.03 MB (estimated generated NumPy/CPU float64 state)
- STRF rate and derivative arrays: 2.27 MB (estimated float32 rate arrays)

### memory: eprop_cupy, batch 5
- output/PSC/loss spike holders: 5.68 MB (estimated generated CuPy/GPU int8 holders)
- input, noise, and target arrays: 18.19 MB (estimated CuPy/GPU float32 inputs/data)
- eligibility and gradient accumulators: 0.19 MB (estimated generated CuPy/GPU float64 accumulators)
- neuron and synapse rolling state: 0.03 MB (estimated generated CuPy/GPU float64 state)
- STRF rate and derivative arrays: 2.27 MB (estimated float32 rate arrays)

### memory: backprop_cpu, batch 20
- input/target/noise tensors: 90.95 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.03 MB (estimated from bin count)
- best-output trace snapshot: 22.74 MB (estimated no_grad trace copy)

### memory: truncated_backprop_cpu, batch 20
- input/target/noise tensors: 90.95 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.03 MB (estimated from bin count)
- best-output trace snapshot: 22.74 MB (estimated no_grad trace copy)

### memory: backprop_cuda, batch 20
- autograd saved timestep graph and intermediates: 1818.61 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 90.95 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.03 MB (estimated from bin count)
- best-output trace snapshot: 22.74 MB (estimated no_grad trace copy)

### memory: truncated_backprop_cuda, batch 20
- autograd saved timestep graph and intermediates: 121.70 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 90.95 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.03 MB (estimated from bin count)
- best-output trace snapshot: 22.74 MB (estimated no_grad trace copy)

### memory: eprop_cpu, batch 20
- output/PSC/loss spike holders: 22.74 MB (estimated generated NumPy/CPU int8 holders)
- input, noise, and target arrays: 69.35 MB (estimated NumPy/CPU float32 inputs/data)
- eligibility and gradient accumulators: 0.78 MB (estimated generated NumPy/CPU float64 accumulators)
- neuron and synapse rolling state: 0.14 MB (estimated generated NumPy/CPU float64 state)
- STRF rate and derivative arrays: 9.09 MB (estimated float32 rate arrays)

### memory: eprop_cupy, batch 20
- output/PSC/loss spike holders: 22.74 MB (estimated generated CuPy/GPU int8 holders)
- input, noise, and target arrays: 69.35 MB (estimated CuPy/GPU float32 inputs/data)
- eligibility and gradient accumulators: 0.78 MB (estimated generated CuPy/GPU float64 accumulators)
- neuron and synapse rolling state: 0.14 MB (estimated generated CuPy/GPU float64 state)
- STRF rate and derivative arrays: 9.09 MB (estimated float32 rate arrays)

### speed: backprop_cpu, batch 5
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### speed: truncated_backprop_cpu, batch 5
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### speed: backprop_cuda, batch 5
- autograd saved timestep graph and intermediates: 974.82 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### speed: truncated_backprop_cuda, batch 5
- autograd saved timestep graph and intermediates: 65.24 MB (measured residual from PyTorch peak allocation)
- input/target/noise tensors: 22.74 MB (measured tensor sizes)
- optimizer state, parameters, and gradients: 0.00 MB (measured/estimated tensor sizes)
- binned PSTH/loss tensors: 0.01 MB (estimated from bin count)
- best-output trace snapshot: 5.68 MB (estimated no_grad trace copy)

### speed: eprop_cpu, batch 5
- output/PSC/loss spike holders: 5.68 MB (estimated generated NumPy/CPU int8 holders)
- input, noise, and target arrays: 18.19 MB (estimated NumPy/CPU float32 inputs/data)
- eligibility and gradient accumulators: 0.19 MB (estimated generated NumPy/CPU float64 accumulators)
- neuron and synapse rolling state: 0.03 MB (estimated generated NumPy/CPU float64 state)
- STRF rate and derivative arrays: 2.27 MB (estimated float32 rate arrays)

### speed: eprop_cupy, batch 5
- output/PSC/loss spike holders: 5.68 MB (estimated generated CuPy/GPU int8 holders)
- input, noise, and target arrays: 18.19 MB (estimated CuPy/GPU float32 inputs/data)
- eligibility and gradient accumulators: 0.19 MB (estimated generated CuPy/GPU float64 accumulators)
- neuron and synapse rolling state: 0.03 MB (estimated generated CuPy/GPU float64 state)
- STRF rate and derivative arrays: 2.27 MB (estimated float32 rate arrays)

## Interpretation Checks

- Full backprop stores a timestep-by-timestep autograd graph until `loss.backward()` finishes.
- Truncated backprop should lower peak graph memory roughly in proportion to the chunk length, but it changes the gradient because state is detached between chunks.
- E-prop does not store a full autograd graph. Its persistent memory should be dominated by generated holders, eligibility accumulators, and input/output arrays.
- `eprop_cpu` is generated from `E-prop/BuildFile/Forwards_Method.py`, the same single-cell path used by `run_main_integrated.py`. `eprop_cupy` is generated by converting that same single-cell body to CuPy.
- The speed plot `speed_fit_outputs.png` compares data PSTH with starting and best outputs.