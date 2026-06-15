# E-prop vs Backprop Benchmark

This benchmark compares six method/backend combinations on cell 7:

- `backprop_cpu`: full PyTorch backprop through the whole stimulus window on CPU.
- `backprop_cuda`: full PyTorch backprop through the whole stimulus window on CUDA.
- `truncated_backprop_cpu`: PyTorch backprop with state detached every `--tbptt-steps` on CPU.
- `truncated_backprop_cuda`: PyTorch backprop with state detached every `--tbptt-steps` on CUDA.
- `eprop_cpu`: generated single-cell E-prop solver using the same `Forwards_Method.py` path as `E-prop/run_main_integrated.py`.
- `eprop_cupy`: CuPy version generated from that same single-cell E-prop solver body.

It uses `JointMDS_Space/LongRunResults.mat` by default and starts from the
stored final parameter set for cell 7.

`LongRunResults.mat` currently has 15 stored batch members. If you request a
larger batch, such as `20`, the benchmark tiles those final parameters to fill
the requested batch size.

## Run Everything

```bash
cd /projectnb/binaural/$USER/src/MouseSpatialGrid/LearningModels/BackPropTorch
python benchmark_eprop_vs_backprop.py
```

The default tests are:

- memory: batch sizes `5` and `20`, one epoch per method.
- speed: batch size `5`, ten epochs per method, with three parameters perturbed.

For the speed test, parameters 4, 6, and 8 are perturbed by `+15%`, `-15%`, and
`+10%` respectively so the optimizers have something to move away from.

Outputs go to `BackPropTorch/benchmark_results/`:

- `benchmark_report.md`
- `benchmark_results.json`
- `benchmark_summary.csv`
- `memory_summary.csv`
- `speed_summary.csv`
- `speed_fit_outputs.png`

## Useful Options

```bash
python benchmark_eprop_vs_backprop.py \
  --backprop-devices cpu cuda \
  --batch-sizes 5 20 \
  --speed-batch-size 5 \
  --speed-epochs 10 \
  --tbptt-steps 2000
```

For a quicker dry run:

```bash
python benchmark_eprop_vs_backprop.py --sim-len 5000 --speed-epochs 2
```

To run only the CUDA-side PyTorch tests plus both E-prop solvers:

```bash
python benchmark_eprop_vs_backprop.py --backprop-devices cuda
```

To skip one E-prop backend:

```bash
python benchmark_eprop_vs_backprop.py --skip-eprop-cpu
python benchmark_eprop_vs_backprop.py --skip-eprop-cupy
```

## Truncated BPTT Test Bench

For focused truncated-BPTT optimization without memory/speed reporting:

```bash
python tbptt_testbench.py --cell-id 7 --batch-size 5 --epochs 25 --device cuda
```

Useful variants:

```bash
python tbptt_testbench.py --cell-id 128 --batch-size 15 --tbptt-steps 200 --step-mode epoch --device cuda
python tbptt_testbench.py --cell-id 7 --batch-size 5 --tbptt-steps 2000 --step-mode chunk --device cuda
```

The default initialization is now random, using the same broad ranges as the old
E-prop `pinit` setup. To restart from `LongRunResults.mat` instead:

```bash
python tbptt_testbench.py --cell-id 7 --batch-size 5 --init-mode longrun --device cuda
python tbptt_testbench.py --cell-id 7 --batch-size 5 --init-mode perturbed-longrun --device cuda
```

`--step-mode epoch` accumulates gradients across truncated chunks and performs
one optimizer step per epoch. This is usually a better first test than stepping
after every quiet chunk, which can bias the model toward silence.

The script prints the full-stimulus PSTH sum-squared error every epoch and
writes CSV/MAT files to `BackPropTorch/tbptt_testbench_results/`.

By default it also saves one PNG per epoch showing the top 5 batches ranked by
full PSTH SSE. Each row contains the data raster, that batch's simulated raster,
and an overlaid data/sim PSTH comparison. The images are written to a sibling
folder named like `tbptt_cell_007_batch_30_chunk_random_epoch_plots/`.

Useful plotting controls:

```bash
python tbptt_testbench.py --plot-top-n 5 --plot-every 1
python tbptt_testbench.py --plot-top-n 0
```

The truncated-BPTT model learns both STRF parameters used by the E-prop parameter
vector:

- `p[0]`: STRF/input gain
- `p[1]`: STRF latency, bounded to `0.001`-`0.05` seconds

Because the PyTorch test bench starts from precomputed E-prop input spikes, these
are implemented as differentiable relative transforms around the starting input:
gain rescales the input and latency applies a fractional temporal shift.

## Notes

PyTorch CUDA memory comes from `torch.cuda.max_memory_allocated` and
`torch.cuda.max_memory_reserved`. CPU memory uses process RSS delta when
`psutil` is available.

`eprop_cupy` memory is reported from the CuPy memory pool plus generated-array
size estimates. `eprop_cpu` uses the same generated single-cell code path but
keeps arrays in NumPy. E-prop does not have a separate backward graph; the
report separates forward/output holders from eligibility-gradient holders
instead.
