# BackPropTorch single-cell prototype

This folder is a side-by-side PyTorch/backprop attempt based on the current
`LearningModels/E-prop` single-cell workflow. It does not edit the existing
E-prop scripts.

## What this does

- Loads one Oliver cell raster from `all_units_info_with_polished_criteria_modified_perf.mat`.
- Uses the same peak tuning-column logic used in recent MATLAB analysis:
  `contra -> 1`, `45 -> 2`, `center -> 3`, `ipsi -> 4`.
- Reuses the current E-prop STRF and input generation code for ON/OFF stimulus
  drive.
- Trains a small PyTorch LIF network with the same broad architecture:
  ON input, OFF input, SOnOff suppressive neuron, and ROn output neuron.
- Uses surrogate-gradient spikes so PyTorch can backpropagate through the
  spiking threshold.
- Lets you freeze or learn a scalar STRF/input gain with `--learn-strf-gain`.

The STRF learning toggle is deliberately conservative. Full STRF kernel
learning is not ported yet; `--learn-strf-gain` trains a differentiable gain on
the generated ON/OFF stimulus drive. That gives us a first memory/gradient test
without re-deriving the whole STRF convolution in PyTorch.

## Setup

Use an environment that has at least:

```bash
python -m pip install torch numpy scipy pyyaml tqdm matplotlib
```

On the SCC, this should live under the lab project space rather than your
course home, for example:

```bash
source /projectnb/binaural/$USER/venvs/eprop/bin/activate
export OLIVER_DATA_DIR=/projectnb/binaural/$USER/src/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting
cd /projectnb/binaural/$USER/src/MouseSpatialGrid/LearningModels/BackPropTorch
```

## Quick smoke test

Start with a short time window on CPU or GPU:

```bash
python run_single_cell_backprop.py --cell-id 7 --epochs 2 --max-steps 2000 --input-mode eprop-spikes
```

For a smoother differentiable input test:

```bash
python run_single_cell_backprop.py --cell-id 7 --epochs 5 --max-steps 5000 --input-mode eprop-rates --learn-strf-gain
```

## Memory-oriented runs

Full BPTT through the whole 29,801-sample stimulus is the most exact test, but
it may keep a large autograd graph in memory:

```bash
python run_single_cell_backprop.py --cell-id 7 --epochs 10 --input-mode eprop-spikes --device cuda
```

If that runs out of memory, use truncated BPTT:

```bash
python run_single_cell_backprop.py --cell-id 7 --epochs 10 --input-mode eprop-spikes --device cuda --tbptt-steps 2000
```

`--tbptt-steps` must be a multiple of `--bin-steps`. The default `--bin-steps`
is `200`, which is 20 ms at `dt = 0.1 ms`.

To measure where time and GPU memory are going, add `--profile`:

```bash
python run_single_cell_backprop.py --cell-id 7 --epochs 2 --input-mode eprop-spikes --device cuda --profile
```

The profile line reports forward, backward, best-output trace capture, optimizer
time, and PyTorch CUDA memory. `peak_allocated` is live tensor memory; `reserved`
and `peak_reserved` include PyTorch's CUDA caching allocator, which is closer to
what `nvidia-smi` shows.

## Useful arguments

- `--cell-id`: 1-based cell number.
- `--focus-mode peak`: use the tuning-type column instead of always column 1.
- `--input-mode eprop-spikes`: use sampled ON/OFF Poisson spikes from the E-prop input pipeline.
- `--input-mode eprop-rates`: use expected spike probability from STRF rates; smoother for early gradient tests.
- `--learn-strf-gain` or `--learn-strf`: make the scalar input/STRF gain trainable.
- `--max-steps`: fit only the first N samples for fast memory tests.
- `--tbptt-steps`: train in chunks and detach state between chunks.
- `--save-trace`: save the final model output raster in the `.mat` result.
- `--profile`: print timing and CUDA memory diagnostics each epoch.

Results are written to `BackPropTorch/results/torch_backprop_cell_###.mat`.

The result file always includes the best output according to the training loss:

- `best_output`: trials x time output raster for the fitted cell.
- `best_loss`: lowest saved PSTH loss.
- `best_epoch`: epoch where that best loss was observed.
- `best_batch_id`: 1-based input batch id used for that fit.
- `best_output_per_cell`, `best_loss_per_cell`, `best_epoch_per_cell`, and
  `best_batch_id_per_cell`: same values with a cell axis for easier MATLAB use.
