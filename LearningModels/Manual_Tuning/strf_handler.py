import os
from scipy.io import savemat
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import yaml
from gen_strf import GenSTRF
import numpy as np
from pathlib import Path

#Function to create he STRFs
def call_strfs():
    args = ArgumentParser()
    args.add_argument('--target_dir', type=str, default='../resampled-stimuli/target', help='directory containing target stimuli')
    args.add_argument('--masker_dir', type=str, default='../resampled-stimuli/masker', help='directory containing masker stimuli')
    parsed_args = args.parse_args()
        
    yaml_path = '../config/config.yaml'
    config = yaml.safe_load(open(yaml_path, 'r'))
    strf_config = config['strf_config']

    lst_target_stim = [os.path.join(parsed_args.target_dir, stim_path) for stim_path in os.listdir(parsed_args.target_dir)]
    lst_masker_stim = [os.path.join(parsed_args.masker_dir, stim_path) for stim_path in os.listdir(parsed_args.masker_dir)]
        
    gen_strfs = GenSTRF(parsed_args, strf_config, lst_target_stim[0])
        
    target_dict = {}
    for count, stim_path in enumerate(lst_target_stim):
                
            fr_on, fr_off = gen_strfs.process_stimulus(stim_path, strf_config['targetlvl'], strf_config['stimGain'])

            #print('here')
            #print(np.shape(fr_on))

            stim_name = os.path.split(stim_path)[-1].split('.')[0]
            target_dict[f'fr_on{count}'] = fr_on
            target_dict[f'fr_off{count}'] = fr_off
        
    masker_dict = {}
    for count2, stim_path in enumerate(lst_masker_stim):
            fr_on, fr_off = gen_strfs.process_stimulus(stim_path, strf_config['maskerlvl'], strf_config['stimGain'])

            stim_name = os.path.split(stim_path)[-1].split('.')[0]
            masker_dict[f'fr{count2}'] = fr_off

    savemat("strfMat.mat", {"targets": target_dict, "maskers": masker_dict}, do_compression=True)


def export_strf_temporal_kernel_svg(output_svg: str = "strf_temporal_kernel.svg"):
        args = ArgumentParser()
        args.add_argument('--target_dir', type=str, default='../resampled-stimuli/target', help='directory containing target stimuli')
        args.add_argument('--masker_dir', type=str, default='../resampled-stimuli/masker', help='directory containing masker stimuli')
        parsed_args = args.parse_args()

        yaml_path = '../config/config.yaml'
        config = yaml.safe_load(open(yaml_path, 'r'))
        strf_config = config['strf_config']

        lst_target_stim = [os.path.join(parsed_args.target_dir, stim_path) for stim_path in os.listdir(parsed_args.target_dir)]
        if len(lst_target_stim) == 0:
                raise RuntimeError(f"No target stimuli found in {parsed_args.target_dir}")

        gen_strfs = GenSTRF(parsed_args, strf_config, lst_target_stim[0])

        t = np.asarray(gen_strfs.strf['t']).squeeze()
        h = np.asarray(gen_strfs.strf['H']).squeeze()

        fig, ax = plt.subplots(figsize=(3.4, 2.0), constrained_layout=True)
        ax.plot(t, h, color='black', linewidth=1.8)
        ax.axhline(0, color='black', linewidth=0.9)

        ax.set_title('STRF Temporal Kernel', fontsize=14, pad=6)
        ax.set_xlabel('Time (ms)', fontsize=12)
        ax.set_ylabel('Amplitude', fontsize=12)
        ax.set_xlim(0.0, 0.1)
        ax.set_xticks([])
        ax.set_yticks([])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)

        output_path = Path(output_svg)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output_path), format='svg', bbox_inches='tight')
        plt.close(fig)

        print(f"Saved STRF temporal kernel SVG: {output_path.resolve()}")


def export_onset_offset_rate_svg(
    output_svg: str = "strf_onset_offset_single_trial.svg",
    stimulus_index: int = 0,
    trial_index: int = 0,
):
        args = ArgumentParser()
        args.add_argument('--target_dir', type=str, default='../resampled-stimuli/target', help='directory containing target stimuli')
        args.add_argument('--masker_dir', type=str, default='../resampled-stimuli/masker', help='directory containing masker stimuli')
        parsed_args = args.parse_args()

        yaml_path = '../config/config.yaml'
        config = yaml.safe_load(open(yaml_path, 'r'))
        strf_config = config['strf_config']

        lst_target_stim = sorted([
                os.path.join(parsed_args.target_dir, stim_path)
                for stim_path in os.listdir(parsed_args.target_dir)
        ])
        if len(lst_target_stim) == 0:
                raise RuntimeError(f"No target stimuli found in {parsed_args.target_dir}")

        if stimulus_index < 0 or stimulus_index >= len(lst_target_stim):
                raise IndexError(f"stimulus_index {stimulus_index} out of range for {len(lst_target_stim)} target stimuli")

        gen_strfs = GenSTRF(parsed_args, strf_config, lst_target_stim[0])
        fr_on, fr_off = gen_strfs.process_stimulus(
                lst_target_stim[stimulus_index],
                strf_config['targetlvl'],
                strf_config['stimGain'],
        )

        def _single_trial_trace(x: np.ndarray, idx: int) -> np.ndarray:
                arr = np.asarray(x).squeeze()
                if arr.ndim == 1:
                        return arr
                if arr.ndim == 2:
                        if arr.shape[0] <= arr.shape[1]:
                                trial_idx = min(max(idx, 0), arr.shape[0] - 1)
                                return arr[trial_idx, :]
                        trial_idx = min(max(idx, 0), arr.shape[1] - 1)
                        return arr[:, trial_idx]
                return arr.reshape(-1)

        on_trace = _single_trial_trace(fr_on, trial_index)
        off_trace = _single_trial_trace(fr_off, trial_index)

        stim_sample_rate = float(gen_strfs.preprocStimParams.get('stimSampleRate', 10000.0))
        dt_ms = 1000.0 / stim_sample_rate
        t_ms = np.arange(on_trace.shape[0], dtype=np.float32) * dt_ms

        fig, ax = plt.subplots(figsize=(3.4, 2.0), constrained_layout=True)
        ax.plot(t_ms, on_trace, color='red', linewidth=1.0, label='Onset')
        ax.plot(t_ms, off_trace, color='green', linewidth=1.0  , label='Offset')

        ax.set_title('Onset/Offset Firing Rate', fontsize=14, pad=6)
        ax.set_xlabel('Time (ms)', fontsize=12)
        ax.set_ylabel('Rate (Hz)', fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.legend(frameon=False, fontsize=10, loc='upper right')

        output_path = Path(output_svg)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output_path), format='svg', bbox_inches='tight')
        plt.close(fig)

        print(f"Saved onset/offset rate SVG: {output_path.resolve()}")
