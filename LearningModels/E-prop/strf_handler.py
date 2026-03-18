import os
from scipy.io import savemat
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import yaml
from gen_strf_Learnable_STRF import GenSTRF
import numpy as np

#Function to create he STRFs
def call_strfs(learnable_params,batch_size):
    base_dir = os.path.dirname(os.path.abspath(__file__))

    args = ArgumentParser()
    args.add_argument('--target_dir', type=str, default=os.path.abspath(os.path.join(base_dir, '..', 'resampled-stimuli', 'target')), help='directory containing target stimuli')
    args.add_argument('--masker_dir', type=str, default=os.path.abspath(os.path.join(base_dir, '..', 'resampled-stimuli', 'masker')), help='directory containing masker stimuli')
    parsed_args = args.parse_args()
        
    yaml_path = os.path.abspath(os.path.join(base_dir, '..', 'config', 'config.yaml'))
    config = yaml.safe_load(open(yaml_path, 'r'))
    strf_config = config['strf_config']

    lst_target_stim = [os.path.join(parsed_args.target_dir, stim_path) for stim_path in os.listdir(parsed_args.target_dir)]
    lst_masker_stim = [os.path.join(parsed_args.masker_dir, stim_path) for stim_path in os.listdir(parsed_args.masker_dir)]
        
    gen_strfs = GenSTRF(parsed_args, strf_config, lst_target_stim[0],learnable_params,batch_size)
        
    target_dict = {}
    for count, stim_path in enumerate(lst_target_stim):

        #This is super hacky, but for now I am really only looking at the first stim, so skip the second one to save time

        if count == 0:

            fr_on, fr_off, fr_on_deriv, fr_off_deriv = gen_strfs.process_stimulus(stim_path, strf_config['targetlvl'], strf_config['stimGain'])

            #print('here')
            #print(np.shape(fr_on))

            stim_name = os.path.split(stim_path)[-1].split('.')[0]
            target_dict[f'fr_on{count}'] = fr_on
            target_dict[f'fr_off{count}'] = fr_off
            target_dict[f'fr_on_deriv{count}'] = fr_on_deriv
            target_dict[f'fr_off_deriv{count}'] = fr_off_deriv
        
        #     masker_dict = {}
        #     for count2, stim_path in enumerate(lst_masker_stim):
        #             fr_on, fr_off = gen_strfs.process_stimulus(stim_path, strf_config['maskerlvl'], strf_config['stimGain'])

        #             stim_name = os.path.split(stim_path)[-1].split('.')[0]
        #             masker_dict[f'fr{count2}'] = fr_off

        #     savemat("strfMat.mat", {"targets": target_dict, "maskers": masker_dict}, do_compression=True)

    return target_dict

