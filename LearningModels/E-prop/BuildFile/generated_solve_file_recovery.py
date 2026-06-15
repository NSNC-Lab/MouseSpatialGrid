# pythran export solve_run(float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:]) -> Tuple[float64[:,:,:,:], float64[:,:,:]]
import numpy as np
import cupy as cp
from BuildFile import calculate_loss_eprop
def solve_run(on_input,off_input,noise_token,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p,recovery_funcs):


    #Declare Variables


    #Transfer inputs to GPU

    on_input = cp.asarray(on_input)
    off_input = cp.asarray(off_input)
    rate_on = cp.asarray(rate_on)
    rate_off = cp.asarray(rate_off)
    rate_on_deriv = cp.asarray(rate_on_deriv)
    rate_off_deriv = cp.asarray(rate_off_deriv)
    noise_token = cp.asarray(noise_token)
    data = cp.asarray(data)
    p = cp.asarray(p)
    recovery_funcs = cp.asarray(recovery_funcs, dtype=cp.float32)
    if recovery_funcs.ndim != 2:
        raise ValueError(f"recovery_funcs must be 2D, got shape {recovery_funcs.shape}")
    if recovery_funcs.shape[0] != 220:
        raise ValueError(f"recovery_funcs first dimension must match num_cells=220, got {recovery_funcs.shape[0]}")
    print(recovery_funcs.shape[1])
    recovery_funcs = cp.clip(recovery_funcs, 0.0, 1.0)
    recovery_max_steps = recovery_funcs.shape[1]
    recovery_cell_index = cp.arange(220).reshape(220, 1, 1, 1)
    On_C = 0.1
    On_g_L = 0.005
    On_E_L = -65
    On_noise = 0
    On_t_ref = 1
    On_E_k = -80
    On_tau_ad = 10
    On_g_inc = 0.0003
    On_Itonic = 0
    On_Imask = cp.ones((1,1))
    On_R = 200.0
    On_tau = 20.0
    On_V_thresh = -47
    On_V_reset = -54
    is_output = 0
    is_noise = 0
    is_input = 1
    On_g_postIC = 0.17
    On_E_exc = 0
    On_netcon = cp.eye(1)
    On_nSYN = 0
    On_noise_E_exc = 0
    On_tauR_N = 0.7
    On_tauD_N = 1.5
    On_noise_scale = 1.9481350796278847
    final_grad_node = 0
    Off_C = 0.1
    Off_g_L = 0.005
    Off_E_L = -65
    Off_noise = 0
    Off_t_ref = 1
    Off_E_k = -80
    Off_tau_ad = 10
    Off_g_inc = 0.0003
    Off_Itonic = 0
    Off_Imask = cp.ones((1,1))
    Off_R = 200.0
    Off_tau = 20.0
    Off_V_thresh = -47
    Off_V_reset = -54
    is_output = 0
    is_noise = 0
    is_input = 1
    Off_g_postIC = 0.17
    Off_E_exc = 0
    Off_netcon = cp.eye(1)
    Off_nSYN = 0
    Off_noise_E_exc = 0
    Off_tauR_N = 0.7
    Off_tauD_N = 1.5
    Off_noise_scale = 1.9481350796278847
    final_grad_node = 0
    SOnOff_C = 0.1
    SOnOff_g_L = 0.01
    SOnOff_E_L = -57
    SOnOff_noise = 0
    SOnOff_t_ref = 0.5
    SOnOff_E_k = -80
    SOnOff_tau_ad = 100
    SOnOff_g_inc = 0.0
    SOnOff_Itonic = 0
    SOnOff_Imask = cp.ones((1,1))
    SOnOff_R = 100.0
    SOnOff_tau = 10.0
    SOnOff_V_thresh = -47
    SOnOff_V_reset = -52
    is_output = 0
    is_noise = 0
    is_input = 0
    SOnOff_g_postIC = 0.17
    SOnOff_E_exc = 0
    SOnOff_netcon = cp.eye(1)
    SOnOff_nSYN = 0
    SOnOff_noise_E_exc = 0
    SOnOff_tauR_N = 0.7
    SOnOff_tauD_N = 1.5
    SOnOff_noise_scale = 1.9481350796278847
    final_grad_node = 0
    ROn_C = 0.1
    ROn_g_L = 0.005
    ROn_E_L = -65
    ROn_noise = 0
    ROn_t_ref = 1
    ROn_E_k = -80
    ROn_tau_ad = 100
    ROn_g_inc = p[2].reshape(220, 5, 1, 1)
    ROn_Itonic = 0
    ROn_Imask = cp.ones((1,1))
    ROn_R = 200.0
    ROn_tau = 20.0
    ROn_V_thresh = -47
    ROn_V_reset = -54
    is_output = 1
    is_noise = 1
    is_input = 0
    ROn_g_postIC = 0.17
    ROn_E_exc = 0
    ROn_netcon = cp.eye(1)
    ROn_nSYN = 0.015
    ROn_noise_E_exc = 0
    ROn_tauR_N = 0.7
    ROn_tauD_N = 1.5
    ROn_noise_scale = 1.9481350796278847
    final_grad_node = 1
    On_ROn_ESYN = 0
    On_ROn_tauD = 1.5
    On_ROn_tauR = 0.7
    On_ROn_PSC_delay = 1
    On_ROn_gSYN = p[3].reshape(220, 5, 1, 1)
    On_ROn_PSC_fF = 0
    On_ROn_PSC_fP = 0.1
    On_ROn_tauF = 180
    On_ROn_tauP = 30
    On_ROn_PSC_maxF = 4
    On_ROn_netcon = cp.eye(1)
    On_ROn_scale = 1.9481350796278847
    Off_ROn_ESYN = 0
    Off_ROn_tauD = 1.5
    Off_ROn_tauR = 0.7
    Off_ROn_PSC_delay = 1
    Off_ROn_gSYN = p[4].reshape(220, 5, 1, 1)
    Off_ROn_PSC_fF = 0
    Off_ROn_PSC_fP = 0.1
    Off_ROn_tauF = 180
    Off_ROn_tauP = 30
    Off_ROn_PSC_maxF = 4
    Off_ROn_netcon = cp.eye(1)
    Off_ROn_scale = 1.9481350796278847
    On_SOnOff_ESYN = 0
    On_SOnOff_tauD = 1
    On_SOnOff_tauR = 0.1
    On_SOnOff_PSC_delay = 3
    On_SOnOff_gSYN = p[5].reshape(220, 5, 1, 1)
    On_SOnOff_PSC_fF = 0
    On_SOnOff_PSC_fP = 0.2
    On_SOnOff_tauF = 180
    On_SOnOff_tauP = 80
    On_SOnOff_PSC_maxF = 4
    On_SOnOff_netcon = cp.eye(1)
    On_SOnOff_scale = 1.2915496650148839
    Off_SOnOff_ESYN = 0
    Off_SOnOff_tauD = 1
    Off_SOnOff_tauR = 0.1
    Off_SOnOff_PSC_delay = 3
    Off_SOnOff_gSYN = p[6].reshape(220, 5, 1, 1)
    Off_SOnOff_PSC_fF = 0
    Off_SOnOff_PSC_fP = 0.0
    Off_SOnOff_tauF = 180
    Off_SOnOff_tauP = 80
    Off_SOnOff_PSC_maxF = 4
    Off_SOnOff_netcon = cp.eye(1)
    Off_SOnOff_scale = 1.2915496650148839
    SOnOff_ROn_ESYN = -80
    SOnOff_ROn_tauD = 4.5
    SOnOff_ROn_tauR = 1
    SOnOff_ROn_PSC_delay = 0.5
    SOnOff_ROn_gSYN = p[7].reshape(220, 5, 1, 1)
    SOnOff_ROn_PSC_fF = 0
    SOnOff_ROn_PSC_fP = 0.5
    SOnOff_ROn_tauF = 180
    SOnOff_ROn_tauP = 120
    SOnOff_ROn_PSC_maxF = 4
    SOnOff_ROn_netcon = cp.eye(1)
    SOnOff_ROn_scale = 1.5368523544529802
    loss_vals = cp.array([0])
    loss_bin_width = 100
    loss_bin_width2 = 300
    grad_On_ROn_accumulate=0
    grad_On_ROn_accumulate2=0
    grad_Off_ROn_accumulate=0
    grad_Off_ROn_accumulate2=0
    grad_On_SOnOff_accumulate=0
    grad_On_SOnOff_accumulate2=0
    grad_Off_SOnOff_accumulate=0
    grad_Off_SOnOff_accumulate2=0
    grad_SOnOff_ROn_accumulate=0
    grad_SOnOff_ROn_accumulate2=0
    grad_strf_gain_accumulate=0
    grad_strf_latency_accumulate=0
    grad_output_adaptation_accumulate=0
    grad_strf_gain_accumulate2=0
    grad_strf_latency_accumulate2=0
    grad_output_adaptation_accumulate2=0

    #Declare Holders

    On_V = cp.ones((220,5,10,1,2)) * On_E_L
    On_g_ad = cp.zeros((220,5,10,1,2))
    On_tspike = cp.ones((220,5,10,1,5)) * -30
    On_buffer_index = cp.ones((220,5,10,1))
    spike_wrt_tau_ad_On = cp.zeros((220,5,10,1))
    spike_wrt_g_inc_On = cp.zeros((220,5,10,1))
    Off_V = cp.ones((220,5,10,1,2)) * Off_E_L
    Off_g_ad = cp.zeros((220,5,10,1,2))
    Off_tspike = cp.ones((220,5,10,1,5)) * -30
    Off_buffer_index = cp.ones((220,5,10,1))
    spike_wrt_tau_ad_Off = cp.zeros((220,5,10,1))
    spike_wrt_g_inc_Off = cp.zeros((220,5,10,1))
    SOnOff_V = cp.ones((220,5,10,1,2)) * SOnOff_E_L
    SOnOff_g_ad = cp.zeros((220,5,10,1,2))
    SOnOff_tspike = cp.ones((220,5,10,1,5)) * -30
    SOnOff_buffer_index = cp.ones((220,5,10,1))
    spike_wrt_tau_ad_SOnOff = cp.zeros((220,5,10,1))
    spike_wrt_g_inc_SOnOff = cp.zeros((220,5,10,1))
    ROn_V = cp.ones((220,5,10,1,2)) * ROn_E_L
    ROn_g_ad = cp.zeros((220,5,10,1,2))
    ROn_tspike = cp.ones((220,5,10,1,5)) * -30
    ROn_buffer_index = cp.ones((220,5,10,1))
    ROn_spikes_holder = cp.zeros((220,5,10,1,29801), dtype=cp.int8)
    On_SOnOff_PSC_s_holder = cp.zeros((220,5,10,1,29801), dtype=cp.int8)
    Off_SOnOff_PSC_s_holder = cp.zeros((220,5,10,1,29801), dtype=cp.int8)
    losses_holder = cp.zeros((220,5,1,1,29801), dtype=cp.int8)
    ROn_noise_sn = cp.zeros((220,5,10,1,2))
    ROn_noise_xn = cp.zeros((220,5,10,1,2))
    spike_wrt_tau_ad_ROn = cp.zeros((220,5,10,1))
    spike_wrt_g_inc_ROn = cp.zeros((220,5,10,1))
    On_ROn_PSC_s = cp.zeros((220,5,10,1,2))
    On_ROn_PSC_x = cp.zeros((220,5,10,1,2))
    On_ROn_PSC_F = cp.ones((220,5,10,1,2))
    On_ROn_PSC_P = cp.ones((220,5,10,1,2))
    On_ROn_PSC_q = cp.ones((220,5,10,1,2))
    spike_wrt_gsyn_On_ROn_accumulate = cp.zeros((220,5,10,1,loss_bin_width))
    grad_On_ROn = cp.zeros((220,5,10,1))
    Off_ROn_PSC_s = cp.zeros((220,5,10,1,2))
    Off_ROn_PSC_x = cp.zeros((220,5,10,1,2))
    Off_ROn_PSC_F = cp.ones((220,5,10,1,2))
    Off_ROn_PSC_P = cp.ones((220,5,10,1,2))
    Off_ROn_PSC_q = cp.ones((220,5,10,1,2))
    spike_wrt_gsyn_Off_ROn_accumulate = cp.zeros((220,5,10,1,loss_bin_width))
    grad_Off_ROn = cp.zeros((220,5,10,1))
    On_SOnOff_PSC_s = cp.zeros((220,5,10,1,2))
    On_SOnOff_PSC_x = cp.zeros((220,5,10,1,2))
    On_SOnOff_PSC_F = cp.ones((220,5,10,1,2))
    On_SOnOff_PSC_P = cp.ones((220,5,10,1,2))
    On_SOnOff_PSC_q = cp.ones((220,5,10,1,2))
    spike_wrt_gsyn_On_SOnOff_accumulate = cp.zeros((220,5,10,1,loss_bin_width))
    grad_On_SOnOff = cp.zeros((220,5,10,1))
    Off_SOnOff_PSC_s = cp.zeros((220,5,10,1,2))
    Off_SOnOff_PSC_x = cp.zeros((220,5,10,1,2))
    Off_SOnOff_PSC_F = cp.ones((220,5,10,1,2))
    Off_SOnOff_PSC_P = cp.ones((220,5,10,1,2))
    Off_SOnOff_PSC_q = cp.ones((220,5,10,1,2))
    spike_wrt_gsyn_Off_SOnOff_accumulate = cp.zeros((220,5,10,1,loss_bin_width))
    grad_Off_SOnOff = cp.zeros((220,5,10,1))
    SOnOff_ROn_PSC_s = cp.zeros((220,5,10,1,2))
    SOnOff_ROn_PSC_x = cp.zeros((220,5,10,1,2))
    SOnOff_ROn_PSC_F = cp.ones((220,5,10,1,2))
    SOnOff_ROn_PSC_P = cp.ones((220,5,10,1,2))
    SOnOff_ROn_PSC_q = cp.ones((220,5,10,1,2))
    spike_wrt_gsyn_SOnOff_ROn_accumulate = cp.zeros((220,5,10,1,loss_bin_width))
    grad_SOnOff_ROn = cp.zeros((220,5,10,1))
    grad_strf_gain = cp.zeros((220,5,10,1))
    grad_strf_latency = cp.zeros((220,5,10,1))
    grad_output_adaptation = cp.zeros((220,5,10,1))

    for timestep,t in enumerate(cp.arange(0,29801*0.1-0.1,0.1)):


        #Declare ODES

        On_V_k1 = (((On_E_L - On_V[:,:,:,:,-1]) - On_R*On_g_ad[:,:,:,:,-1]*(On_V[:,:,:,:,-1]-On_E_k) - On_R*On_g_postIC*on_input[:,:,:,:,timestep]*On_netcon*(On_V[:,:,:,:,-1]-On_E_exc) + On_R*On_Itonic*On_Imask) / On_tau)
        On_g_ad_k1 = -On_g_ad[:,:,:,:,-1] / On_tau_ad
        Off_V_k1 = (((Off_E_L - Off_V[:,:,:,:,-1]) - Off_R*Off_g_ad[:,:,:,:,-1]*(Off_V[:,:,:,:,-1]-Off_E_k) - Off_R*Off_g_postIC*off_input[:,:,:,:,timestep]*Off_netcon*(Off_V[:,:,:,:,-1]-Off_E_exc) + Off_R*Off_Itonic*Off_Imask) / Off_tau)
        Off_g_ad_k1 = -Off_g_ad[:,:,:,:,-1] / Off_tau_ad
        SOnOff_V_k1 = (((SOnOff_E_L - SOnOff_V[:,:,:,:,-1]) - SOnOff_R*SOnOff_g_ad[:,:,:,:,-1]*(SOnOff_V[:,:,:,:,-1]-SOnOff_E_k) - SOnOff_R*(On_SOnOff_gSYN*On_SOnOff_PSC_s[:,:,:,:,-1]*On_SOnOff_netcon*(SOnOff_V[:,:,:,:,-1]-On_SOnOff_ESYN) +Off_SOnOff_gSYN*Off_SOnOff_PSC_s[:,:,:,:,-1]*Off_SOnOff_netcon*(SOnOff_V[:,:,:,:,-1]-Off_SOnOff_ESYN) ) + SOnOff_R*SOnOff_Itonic*SOnOff_Imask) / SOnOff_tau)
        SOnOff_g_ad_k1 = -SOnOff_g_ad[:,:,:,:,-1] / SOnOff_tau_ad
        ROn_V_k1 = (((ROn_E_L - ROn_V[:,:,:,:,-1]) - ROn_R*ROn_g_ad[:,:,:,:,-1]*(ROn_V[:,:,:,:,-1]-ROn_E_k) - ROn_R*(On_ROn_gSYN*On_ROn_PSC_s[:,:,:,:,-1]*On_ROn_netcon*(ROn_V[:,:,:,:,-1]-On_ROn_ESYN) +Off_ROn_gSYN*Off_ROn_PSC_s[:,:,:,:,-1]*Off_ROn_netcon*(ROn_V[:,:,:,:,-1]-Off_ROn_ESYN) +SOnOff_ROn_gSYN*SOnOff_ROn_PSC_s[:,:,:,:,-1]*SOnOff_ROn_netcon*(ROn_V[:,:,:,:,-1]-SOnOff_ROn_ESYN) ) + ROn_R*ROn_Itonic*ROn_Imask) / ROn_tau) + ((-ROn_R * ROn_nSYN * ROn_noise_sn[:,:,:,:,-1]*(ROn_V[:,:,:,:,-1]-ROn_noise_E_exc)) / ROn_tau)
        ROn_noise_sn_k1 = (ROn_noise_scale * ROn_noise_xn[:,:,:,:,-1] - ROn_noise_sn[:,:,:,:,-1]) / ROn_tauR_N
        ROn_noise_xn_k1 = -(ROn_noise_xn[:,:,:,:,-1]/ROn_tauD_N) + noise_token[:,:,:,:,timestep]/0.1
        ROn_g_ad_k1 = -ROn_g_ad[:,:,:,:,-1] / ROn_tau_ad
        On_ROn_PSC_s_k1 = (On_ROn_scale*On_ROn_PSC_x[:,:,:,:,-1] - On_ROn_PSC_s[:,:,:,:,-1]) / On_ROn_tauR
        On_ROn_PSC_x_k1 = -On_ROn_PSC_x[:,:,:,:,-1]/On_ROn_tauD
        On_ROn_PSC_F_k1 = (1 - On_ROn_PSC_F[:,:,:,:,-1])/On_ROn_tauF
        On_ROn_PSC_P_k1 = (1 - On_ROn_PSC_P[:,:,:,:,-1])/On_ROn_tauP
        On_ROn_PSC_q_k1 = 0
        Off_ROn_PSC_s_k1 = (Off_ROn_scale*Off_ROn_PSC_x[:,:,:,:,-1] - Off_ROn_PSC_s[:,:,:,:,-1]) / Off_ROn_tauR
        Off_ROn_PSC_x_k1 = -Off_ROn_PSC_x[:,:,:,:,-1]/Off_ROn_tauD
        Off_ROn_PSC_F_k1 = (1 - Off_ROn_PSC_F[:,:,:,:,-1])/Off_ROn_tauF
        Off_ROn_PSC_P_k1 = (1 - Off_ROn_PSC_P[:,:,:,:,-1])/Off_ROn_tauP
        Off_ROn_PSC_q_k1 = 0
        On_SOnOff_PSC_s_k1 = (On_SOnOff_scale*On_SOnOff_PSC_x[:,:,:,:,-1] - On_SOnOff_PSC_s[:,:,:,:,-1]) / On_SOnOff_tauR
        On_SOnOff_PSC_x_k1 = -On_SOnOff_PSC_x[:,:,:,:,-1]/On_SOnOff_tauD
        On_SOnOff_PSC_F_k1 = (1 - On_SOnOff_PSC_F[:,:,:,:,-1])/On_SOnOff_tauF
        On_SOnOff_PSC_P_k1 = (1 - On_SOnOff_PSC_P[:,:,:,:,-1])/On_SOnOff_tauP
        On_SOnOff_PSC_q_k1 = 0
        Off_SOnOff_PSC_s_k1 = (Off_SOnOff_scale*Off_SOnOff_PSC_x[:,:,:,:,-1] - Off_SOnOff_PSC_s[:,:,:,:,-1]) / Off_SOnOff_tauR
        Off_SOnOff_PSC_x_k1 = -Off_SOnOff_PSC_x[:,:,:,:,-1]/Off_SOnOff_tauD
        Off_SOnOff_PSC_F_k1 = (1 - Off_SOnOff_PSC_F[:,:,:,:,-1])/Off_SOnOff_tauF
        Off_SOnOff_PSC_P_k1 = (1 - Off_SOnOff_PSC_P[:,:,:,:,-1])/Off_SOnOff_tauP
        Off_SOnOff_PSC_q_k1 = 0
        SOnOff_ROn_PSC_s_k1 = (SOnOff_ROn_scale*SOnOff_ROn_PSC_x[:,:,:,:,-1] - SOnOff_ROn_PSC_s[:,:,:,:,-1]) / SOnOff_ROn_tauR
        SOnOff_ROn_PSC_x_k1 = -SOnOff_ROn_PSC_x[:,:,:,:,-1]/SOnOff_ROn_tauD
        SOnOff_ROn_PSC_F_k1 = (1 - SOnOff_ROn_PSC_F[:,:,:,:,-1])/SOnOff_ROn_tauF
        SOnOff_ROn_PSC_P_k1 = (1 - SOnOff_ROn_PSC_P[:,:,:,:,-1])/SOnOff_ROn_tauP
        SOnOff_ROn_PSC_q_k1 = 0

        #Declare State Updates

        On_V[:,:,:,:,-2] = On_V[:,:,:,:,-1]
        On_V[:,:,:,:,-1] = On_V[:,:,:,:,-1] + 0.1*On_V_k1
        On_g_ad[:,:,:,:,-2] = On_g_ad[:,:,:,:,-1]
        On_g_ad[:,:,:,:,-1] = On_g_ad[:,:,:,:,-1] + 0.1*On_g_ad_k1
        Off_V[:,:,:,:,-2] = Off_V[:,:,:,:,-1]
        Off_V[:,:,:,:,-1] = Off_V[:,:,:,:,-1] + 0.1*Off_V_k1
        Off_g_ad[:,:,:,:,-2] = Off_g_ad[:,:,:,:,-1]
        Off_g_ad[:,:,:,:,-1] = Off_g_ad[:,:,:,:,-1] + 0.1*Off_g_ad_k1
        SOnOff_V[:,:,:,:,-2] = SOnOff_V[:,:,:,:,-1]
        SOnOff_V[:,:,:,:,-1] = SOnOff_V[:,:,:,:,-1] + 0.1*SOnOff_V_k1
        SOnOff_g_ad[:,:,:,:,-2] = SOnOff_g_ad[:,:,:,:,-1]
        SOnOff_g_ad[:,:,:,:,-1] = SOnOff_g_ad[:,:,:,:,-1] + 0.1*SOnOff_g_ad_k1
        ROn_V[:,:,:,:,-2] = ROn_V[:,:,:,:,-1]
        ROn_V[:,:,:,:,-1] = ROn_V[:,:,:,:,-1] + 0.1*ROn_V_k1
        ROn_g_ad[:,:,:,:,-2] = ROn_g_ad[:,:,:,:,-1]
        ROn_g_ad[:,:,:,:,-1] = ROn_g_ad[:,:,:,:,-1] + 0.1*ROn_g_ad_k1
        ROn_noise_sn[:,:,:,:,-2] = ROn_noise_sn[:,:,:,:,-1]
        ROn_noise_sn[:,:,:,:,-1] = ROn_noise_sn[:,:,:,:,-1] + 0.1*ROn_noise_sn_k1
        ROn_noise_xn[:,:,:,:,-2] = ROn_noise_xn[:,:,:,:,-1]
        ROn_noise_xn[:,:,:,:,-1] = ROn_noise_xn[:,:,:,:,-1] + 0.1*ROn_noise_xn_k1
        On_ROn_PSC_s[:,:,:,:,-2] = On_ROn_PSC_s[:,:,:,:,-1]
        On_ROn_PSC_s[:,:,:,:,-1] = On_ROn_PSC_s[:,:,:,:,-1] + 0.1*On_ROn_PSC_s_k1
        On_ROn_PSC_x[:,:,:,:,-2] = On_ROn_PSC_x[:,:,:,:,-1]
        On_ROn_PSC_x[:,:,:,:,-1] = On_ROn_PSC_x[:,:,:,:,-1] + 0.1*On_ROn_PSC_x_k1
        On_ROn_PSC_F[:,:,:,:,-2] = On_ROn_PSC_F[:,:,:,:,-1]
        On_ROn_PSC_F[:,:,:,:,-1] = On_ROn_PSC_F[:,:,:,:,-1] + 0.1*On_ROn_PSC_F_k1
        On_ROn_PSC_P[:,:,:,:,-2] = On_ROn_PSC_P[:,:,:,:,-1]
        On_ROn_PSC_P[:,:,:,:,-1] = On_ROn_PSC_P[:,:,:,:,-1] + 0.1*On_ROn_PSC_P_k1
        On_ROn_PSC_q[:,:,:,:,-2] = On_ROn_PSC_q[:,:,:,:,-1]
        On_ROn_PSC_q[:,:,:,:,-1] = On_ROn_PSC_q[:,:,:,:,-1] + 0.1*On_ROn_PSC_q_k1
        Off_ROn_PSC_s[:,:,:,:,-2] = Off_ROn_PSC_s[:,:,:,:,-1]
        Off_ROn_PSC_s[:,:,:,:,-1] = Off_ROn_PSC_s[:,:,:,:,-1] + 0.1*Off_ROn_PSC_s_k1
        Off_ROn_PSC_x[:,:,:,:,-2] = Off_ROn_PSC_x[:,:,:,:,-1]
        Off_ROn_PSC_x[:,:,:,:,-1] = Off_ROn_PSC_x[:,:,:,:,-1] + 0.1*Off_ROn_PSC_x_k1
        Off_ROn_PSC_F[:,:,:,:,-2] = Off_ROn_PSC_F[:,:,:,:,-1]
        Off_ROn_PSC_F[:,:,:,:,-1] = Off_ROn_PSC_F[:,:,:,:,-1] + 0.1*Off_ROn_PSC_F_k1
        Off_ROn_PSC_P[:,:,:,:,-2] = Off_ROn_PSC_P[:,:,:,:,-1]
        Off_ROn_PSC_P[:,:,:,:,-1] = Off_ROn_PSC_P[:,:,:,:,-1] + 0.1*Off_ROn_PSC_P_k1
        Off_ROn_PSC_q[:,:,:,:,-2] = Off_ROn_PSC_q[:,:,:,:,-1]
        Off_ROn_PSC_q[:,:,:,:,-1] = Off_ROn_PSC_q[:,:,:,:,-1] + 0.1*Off_ROn_PSC_q_k1
        On_SOnOff_PSC_s[:,:,:,:,-2] = On_SOnOff_PSC_s[:,:,:,:,-1]
        On_SOnOff_PSC_s[:,:,:,:,-1] = On_SOnOff_PSC_s[:,:,:,:,-1] + 0.1*On_SOnOff_PSC_s_k1
        On_SOnOff_PSC_x[:,:,:,:,-2] = On_SOnOff_PSC_x[:,:,:,:,-1]
        On_SOnOff_PSC_x[:,:,:,:,-1] = On_SOnOff_PSC_x[:,:,:,:,-1] + 0.1*On_SOnOff_PSC_x_k1
        On_SOnOff_PSC_F[:,:,:,:,-2] = On_SOnOff_PSC_F[:,:,:,:,-1]
        On_SOnOff_PSC_F[:,:,:,:,-1] = On_SOnOff_PSC_F[:,:,:,:,-1] + 0.1*On_SOnOff_PSC_F_k1
        On_SOnOff_PSC_P[:,:,:,:,-2] = On_SOnOff_PSC_P[:,:,:,:,-1]
        On_SOnOff_PSC_P[:,:,:,:,-1] = On_SOnOff_PSC_P[:,:,:,:,-1] + 0.1*On_SOnOff_PSC_P_k1
        On_SOnOff_PSC_q[:,:,:,:,-2] = On_SOnOff_PSC_q[:,:,:,:,-1]
        On_SOnOff_PSC_q[:,:,:,:,-1] = On_SOnOff_PSC_q[:,:,:,:,-1] + 0.1*On_SOnOff_PSC_q_k1
        Off_SOnOff_PSC_s[:,:,:,:,-2] = Off_SOnOff_PSC_s[:,:,:,:,-1]
        Off_SOnOff_PSC_s[:,:,:,:,-1] = Off_SOnOff_PSC_s[:,:,:,:,-1] + 0.1*Off_SOnOff_PSC_s_k1
        Off_SOnOff_PSC_x[:,:,:,:,-2] = Off_SOnOff_PSC_x[:,:,:,:,-1]
        Off_SOnOff_PSC_x[:,:,:,:,-1] = Off_SOnOff_PSC_x[:,:,:,:,-1] + 0.1*Off_SOnOff_PSC_x_k1
        Off_SOnOff_PSC_F[:,:,:,:,-2] = Off_SOnOff_PSC_F[:,:,:,:,-1]
        Off_SOnOff_PSC_F[:,:,:,:,-1] = Off_SOnOff_PSC_F[:,:,:,:,-1] + 0.1*Off_SOnOff_PSC_F_k1
        Off_SOnOff_PSC_P[:,:,:,:,-2] = Off_SOnOff_PSC_P[:,:,:,:,-1]
        Off_SOnOff_PSC_P[:,:,:,:,-1] = Off_SOnOff_PSC_P[:,:,:,:,-1] + 0.1*Off_SOnOff_PSC_P_k1
        Off_SOnOff_PSC_q[:,:,:,:,-2] = Off_SOnOff_PSC_q[:,:,:,:,-1]
        Off_SOnOff_PSC_q[:,:,:,:,-1] = Off_SOnOff_PSC_q[:,:,:,:,-1] + 0.1*Off_SOnOff_PSC_q_k1
        SOnOff_ROn_PSC_s[:,:,:,:,-2] = SOnOff_ROn_PSC_s[:,:,:,:,-1]
        SOnOff_ROn_PSC_s[:,:,:,:,-1] = SOnOff_ROn_PSC_s[:,:,:,:,-1] + 0.1*SOnOff_ROn_PSC_s_k1
        SOnOff_ROn_PSC_x[:,:,:,:,-2] = SOnOff_ROn_PSC_x[:,:,:,:,-1]
        SOnOff_ROn_PSC_x[:,:,:,:,-1] = SOnOff_ROn_PSC_x[:,:,:,:,-1] + 0.1*SOnOff_ROn_PSC_x_k1
        SOnOff_ROn_PSC_F[:,:,:,:,-2] = SOnOff_ROn_PSC_F[:,:,:,:,-1]
        SOnOff_ROn_PSC_F[:,:,:,:,-1] = SOnOff_ROn_PSC_F[:,:,:,:,-1] + 0.1*SOnOff_ROn_PSC_F_k1
        SOnOff_ROn_PSC_P[:,:,:,:,-2] = SOnOff_ROn_PSC_P[:,:,:,:,-1]
        SOnOff_ROn_PSC_P[:,:,:,:,-1] = SOnOff_ROn_PSC_P[:,:,:,:,-1] + 0.1*SOnOff_ROn_PSC_P_k1
        SOnOff_ROn_PSC_q[:,:,:,:,-2] = SOnOff_ROn_PSC_q[:,:,:,:,-1]
        SOnOff_ROn_PSC_q[:,:,:,:,-1] = SOnOff_ROn_PSC_q[:,:,:,:,-1] + 0.1*SOnOff_ROn_PSC_q_k1
        #Compute Phi

        psi_On = (1 - cp.tanh(On_V[:,:,:,:,-1] - On_V_thresh)**2)
        psi_Off = (1 - cp.tanh(Off_V[:,:,:,:,-1] - Off_V_thresh)**2)
        psi_ROn = (1 - cp.tanh(ROn_V[:,:,:,:,-1] - ROn_V_thresh)**2)
        psi_SOnOff = (1 - cp.tanh(SOnOff_V[:,:,:,:,-1] - SOnOff_V_thresh)**2)

        eligibility_On_ROn = 0.1*(psi_ROn*(-ROn_R) * On_ROn_PSC_s[:,:,:,:,-1]*On_ROn_netcon*(ROn_V[:,:,:,:,-1]-On_ROn_ESYN))/ROn_tau
        eligibility_Off_ROn = 0.1*(psi_ROn*(-ROn_R) * Off_ROn_PSC_s[:,:,:,:,-1]*Off_ROn_netcon*(ROn_V[:,:,:,:,-1]-Off_ROn_ESYN))/ROn_tau
        eligibility_On_SOnOff = 0.1*(psi_SOnOff*(-SOnOff_R) * On_SOnOff_PSC_s[:,:,:,:,-1]*On_SOnOff_netcon*(SOnOff_V[:,:,:,:,-1]-On_SOnOff_ESYN))/SOnOff_tau
        eligibility_Off_SOnOff = 0.1*(psi_SOnOff*(-SOnOff_R) * Off_SOnOff_PSC_s[:,:,:,:,-1]*Off_SOnOff_netcon*(SOnOff_V[:,:,:,:,-1]-Off_SOnOff_ESYN))/SOnOff_tau
        eligibility_SOnOff_ROn = 0.1*(psi_ROn*(-ROn_R) * SOnOff_ROn_PSC_s[:,:,:,:,-1]*SOnOff_ROn_netcon*(ROn_V[:,:,:,:,-1]-SOnOff_ROn_ESYN))/ROn_tau
        eligibility_strf_gain = ((psi_On * (-On_R) * On_g_postIC * 30 * (1 - cp.tanh(0.0001*rate_on[timestep,:,:][:,:,None,None]-1.5)**2) * On_netcon * (On_V[:,:,:,:,-1]-On_E_exc)) / On_tau) + ((psi_Off * (-Off_R) * Off_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_off[timestep,:,:][:,:,None,None]-1.5)**2) * Off_netcon * (Off_V[:,:,:,:,-1]-Off_E_exc)) / Off_tau)
        eligibility_strf_latency = ((psi_On * (-On_R) * On_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_on_deriv[timestep,:,:][:,:,None,None]-1.5)**2) * On_netcon * (On_V[:,:,:,:,-1]-On_E_exc)) / On_tau) + ((psi_Off * (-Off_R) * Off_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_off_deriv[timestep,:,:][:,:,None,None]-1.5)**2) * Off_netcon * (Off_V[:,:,:,:,-1]-Off_E_exc)) / Off_tau)
        eligibility_output_adaptation = (-ROn_R*(((1+cp.tanh(ROn_V[:,:,:,:,-1]-ROn_V_thresh))/2)*((1+cp.tanh(-(ROn_V[:,:,:,:,-2]-ROn_V_thresh)))/2))*0.1*(ROn_V[:,:,:,:,-1]-ROn_E_k))/ROn_tau

        #Compute Bk

        Bk = -SOnOff_ROn_gSYN*(ROn_V_thresh - SOnOff_ROn_ESYN)


        #Declare Conditionals

        On_last_spike = cp.max(On_tspike, axis=-1)
        On_steps_since_spike = cp.rint((t - On_last_spike) / 0.1).astype(cp.int32)
        On_steps_since_spike = cp.maximum(On_steps_since_spike, 0)
        On_recovery_idx = cp.minimum(On_steps_since_spike, recovery_max_steps - 1)
        On_recovery_prob = recovery_funcs[recovery_cell_index, On_recovery_idx]
        On_recovery_allowed = (On_steps_since_spike >= recovery_max_steps) | (cp.random.random(On_steps_since_spike.shape) <= On_recovery_prob)
        On_mask_2b = (On_V[:,:,:,:,-1] >= On_V_thresh) & (~On_recovery_allowed)
        if cp.any(On_mask_2b).item():
            On_spikers_2b = cp.where(On_mask_2b)
            On_V[On_spikers_2b + (-2,)] = On_V[On_spikers_2b + (-1,)]
            On_V[On_spikers_2b + (-1,)] = On_V_reset
        On_mask_1 = ((On_V[:,:,:,:,-1] >= On_V_thresh) & (On_V[:,:,:,:,-2] < On_V_thresh) & On_recovery_allowed).astype(cp.int8)
        if cp.any(On_mask_1).item():
            On_spikers_1 = cp.where(On_mask_1)
            On_tspike[On_spikers_1 + (On_buffer_index[On_spikers_1].astype(cp.int8)-1,)] = t
            On_buffer_index[On_spikers_1] = On_buffer_index[On_spikers_1] % 5 + 1
        On_mask_2a = (On_V[:,:,:,:,-1] >= On_V_thresh)
        if cp.any(On_mask_2a).item():
            On_spikers_2a = cp.where(On_mask_2a)
            On_V[On_spikers_2a + (-2,)] = On_V[On_spikers_2a + (-1,)]
            On_V[On_spikers_2a + (-1,)] = On_V_reset
            On_g_ad[On_spikers_2a + (-2,)] = On_g_ad[On_spikers_2a + (-1,)]
            On_g_ad[On_spikers_2a + (-1,)] = On_g_ad[On_spikers_2a + (-1,)] + On_g_inc
        Off_last_spike = cp.max(Off_tspike, axis=-1)
        Off_steps_since_spike = cp.rint((t - Off_last_spike) / 0.1).astype(cp.int32)
        Off_steps_since_spike = cp.maximum(Off_steps_since_spike, 0)
        Off_recovery_idx = cp.minimum(Off_steps_since_spike, recovery_max_steps - 1)
        Off_recovery_prob = recovery_funcs[recovery_cell_index, Off_recovery_idx]
        Off_recovery_allowed = (Off_steps_since_spike >= recovery_max_steps) | (cp.random.random(Off_steps_since_spike.shape) <= Off_recovery_prob)
        Off_mask_2b = (Off_V[:,:,:,:,-1] >= Off_V_thresh) & (~Off_recovery_allowed)
        if cp.any(Off_mask_2b).item():
            Off_spikers_2b = cp.where(Off_mask_2b)
            Off_V[Off_spikers_2b + (-2,)] = Off_V[Off_spikers_2b + (-1,)]
            Off_V[Off_spikers_2b + (-1,)] = Off_V_reset
        Off_mask_1 = ((Off_V[:,:,:,:,-1] >= Off_V_thresh) & (Off_V[:,:,:,:,-2] < Off_V_thresh) & Off_recovery_allowed).astype(cp.int8)
        if cp.any(Off_mask_1).item():
            Off_spikers_1 = cp.where(Off_mask_1)
            Off_tspike[Off_spikers_1 + (Off_buffer_index[Off_spikers_1].astype(cp.int8)-1,)] = t
            Off_buffer_index[Off_spikers_1] = Off_buffer_index[Off_spikers_1] % 5 + 1
        Off_mask_2a = (Off_V[:,:,:,:,-1] >= Off_V_thresh)
        if cp.any(Off_mask_2a).item():
            Off_spikers_2a = cp.where(Off_mask_2a)
            Off_V[Off_spikers_2a + (-2,)] = Off_V[Off_spikers_2a + (-1,)]
            Off_V[Off_spikers_2a + (-1,)] = Off_V_reset
            Off_g_ad[Off_spikers_2a + (-2,)] = Off_g_ad[Off_spikers_2a + (-1,)]
            Off_g_ad[Off_spikers_2a + (-1,)] = Off_g_ad[Off_spikers_2a + (-1,)] + Off_g_inc
        SOnOff_last_spike = cp.max(SOnOff_tspike, axis=-1)
        SOnOff_steps_since_spike = cp.rint((t - SOnOff_last_spike) / 0.1).astype(cp.int32)
        SOnOff_steps_since_spike = cp.maximum(SOnOff_steps_since_spike, 0)
        SOnOff_recovery_idx = cp.minimum(SOnOff_steps_since_spike, recovery_max_steps - 1)
        SOnOff_recovery_prob = recovery_funcs[recovery_cell_index, SOnOff_recovery_idx]
        SOnOff_recovery_allowed = (SOnOff_steps_since_spike >= recovery_max_steps) | (cp.random.random(SOnOff_steps_since_spike.shape) <= SOnOff_recovery_prob)
        SOnOff_mask_2b = (SOnOff_V[:,:,:,:,-1] >= SOnOff_V_thresh) & (~SOnOff_recovery_allowed)
        if cp.any(SOnOff_mask_2b).item():
            SOnOff_spikers_2b = cp.where(SOnOff_mask_2b)
            SOnOff_V[SOnOff_spikers_2b + (-2,)] = SOnOff_V[SOnOff_spikers_2b + (-1,)]
            SOnOff_V[SOnOff_spikers_2b + (-1,)] = SOnOff_V_reset
        SOnOff_mask_1 = ((SOnOff_V[:,:,:,:,-1] >= SOnOff_V_thresh) & (SOnOff_V[:,:,:,:,-2] < SOnOff_V_thresh) & SOnOff_recovery_allowed).astype(cp.int8)
        if cp.any(SOnOff_mask_1).item():
            SOnOff_spikers_1 = cp.where(SOnOff_mask_1)
            SOnOff_tspike[SOnOff_spikers_1 + (SOnOff_buffer_index[SOnOff_spikers_1].astype(cp.int8)-1,)] = t
            SOnOff_buffer_index[SOnOff_spikers_1] = SOnOff_buffer_index[SOnOff_spikers_1] % 5 + 1
        SOnOff_mask_2a = (SOnOff_V[:,:,:,:,-1] >= SOnOff_V_thresh)
        if cp.any(SOnOff_mask_2a).item():
            SOnOff_spikers_2a = cp.where(SOnOff_mask_2a)
            SOnOff_V[SOnOff_spikers_2a + (-2,)] = SOnOff_V[SOnOff_spikers_2a + (-1,)]
            SOnOff_V[SOnOff_spikers_2a + (-1,)] = SOnOff_V_reset
            SOnOff_g_ad[SOnOff_spikers_2a + (-2,)] = SOnOff_g_ad[SOnOff_spikers_2a + (-1,)]
            SOnOff_g_ad[SOnOff_spikers_2a + (-1,)] = SOnOff_g_ad[SOnOff_spikers_2a + (-1,)] + SOnOff_g_inc
        ROn_last_spike = cp.max(ROn_tspike, axis=-1)
        ROn_steps_since_spike = cp.rint((t - ROn_last_spike) / 0.1).astype(cp.int32)
        ROn_steps_since_spike = cp.maximum(ROn_steps_since_spike, 0)
        ROn_recovery_idx = cp.minimum(ROn_steps_since_spike, recovery_max_steps - 1)
        ROn_recovery_prob = recovery_funcs[recovery_cell_index, ROn_recovery_idx]
        ROn_recovery_allowed = (ROn_steps_since_spike >= recovery_max_steps) | (cp.random.random(ROn_steps_since_spike.shape) <= ROn_recovery_prob)
        ROn_mask_2b = (ROn_V[:,:,:,:,-1] >= ROn_V_thresh) & (~ROn_recovery_allowed)
        if cp.any(ROn_mask_2b).item():
            ROn_spikers_2b = cp.where(ROn_mask_2b)
            ROn_V[ROn_spikers_2b + (-2,)] = ROn_V[ROn_spikers_2b + (-1,)]
            ROn_V[ROn_spikers_2b + (-1,)] = ROn_V_reset
        ROn_mask_1 = ((ROn_V[:,:,:,:,-1] >= ROn_V_thresh) & (ROn_V[:,:,:,:,-2] < ROn_V_thresh) & ROn_recovery_allowed).astype(cp.int8)
        ROn_spikes_holder[:,:,:,:,timestep] = ROn_mask_1
        if cp.any(ROn_mask_1).item():
            ROn_spikers_1 = cp.where(ROn_mask_1)
            ROn_tspike[ROn_spikers_1 + (ROn_buffer_index[ROn_spikers_1].astype(cp.int8)-1,)] = t
            ROn_buffer_index[ROn_spikers_1] = ROn_buffer_index[ROn_spikers_1] % 5 + 1
        ROn_mask_2a = (ROn_V[:,:,:,:,-1] >= ROn_V_thresh)
        if cp.any(ROn_mask_2a).item():
            ROn_spikers_2a = cp.where(ROn_mask_2a)
            ROn_V[ROn_spikers_2a + (-2,)] = ROn_V[ROn_spikers_2a + (-1,)]
            ROn_V[ROn_spikers_2a + (-1,)] = ROn_V_reset
            ROn_g_ad[ROn_spikers_2a + (-2,)] = ROn_g_ad[ROn_spikers_2a + (-1,)]
            ROn_g_ad[ROn_spikers_2a + (-1,)] = ROn_g_ad[ROn_spikers_2a + (-1,)] + ROn_g_inc[ROn_spikers_2a[0],ROn_spikers_2a[1],0,0]
        On_ROn_mask_3 = cp.any((t == (On_tspike + On_ROn_PSC_delay)), axis=-1)
        if cp.any(On_ROn_mask_3).item():
            On_ROn_spikers_3 = cp.where(On_ROn_mask_3)
            On_ROn_PSC_x[On_ROn_spikers_3 + (-2,)] = On_ROn_PSC_x[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_q[On_ROn_spikers_3 + (-2,)] = On_ROn_PSC_q[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_F[On_ROn_spikers_3 + (-2,)] = On_ROn_PSC_F[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_P[On_ROn_spikers_3 + (-2,)] = On_ROn_PSC_P[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_x[On_ROn_spikers_3 + (-1,)] = On_ROn_PSC_x[On_ROn_spikers_3 + (-1,)] + On_ROn_PSC_q[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_q[On_ROn_spikers_3 + (-1,)] = On_ROn_PSC_F[On_ROn_spikers_3 + (-1,)] * On_ROn_PSC_P[On_ROn_spikers_3 + (-1,)]
            On_ROn_PSC_F[On_ROn_spikers_3 + (-1,)] = On_ROn_PSC_F[On_ROn_spikers_3 + (-1,)] + On_ROn_PSC_fF *(On_ROn_PSC_maxF - On_ROn_PSC_F[On_ROn_spikers_3 + (-1,)])
            On_ROn_PSC_P[On_ROn_spikers_3 + (-1,)] = On_ROn_PSC_P[On_ROn_spikers_3 + (-1,)] * (1 - On_ROn_PSC_fP)
        Off_ROn_mask_3 = cp.any((t == (Off_tspike + Off_ROn_PSC_delay)), axis=-1)
        if cp.any(Off_ROn_mask_3).item():
            Off_ROn_spikers_3 = cp.where(Off_ROn_mask_3)
            Off_ROn_PSC_x[Off_ROn_spikers_3 + (-2,)] = Off_ROn_PSC_x[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_q[Off_ROn_spikers_3 + (-2,)] = Off_ROn_PSC_q[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_F[Off_ROn_spikers_3 + (-2,)] = Off_ROn_PSC_F[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_P[Off_ROn_spikers_3 + (-2,)] = Off_ROn_PSC_P[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_x[Off_ROn_spikers_3 + (-1,)] = Off_ROn_PSC_x[Off_ROn_spikers_3 + (-1,)] + Off_ROn_PSC_q[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_q[Off_ROn_spikers_3 + (-1,)] = Off_ROn_PSC_F[Off_ROn_spikers_3 + (-1,)] * Off_ROn_PSC_P[Off_ROn_spikers_3 + (-1,)]
            Off_ROn_PSC_F[Off_ROn_spikers_3 + (-1,)] = Off_ROn_PSC_F[Off_ROn_spikers_3 + (-1,)] + Off_ROn_PSC_fF *(Off_ROn_PSC_maxF - Off_ROn_PSC_F[Off_ROn_spikers_3 + (-1,)])
            Off_ROn_PSC_P[Off_ROn_spikers_3 + (-1,)] = Off_ROn_PSC_P[Off_ROn_spikers_3 + (-1,)] * (1 - Off_ROn_PSC_fP)
        On_SOnOff_mask_3 = cp.any((t == (On_tspike + On_SOnOff_PSC_delay)), axis=-1)
        if cp.any(On_SOnOff_mask_3).item():
            On_SOnOff_spikers_3 = cp.where(On_SOnOff_mask_3)
            On_SOnOff_PSC_x[On_SOnOff_spikers_3 + (-2,)] = On_SOnOff_PSC_x[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_q[On_SOnOff_spikers_3 + (-2,)] = On_SOnOff_PSC_q[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-2,)] = On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_P[On_SOnOff_spikers_3 + (-2,)] = On_SOnOff_PSC_P[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_x[On_SOnOff_spikers_3 + (-1,)] = On_SOnOff_PSC_x[On_SOnOff_spikers_3 + (-1,)] + On_SOnOff_PSC_q[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_q[On_SOnOff_spikers_3 + (-1,)] = On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-1,)] * On_SOnOff_PSC_P[On_SOnOff_spikers_3 + (-1,)]
            On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-1,)] = On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-1,)] + On_SOnOff_PSC_fF *(On_SOnOff_PSC_maxF - On_SOnOff_PSC_F[On_SOnOff_spikers_3 + (-1,)])
            On_SOnOff_PSC_P[On_SOnOff_spikers_3 + (-1,)] = On_SOnOff_PSC_P[On_SOnOff_spikers_3 + (-1,)] * (1 - On_SOnOff_PSC_fP)
        Off_SOnOff_mask_3 = cp.any((t == (Off_tspike + Off_SOnOff_PSC_delay)), axis=-1)
        if cp.any(Off_SOnOff_mask_3).item():
            Off_SOnOff_spikers_3 = cp.where(Off_SOnOff_mask_3)
            Off_SOnOff_PSC_x[Off_SOnOff_spikers_3 + (-2,)] = Off_SOnOff_PSC_x[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_q[Off_SOnOff_spikers_3 + (-2,)] = Off_SOnOff_PSC_q[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-2,)] = Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_P[Off_SOnOff_spikers_3 + (-2,)] = Off_SOnOff_PSC_P[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_x[Off_SOnOff_spikers_3 + (-1,)] = Off_SOnOff_PSC_x[Off_SOnOff_spikers_3 + (-1,)] + Off_SOnOff_PSC_q[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_q[Off_SOnOff_spikers_3 + (-1,)] = Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-1,)] * Off_SOnOff_PSC_P[Off_SOnOff_spikers_3 + (-1,)]
            Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-1,)] = Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-1,)] + Off_SOnOff_PSC_fF *(Off_SOnOff_PSC_maxF - Off_SOnOff_PSC_F[Off_SOnOff_spikers_3 + (-1,)])
            Off_SOnOff_PSC_P[Off_SOnOff_spikers_3 + (-1,)] = Off_SOnOff_PSC_P[Off_SOnOff_spikers_3 + (-1,)] * (1 - Off_SOnOff_PSC_fP)
        SOnOff_ROn_mask_3 = cp.any((t == (SOnOff_tspike + SOnOff_ROn_PSC_delay)), axis=-1)
        if cp.any(SOnOff_ROn_mask_3).item():
            SOnOff_ROn_spikers_3 = cp.where(SOnOff_ROn_mask_3)
            SOnOff_ROn_PSC_x[SOnOff_ROn_spikers_3 + (-2,)] = SOnOff_ROn_PSC_x[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_q[SOnOff_ROn_spikers_3 + (-2,)] = SOnOff_ROn_PSC_q[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-2,)] = SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_P[SOnOff_ROn_spikers_3 + (-2,)] = SOnOff_ROn_PSC_P[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_x[SOnOff_ROn_spikers_3 + (-1,)] = SOnOff_ROn_PSC_x[SOnOff_ROn_spikers_3 + (-1,)] + SOnOff_ROn_PSC_q[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_q[SOnOff_ROn_spikers_3 + (-1,)] = SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-1,)] * SOnOff_ROn_PSC_P[SOnOff_ROn_spikers_3 + (-1,)]
            SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-1,)] = SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-1,)] + SOnOff_ROn_PSC_fF *(SOnOff_ROn_PSC_maxF - SOnOff_ROn_PSC_F[SOnOff_ROn_spikers_3 + (-1,)])
            SOnOff_ROn_PSC_P[SOnOff_ROn_spikers_3 + (-1,)] = SOnOff_ROn_PSC_P[SOnOff_ROn_spikers_3 + (-1,)] * (1 - SOnOff_ROn_PSC_fP)
        #Compute Gradients

        grad_On_ROn_accumulate += eligibility_On_ROn
        grad_On_ROn_accumulate2 += eligibility_On_ROn
        grad_Off_ROn_accumulate += eligibility_Off_ROn
        grad_Off_ROn_accumulate2 += eligibility_Off_ROn
        grad_On_SOnOff_accumulate += eligibility_On_SOnOff
        grad_On_SOnOff_accumulate2 += eligibility_On_SOnOff
        grad_Off_SOnOff_accumulate += eligibility_Off_SOnOff
        grad_Off_SOnOff_accumulate2 += eligibility_Off_SOnOff
        grad_SOnOff_ROn_accumulate += eligibility_SOnOff_ROn
        grad_SOnOff_ROn_accumulate2 += eligibility_SOnOff_ROn
        grad_strf_gain_accumulate += eligibility_strf_gain
        grad_strf_latency_accumulate += eligibility_strf_latency
        grad_output_adaptation_accumulate += eligibility_output_adaptation
        grad_strf_gain_accumulate2 += eligibility_strf_gain
        grad_strf_latency_accumulate2 += eligibility_strf_latency
        grad_output_adaptation_accumulate2 += eligibility_output_adaptation

        #grab loss


        if timestep % loss_bin_width == 0 and timestep != 0:
            sim_bin = cp.sum(cp.squeeze(cp.sum(ROn_spikes_holder[:,:,:,:,timestep-loss_bin_width:timestep], axis=-1)),axis=-1)

            data_bin = cp.sum(cp.sum(cp.squeeze(data[:,:,timestep-loss_bin_width:timestep,:], axis=-1), axis=-1),axis=-1)

            loss_deriv = 2.0*(sim_bin - data_bin[:,None])
            grad_On_ROn += grad_On_ROn_accumulate*loss_deriv[:,:,None,None]
            grad_On_ROn_accumulate = 0
            grad_Off_ROn += grad_Off_ROn_accumulate*loss_deriv[:,:,None,None]
            grad_Off_ROn_accumulate = 0
            grad_On_SOnOff += grad_On_SOnOff_accumulate*Bk*loss_deriv[:,:,None,None]
            grad_On_SOnOff_accumulate = 0
            grad_Off_SOnOff += grad_Off_SOnOff_accumulate*Bk*loss_deriv[:,:,None,None]
            grad_Off_SOnOff_accumulate = 0
            grad_SOnOff_ROn += grad_SOnOff_ROn_accumulate*loss_deriv[:,:,None,None]
            grad_SOnOff_ROn_accumulate = 0
            grad_strf_gain += grad_strf_gain_accumulate*loss_deriv[:,:,None,None]
            grad_strf_gain_accumulate = 0
            grad_strf_latency += grad_strf_latency_accumulate*loss_deriv[:,:,None,None]
            grad_strf_latency_accumulate = 0
            grad_output_adaptation += grad_output_adaptation_accumulate*loss_deriv[:,:,None,None]
            grad_output_adaptation_accumulate = 0

        if timestep % loss_bin_width2 == 0 and timestep != 0:
            sim_bin = cp.sum(cp.squeeze(cp.sum(ROn_spikes_holder[:,:,:,:,timestep-loss_bin_width2:timestep], axis=-1)),axis=-1)

            data_bin = cp.sum(cp.sum(cp.squeeze(data[:,:,timestep-loss_bin_width2:timestep,:], axis=-1), axis=-1),axis=-1)

            loss_deriv = 2.0*(sim_bin - data_bin[:,None])
            grad_On_ROn += grad_On_ROn_accumulate2*loss_deriv[:,:,None,None]
            grad_On_ROn_accumulate2 = 0
            grad_Off_ROn += grad_Off_ROn_accumulate2*loss_deriv[:,:,None,None]
            grad_Off_ROn_accumulate2 = 0
            grad_On_SOnOff += grad_On_SOnOff_accumulate2*Bk*loss_deriv[:,:,None,None]
            grad_On_SOnOff_accumulate2 = 0
            grad_Off_SOnOff += grad_Off_SOnOff_accumulate2*Bk*loss_deriv[:,:,None,None]
            grad_Off_SOnOff_accumulate2 = 0
            grad_SOnOff_ROn += grad_SOnOff_ROn_accumulate2*loss_deriv[:,:,None,None]
            grad_SOnOff_ROn_accumulate2 = 0
            grad_strf_gain += grad_strf_gain_accumulate2*loss_deriv[:,:,None,None]
            grad_strf_gain_accumulate2 = 0
            grad_strf_latency += grad_strf_latency_accumulate2*loss_deriv[:,:,None,None]
            grad_strf_latency_accumulate2 = 0
            grad_output_adaptation += grad_output_adaptation_accumulate2*loss_deriv[:,:,None,None]
            grad_output_adaptation_accumulate2 = 0
    grads = cp.sum(cp.stack([grad_strf_gain,grad_strf_latency,grad_output_adaptation,grad_On_ROn,grad_Off_ROn,grad_On_SOnOff,grad_Off_SOnOff,grad_SOnOff_ROn], axis = 0), axis = 3)



    return cp.asnumpy(ROn_spikes_holder), cp.asnumpy(grads), cp.asnumpy(On_SOnOff_PSC_s_holder), cp.asnumpy(Off_SOnOff_PSC_s_holder), cp.asnumpy(losses_holder)