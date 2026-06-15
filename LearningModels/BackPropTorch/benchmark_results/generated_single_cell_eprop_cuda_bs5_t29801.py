import numpy as np
import cupy as cp
from BuildFile import calculate_loss_eprop
def solve_run(on_input,off_input,noise_token,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p):


    #Transfer inputs to GPU
    on_input = cp.asarray(on_input)
    off_input = cp.asarray(off_input)
    noise_token = cp.asarray(noise_token)
    rate_on = cp.asarray(rate_on)
    rate_off = cp.asarray(rate_off)
    rate_on_deriv = cp.asarray(rate_on_deriv)
    rate_off_deriv = cp.asarray(rate_off_deriv)
    data = cp.asarray(data)
    p = cp.asarray(p)


    #Declare Variables

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
    ROn_t_ref = 4
    ROn_E_k = -80
    ROn_tau_ad = 100
    ROn_g_inc = p[2].reshape(5, 1, 1)
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
    On_ROn_gSYN = p[3].reshape(5, 1, 1)
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
    Off_ROn_gSYN = p[4].reshape(5, 1, 1)
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
    On_SOnOff_gSYN = p[5].reshape(5, 1, 1)
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
    Off_SOnOff_gSYN = p[6].reshape(5, 1, 1)
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
    SOnOff_ROn_gSYN = p[7].reshape(5, 1, 1)
    SOnOff_ROn_PSC_fF = 0
    SOnOff_ROn_PSC_fP = 0.5
    SOnOff_ROn_tauF = 180
    SOnOff_ROn_tauP = 120
    SOnOff_ROn_PSC_maxF = 4
    SOnOff_ROn_netcon = cp.eye(1)
    SOnOff_ROn_scale = 1.5368523544529802
    loss_vals = cp.array([0])
    loss_bin_width = 100
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

    On_V = cp.ones((5,10,1,2)) * On_E_L
    On_g_ad = cp.zeros((5,10,1,2))
    On_tspike = cp.ones((5,10,1,5)) * -30
    On_buffer_index = cp.ones((5,10,1))
    spike_wrt_tau_ad_On = cp.zeros((5,10,1))
    spike_wrt_g_inc_On = cp.zeros((5,10,1))
    Off_V = cp.ones((5,10,1,2)) * Off_E_L
    Off_g_ad = cp.zeros((5,10,1,2))
    Off_tspike = cp.ones((5,10,1,5)) * -30
    Off_buffer_index = cp.ones((5,10,1))
    spike_wrt_tau_ad_Off = cp.zeros((5,10,1))
    spike_wrt_g_inc_Off = cp.zeros((5,10,1))
    SOnOff_V = cp.ones((5,10,1,2)) * SOnOff_E_L
    SOnOff_g_ad = cp.zeros((5,10,1,2))
    SOnOff_tspike = cp.ones((5,10,1,5)) * -30
    SOnOff_buffer_index = cp.ones((5,10,1))
    spike_wrt_tau_ad_SOnOff = cp.zeros((5,10,1))
    spike_wrt_g_inc_SOnOff = cp.zeros((5,10,1))
    ROn_V = cp.ones((5,10,1,2)) * ROn_E_L
    ROn_g_ad = cp.zeros((5,10,1,2))
    ROn_tspike = cp.ones((5,10,1,5)) * -30
    ROn_buffer_index = cp.ones((5,10,1))
    ROn_spikes_holder = cp.zeros((5,10,1,29801), dtype=cp.int8)
    On_SOnOff_PSC_s_holder = cp.zeros((5,10,1,29801), dtype=cp.int8)
    Off_SOnOff_PSC_s_holder = cp.zeros((5,10,1,29801), dtype=cp.int8)
    losses_holder = cp.zeros((5,1,1,29801), dtype=cp.int8)
    ROn_noise_sn = cp.zeros((5,10,1,2))
    ROn_noise_xn = cp.zeros((5,10,1,2))
    spike_wrt_tau_ad_ROn = cp.zeros((5,10,1))
    spike_wrt_g_inc_ROn = cp.zeros((5,10,1))
    On_ROn_PSC_s = cp.zeros((5,10,1,2))
    On_ROn_PSC_x = cp.zeros((5,10,1,2))
    On_ROn_PSC_F = cp.ones((5,10,1,2))
    On_ROn_PSC_P = cp.ones((5,10,1,2))
    On_ROn_PSC_q = cp.ones((5,10,1,2))
    spike_wrt_gsyn_On_ROn_accumulate = cp.zeros((5,10,1,loss_bin_width))
    grad_On_ROn = cp.zeros((5,10,1))
    Off_ROn_PSC_s = cp.zeros((5,10,1,2))
    Off_ROn_PSC_x = cp.zeros((5,10,1,2))
    Off_ROn_PSC_F = cp.ones((5,10,1,2))
    Off_ROn_PSC_P = cp.ones((5,10,1,2))
    Off_ROn_PSC_q = cp.ones((5,10,1,2))
    spike_wrt_gsyn_Off_ROn_accumulate = cp.zeros((5,10,1,loss_bin_width))
    grad_Off_ROn = cp.zeros((5,10,1))
    On_SOnOff_PSC_s = cp.zeros((5,10,1,2))
    On_SOnOff_PSC_x = cp.zeros((5,10,1,2))
    On_SOnOff_PSC_F = cp.ones((5,10,1,2))
    On_SOnOff_PSC_P = cp.ones((5,10,1,2))
    On_SOnOff_PSC_q = cp.ones((5,10,1,2))
    spike_wrt_gsyn_On_SOnOff_accumulate = cp.zeros((5,10,1,loss_bin_width))
    grad_On_SOnOff = cp.zeros((5,10,1))
    Off_SOnOff_PSC_s = cp.zeros((5,10,1,2))
    Off_SOnOff_PSC_x = cp.zeros((5,10,1,2))
    Off_SOnOff_PSC_F = cp.ones((5,10,1,2))
    Off_SOnOff_PSC_P = cp.ones((5,10,1,2))
    Off_SOnOff_PSC_q = cp.ones((5,10,1,2))
    spike_wrt_gsyn_Off_SOnOff_accumulate = cp.zeros((5,10,1,loss_bin_width))
    grad_Off_SOnOff = cp.zeros((5,10,1))
    SOnOff_ROn_PSC_s = cp.zeros((5,10,1,2))
    SOnOff_ROn_PSC_x = cp.zeros((5,10,1,2))
    SOnOff_ROn_PSC_F = cp.ones((5,10,1,2))
    SOnOff_ROn_PSC_P = cp.ones((5,10,1,2))
    SOnOff_ROn_PSC_q = cp.ones((5,10,1,2))
    spike_wrt_gsyn_SOnOff_ROn_accumulate = cp.zeros((5,10,1,loss_bin_width))
    grad_SOnOff_ROn = cp.zeros((5,10,1))
    grad_strf_gain = cp.zeros((5,10,1))
    grad_strf_latency = cp.zeros((5,10,1))
    grad_output_adaptation = cp.zeros((5,10,1))

    for timestep,t in enumerate(cp.arange(0,29801*0.1-0.1,0.1)):


        #Declare ODES

        On_V_k1 = (((On_E_L - On_V[:,:,:,-1]) - On_R*On_g_ad[:,:,:,-1]*(On_V[:,:,:,-1]-On_E_k) - On_R*On_g_postIC*on_input[:,:,:,timestep]*On_netcon*(On_V[:,:,:,-1]-On_E_exc) + On_R*On_Itonic*On_Imask) / On_tau)
        On_g_ad_k1 = -On_g_ad[:,:,:,-1] / On_tau_ad
        Off_V_k1 = (((Off_E_L - Off_V[:,:,:,-1]) - Off_R*Off_g_ad[:,:,:,-1]*(Off_V[:,:,:,-1]-Off_E_k) - Off_R*Off_g_postIC*off_input[:,:,:,timestep]*Off_netcon*(Off_V[:,:,:,-1]-Off_E_exc) + Off_R*Off_Itonic*Off_Imask) / Off_tau)
        Off_g_ad_k1 = -Off_g_ad[:,:,:,-1] / Off_tau_ad
        SOnOff_V_k1 = (((SOnOff_E_L - SOnOff_V[:,:,:,-1]) - SOnOff_R*SOnOff_g_ad[:,:,:,-1]*(SOnOff_V[:,:,:,-1]-SOnOff_E_k) - SOnOff_R*(On_SOnOff_gSYN*On_SOnOff_PSC_s[:,:,:,-1]*On_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-On_SOnOff_ESYN) +Off_SOnOff_gSYN*Off_SOnOff_PSC_s[:,:,:,-1]*Off_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-Off_SOnOff_ESYN) ) + SOnOff_R*SOnOff_Itonic*SOnOff_Imask) / SOnOff_tau)
        SOnOff_g_ad_k1 = -SOnOff_g_ad[:,:,:,-1] / SOnOff_tau_ad
        ROn_V_k1 = (((ROn_E_L - ROn_V[:,:,:,-1]) - ROn_R*ROn_g_ad[:,:,:,-1]*(ROn_V[:,:,:,-1]-ROn_E_k) - ROn_R*(On_ROn_gSYN*On_ROn_PSC_s[:,:,:,-1]*On_ROn_netcon*(ROn_V[:,:,:,-1]-On_ROn_ESYN) +Off_ROn_gSYN*Off_ROn_PSC_s[:,:,:,-1]*Off_ROn_netcon*(ROn_V[:,:,:,-1]-Off_ROn_ESYN) +SOnOff_ROn_gSYN*SOnOff_ROn_PSC_s[:,:,:,-1]*SOnOff_ROn_netcon*(ROn_V[:,:,:,-1]-SOnOff_ROn_ESYN) ) + ROn_R*ROn_Itonic*ROn_Imask) / ROn_tau) + ((-ROn_R * ROn_nSYN * ROn_noise_sn[:,:,:,-1]*(ROn_V[:,:,:,-1]-ROn_noise_E_exc)) / ROn_tau)
        ROn_noise_sn_k1 = (ROn_noise_scale * ROn_noise_xn[:,:,:,-1] - ROn_noise_sn[:,:,:,-1]) / ROn_tauR_N
        ROn_noise_xn_k1 = -(ROn_noise_xn[:,:,:,-1]/ROn_tauD_N) + noise_token[:,:,timestep,:]/0.1
        ROn_g_ad_k1 = -ROn_g_ad[:,:,:,-1] / ROn_tau_ad
        On_ROn_PSC_s_k1 = (On_ROn_scale*On_ROn_PSC_x[:,:,:,-1] - On_ROn_PSC_s[:,:,:,-1]) / On_ROn_tauR
        On_ROn_PSC_x_k1 = -On_ROn_PSC_x[:,:,:,-1]/On_ROn_tauD
        On_ROn_PSC_F_k1 = (1 - On_ROn_PSC_F[:,:,:,-1])/On_ROn_tauF
        On_ROn_PSC_P_k1 = (1 - On_ROn_PSC_P[:,:,:,-1])/On_ROn_tauP
        On_ROn_PSC_q_k1 = 0
        Off_ROn_PSC_s_k1 = (Off_ROn_scale*Off_ROn_PSC_x[:,:,:,-1] - Off_ROn_PSC_s[:,:,:,-1]) / Off_ROn_tauR
        Off_ROn_PSC_x_k1 = -Off_ROn_PSC_x[:,:,:,-1]/Off_ROn_tauD
        Off_ROn_PSC_F_k1 = (1 - Off_ROn_PSC_F[:,:,:,-1])/Off_ROn_tauF
        Off_ROn_PSC_P_k1 = (1 - Off_ROn_PSC_P[:,:,:,-1])/Off_ROn_tauP
        Off_ROn_PSC_q_k1 = 0
        On_SOnOff_PSC_s_k1 = (On_SOnOff_scale*On_SOnOff_PSC_x[:,:,:,-1] - On_SOnOff_PSC_s[:,:,:,-1]) / On_SOnOff_tauR
        On_SOnOff_PSC_x_k1 = -On_SOnOff_PSC_x[:,:,:,-1]/On_SOnOff_tauD
        On_SOnOff_PSC_F_k1 = (1 - On_SOnOff_PSC_F[:,:,:,-1])/On_SOnOff_tauF
        On_SOnOff_PSC_P_k1 = (1 - On_SOnOff_PSC_P[:,:,:,-1])/On_SOnOff_tauP
        On_SOnOff_PSC_q_k1 = 0
        Off_SOnOff_PSC_s_k1 = (Off_SOnOff_scale*Off_SOnOff_PSC_x[:,:,:,-1] - Off_SOnOff_PSC_s[:,:,:,-1]) / Off_SOnOff_tauR
        Off_SOnOff_PSC_x_k1 = -Off_SOnOff_PSC_x[:,:,:,-1]/Off_SOnOff_tauD
        Off_SOnOff_PSC_F_k1 = (1 - Off_SOnOff_PSC_F[:,:,:,-1])/Off_SOnOff_tauF
        Off_SOnOff_PSC_P_k1 = (1 - Off_SOnOff_PSC_P[:,:,:,-1])/Off_SOnOff_tauP
        Off_SOnOff_PSC_q_k1 = 0
        SOnOff_ROn_PSC_s_k1 = (SOnOff_ROn_scale*SOnOff_ROn_PSC_x[:,:,:,-1] - SOnOff_ROn_PSC_s[:,:,:,-1]) / SOnOff_ROn_tauR
        SOnOff_ROn_PSC_x_k1 = -SOnOff_ROn_PSC_x[:,:,:,-1]/SOnOff_ROn_tauD
        SOnOff_ROn_PSC_F_k1 = (1 - SOnOff_ROn_PSC_F[:,:,:,-1])/SOnOff_ROn_tauF
        SOnOff_ROn_PSC_P_k1 = (1 - SOnOff_ROn_PSC_P[:,:,:,-1])/SOnOff_ROn_tauP
        SOnOff_ROn_PSC_q_k1 = 0

        #Declare State Updates

        On_V[:,:,:,-2] = On_V[:,:,:,-1]
        On_V[:,:,:,-1] = On_V[:,:,:,-1] + 0.1*On_V_k1
        On_g_ad[:,:,:,-2] = On_g_ad[:,:,:,-1]
        On_g_ad[:,:,:,-1] = On_g_ad[:,:,:,-1] + 0.1*On_g_ad_k1
        Off_V[:,:,:,-2] = Off_V[:,:,:,-1]
        Off_V[:,:,:,-1] = Off_V[:,:,:,-1] + 0.1*Off_V_k1
        Off_g_ad[:,:,:,-2] = Off_g_ad[:,:,:,-1]
        Off_g_ad[:,:,:,-1] = Off_g_ad[:,:,:,-1] + 0.1*Off_g_ad_k1
        SOnOff_V[:,:,:,-2] = SOnOff_V[:,:,:,-1]
        SOnOff_V[:,:,:,-1] = SOnOff_V[:,:,:,-1] + 0.1*SOnOff_V_k1
        SOnOff_g_ad[:,:,:,-2] = SOnOff_g_ad[:,:,:,-1]
        SOnOff_g_ad[:,:,:,-1] = SOnOff_g_ad[:,:,:,-1] + 0.1*SOnOff_g_ad_k1
        ROn_V[:,:,:,-2] = ROn_V[:,:,:,-1]
        ROn_V[:,:,:,-1] = ROn_V[:,:,:,-1] + 0.1*ROn_V_k1
        ROn_g_ad[:,:,:,-2] = ROn_g_ad[:,:,:,-1]
        ROn_g_ad[:,:,:,-1] = ROn_g_ad[:,:,:,-1] + 0.1*ROn_g_ad_k1
        ROn_noise_sn[:,:,:,-2] = ROn_noise_sn[:,:,:,-1]
        ROn_noise_sn[:,:,:,-1] = ROn_noise_sn[:,:,:,-1] + 0.1*ROn_noise_sn_k1
        ROn_noise_xn[:,:,:,-2] = ROn_noise_xn[:,:,:,-1]
        ROn_noise_xn[:,:,:,-1] = ROn_noise_xn[:,:,:,-1] + 0.1*ROn_noise_xn_k1
        On_ROn_PSC_s[:,:,:,-2] = On_ROn_PSC_s[:,:,:,-1]
        On_ROn_PSC_s[:,:,:,-1] = On_ROn_PSC_s[:,:,:,-1] + 0.1*On_ROn_PSC_s_k1
        On_ROn_PSC_x[:,:,:,-2] = On_ROn_PSC_x[:,:,:,-1]
        On_ROn_PSC_x[:,:,:,-1] = On_ROn_PSC_x[:,:,:,-1] + 0.1*On_ROn_PSC_x_k1
        On_ROn_PSC_F[:,:,:,-2] = On_ROn_PSC_F[:,:,:,-1]
        On_ROn_PSC_F[:,:,:,-1] = On_ROn_PSC_F[:,:,:,-1] + 0.1*On_ROn_PSC_F_k1
        On_ROn_PSC_P[:,:,:,-2] = On_ROn_PSC_P[:,:,:,-1]
        On_ROn_PSC_P[:,:,:,-1] = On_ROn_PSC_P[:,:,:,-1] + 0.1*On_ROn_PSC_P_k1
        On_ROn_PSC_q[:,:,:,-2] = On_ROn_PSC_q[:,:,:,-1]
        On_ROn_PSC_q[:,:,:,-1] = On_ROn_PSC_q[:,:,:,-1] + 0.1*On_ROn_PSC_q_k1
        Off_ROn_PSC_s[:,:,:,-2] = Off_ROn_PSC_s[:,:,:,-1]
        Off_ROn_PSC_s[:,:,:,-1] = Off_ROn_PSC_s[:,:,:,-1] + 0.1*Off_ROn_PSC_s_k1
        Off_ROn_PSC_x[:,:,:,-2] = Off_ROn_PSC_x[:,:,:,-1]
        Off_ROn_PSC_x[:,:,:,-1] = Off_ROn_PSC_x[:,:,:,-1] + 0.1*Off_ROn_PSC_x_k1
        Off_ROn_PSC_F[:,:,:,-2] = Off_ROn_PSC_F[:,:,:,-1]
        Off_ROn_PSC_F[:,:,:,-1] = Off_ROn_PSC_F[:,:,:,-1] + 0.1*Off_ROn_PSC_F_k1
        Off_ROn_PSC_P[:,:,:,-2] = Off_ROn_PSC_P[:,:,:,-1]
        Off_ROn_PSC_P[:,:,:,-1] = Off_ROn_PSC_P[:,:,:,-1] + 0.1*Off_ROn_PSC_P_k1
        Off_ROn_PSC_q[:,:,:,-2] = Off_ROn_PSC_q[:,:,:,-1]
        Off_ROn_PSC_q[:,:,:,-1] = Off_ROn_PSC_q[:,:,:,-1] + 0.1*Off_ROn_PSC_q_k1
        On_SOnOff_PSC_s[:,:,:,-2] = On_SOnOff_PSC_s[:,:,:,-1]
        On_SOnOff_PSC_s[:,:,:,-1] = On_SOnOff_PSC_s[:,:,:,-1] + 0.1*On_SOnOff_PSC_s_k1
        On_SOnOff_PSC_x[:,:,:,-2] = On_SOnOff_PSC_x[:,:,:,-1]
        On_SOnOff_PSC_x[:,:,:,-1] = On_SOnOff_PSC_x[:,:,:,-1] + 0.1*On_SOnOff_PSC_x_k1
        On_SOnOff_PSC_F[:,:,:,-2] = On_SOnOff_PSC_F[:,:,:,-1]
        On_SOnOff_PSC_F[:,:,:,-1] = On_SOnOff_PSC_F[:,:,:,-1] + 0.1*On_SOnOff_PSC_F_k1
        On_SOnOff_PSC_P[:,:,:,-2] = On_SOnOff_PSC_P[:,:,:,-1]
        On_SOnOff_PSC_P[:,:,:,-1] = On_SOnOff_PSC_P[:,:,:,-1] + 0.1*On_SOnOff_PSC_P_k1
        On_SOnOff_PSC_q[:,:,:,-2] = On_SOnOff_PSC_q[:,:,:,-1]
        On_SOnOff_PSC_q[:,:,:,-1] = On_SOnOff_PSC_q[:,:,:,-1] + 0.1*On_SOnOff_PSC_q_k1
        Off_SOnOff_PSC_s[:,:,:,-2] = Off_SOnOff_PSC_s[:,:,:,-1]
        Off_SOnOff_PSC_s[:,:,:,-1] = Off_SOnOff_PSC_s[:,:,:,-1] + 0.1*Off_SOnOff_PSC_s_k1
        Off_SOnOff_PSC_x[:,:,:,-2] = Off_SOnOff_PSC_x[:,:,:,-1]
        Off_SOnOff_PSC_x[:,:,:,-1] = Off_SOnOff_PSC_x[:,:,:,-1] + 0.1*Off_SOnOff_PSC_x_k1
        Off_SOnOff_PSC_F[:,:,:,-2] = Off_SOnOff_PSC_F[:,:,:,-1]
        Off_SOnOff_PSC_F[:,:,:,-1] = Off_SOnOff_PSC_F[:,:,:,-1] + 0.1*Off_SOnOff_PSC_F_k1
        Off_SOnOff_PSC_P[:,:,:,-2] = Off_SOnOff_PSC_P[:,:,:,-1]
        Off_SOnOff_PSC_P[:,:,:,-1] = Off_SOnOff_PSC_P[:,:,:,-1] + 0.1*Off_SOnOff_PSC_P_k1
        Off_SOnOff_PSC_q[:,:,:,-2] = Off_SOnOff_PSC_q[:,:,:,-1]
        Off_SOnOff_PSC_q[:,:,:,-1] = Off_SOnOff_PSC_q[:,:,:,-1] + 0.1*Off_SOnOff_PSC_q_k1
        SOnOff_ROn_PSC_s[:,:,:,-2] = SOnOff_ROn_PSC_s[:,:,:,-1]
        SOnOff_ROn_PSC_s[:,:,:,-1] = SOnOff_ROn_PSC_s[:,:,:,-1] + 0.1*SOnOff_ROn_PSC_s_k1
        SOnOff_ROn_PSC_x[:,:,:,-2] = SOnOff_ROn_PSC_x[:,:,:,-1]
        SOnOff_ROn_PSC_x[:,:,:,-1] = SOnOff_ROn_PSC_x[:,:,:,-1] + 0.1*SOnOff_ROn_PSC_x_k1
        SOnOff_ROn_PSC_F[:,:,:,-2] = SOnOff_ROn_PSC_F[:,:,:,-1]
        SOnOff_ROn_PSC_F[:,:,:,-1] = SOnOff_ROn_PSC_F[:,:,:,-1] + 0.1*SOnOff_ROn_PSC_F_k1
        SOnOff_ROn_PSC_P[:,:,:,-2] = SOnOff_ROn_PSC_P[:,:,:,-1]
        SOnOff_ROn_PSC_P[:,:,:,-1] = SOnOff_ROn_PSC_P[:,:,:,-1] + 0.1*SOnOff_ROn_PSC_P_k1
        SOnOff_ROn_PSC_q[:,:,:,-2] = SOnOff_ROn_PSC_q[:,:,:,-1]
        SOnOff_ROn_PSC_q[:,:,:,-1] = SOnOff_ROn_PSC_q[:,:,:,-1] + 0.1*SOnOff_ROn_PSC_q_k1
        #Compute Phi

        psi_On = (1 - cp.tanh(On_V[:,:,:,-1] - On_V_thresh)**2)
        psi_Off = (1 - cp.tanh(Off_V[:,:,:,-1] - Off_V_thresh)**2)
        psi_ROn = (1 - cp.tanh(ROn_V[:,:,:,-1] - ROn_V_thresh)**2)
        psi_SOnOff = (1 - cp.tanh(SOnOff_V[:,:,:,-1] - SOnOff_V_thresh)**2)

        eligibility_On_ROn = 0.1*(psi_ROn*(-ROn_R) * On_ROn_PSC_s[:,:,:,-1]*On_ROn_netcon*(ROn_V[:,:,:,-1]-On_ROn_ESYN))/ROn_tau
        eligibility_Off_ROn = 0.1*(psi_ROn*(-ROn_R) * Off_ROn_PSC_s[:,:,:,-1]*Off_ROn_netcon*(ROn_V[:,:,:,-1]-Off_ROn_ESYN))/ROn_tau
        eligibility_On_SOnOff = 0.1*(psi_SOnOff*(-SOnOff_R) * On_SOnOff_PSC_s[:,:,:,-1]*On_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-On_SOnOff_ESYN))/SOnOff_tau
        eligibility_Off_SOnOff = 0.1*(psi_SOnOff*(-SOnOff_R) * Off_SOnOff_PSC_s[:,:,:,-1]*Off_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-Off_SOnOff_ESYN))/SOnOff_tau
        eligibility_SOnOff_ROn = 0.1*(psi_ROn*(-ROn_R) * SOnOff_ROn_PSC_s[:,:,:,-1]*SOnOff_ROn_netcon*(ROn_V[:,:,:,-1]-SOnOff_ROn_ESYN))/ROn_tau
        eligibility_strf_gain = ((psi_On * (-On_R) * On_g_postIC * 30 * (1 - cp.tanh(0.0001*rate_on[timestep,:][:,None,None]-1.5)**2) * On_netcon * (On_V[:,:,:,-1]-On_E_exc)) / On_tau) + ((psi_Off * (-Off_R) * Off_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_off[timestep,:][:,None,None]-1.5)**2) * Off_netcon * (Off_V[:,:,:,-1]-Off_E_exc)) / Off_tau)
        eligibility_strf_latency = ((psi_On * (-On_R) * On_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_on_deriv[timestep,:][:,None,None]-1.5)**2) * On_netcon * (On_V[:,:,:,-1]-On_E_exc)) / On_tau) + ((psi_Off * (-Off_R) * Off_g_postIC * 30 *(1 - cp.tanh(0.0001*rate_off_deriv[timestep,:][:,None,None]-1.5)**2) * Off_netcon * (Off_V[:,:,:,-1]-Off_E_exc)) / Off_tau)
        eligibility_output_adaptation = (-ROn_R*(((1+cp.tanh(ROn_V[:,:,:,-1]-ROn_V_thresh))/2)*((1+cp.tanh(-(ROn_V[:,:,:,-2]-ROn_V_thresh)))/2))*0.1*(ROn_V[:,:,:,-1]-ROn_E_k))/ROn_tau

        #Compute Bk

        Bk = -SOnOff_ROn_gSYN*(ROn_V_thresh - SOnOff_ROn_ESYN)


        #Declare Conditionals

        On_mask = ((On_V[:,:,:,-1] >= On_V_thresh) & (On_V[:,:,:,-2] < On_V_thresh))
        On_V[:,:,:,-2] = cp.where(On_mask,On_V[:,:,:,-1], On_V[:,:,:,-2])
        On_V[:,:,:,-1] = cp.where(On_mask,On_V_reset, On_V[:,:,:,-1])
        On_g_ad[:,:,:,-2] = cp.where(On_mask,On_g_ad[:,:,:,-1], On_g_ad[:,:,:,-2])
        On_g_ad[:,:,:,-1] = cp.where(On_mask,On_g_ad[:,:,:,-1]+On_g_inc,On_g_ad[:,:,:,-1])
        B_On, Tr_On, N_On = On_mask.shape
        b_On, tr_On, n_On = cp.where(On_mask != 0)
        flat_On = (b_On*Tr_On+tr_On) * N_On + n_On
        tspike_flat_On = On_tspike.reshape(B_On*Tr_On*N_On * 5)
        buffer_flat_On = On_buffer_index.reshape(B_On*Tr_On*N_On)
        row_On = ((buffer_flat_On[flat_On]-1) % 5)
        lin_On = (flat_On*5 + row_On).astype(cp.int64)
        tspike_flat_On[lin_On] = t
        mask_flat_On = (On_mask.reshape(B_On*Tr_On*N_On)).astype(cp.int64)
        buffer_flat_On[:] = ((buffer_flat_On - 1) + mask_flat_On) % 5 + 1
        On_tspike = tspike_flat_On.reshape(B_On,Tr_On,N_On,5)
        On_buffer_index = buffer_flat_On.reshape(B_On,Tr_On,N_On)
        t4On = t + cp.zeros_like(On_tspike)
        tref4On = On_t_ref + cp.zeros_like(On_tspike)
        cmpOn = t4On <= (On_tspike + tref4On)
        On_mask_ref = cp.any(cmpOn, axis=3)
        On_V[:,:,:,-2] = cp.where(On_mask_ref,On_V[:,:,:,-1], On_V[:,:,:,-2])
        On_V[:,:,:,-1] = cp.where(On_mask_ref, On_V_reset,On_V[:,:,:,-1])
        Off_mask = ((Off_V[:,:,:,-1] >= Off_V_thresh) & (Off_V[:,:,:,-2] < Off_V_thresh))
        Off_V[:,:,:,-2] = cp.where(Off_mask,Off_V[:,:,:,-1], Off_V[:,:,:,-2])
        Off_V[:,:,:,-1] = cp.where(Off_mask,Off_V_reset, Off_V[:,:,:,-1])
        Off_g_ad[:,:,:,-2] = cp.where(Off_mask,Off_g_ad[:,:,:,-1], Off_g_ad[:,:,:,-2])
        Off_g_ad[:,:,:,-1] = cp.where(Off_mask,Off_g_ad[:,:,:,-1]+Off_g_inc,Off_g_ad[:,:,:,-1])
        B_Off, Tr_Off, N_Off = Off_mask.shape
        b_Off, tr_Off, n_Off = cp.where(Off_mask != 0)
        flat_Off = (b_Off*Tr_Off+tr_Off) * N_Off + n_Off
        tspike_flat_Off = Off_tspike.reshape(B_Off*Tr_Off*N_Off * 5)
        buffer_flat_Off = Off_buffer_index.reshape(B_Off*Tr_Off*N_Off)
        row_Off = ((buffer_flat_Off[flat_Off]-1) % 5)
        lin_Off = (flat_Off*5 + row_Off).astype(cp.int64)
        tspike_flat_Off[lin_Off] = t
        mask_flat_Off = (Off_mask.reshape(B_Off*Tr_Off*N_Off)).astype(cp.int64)
        buffer_flat_Off[:] = ((buffer_flat_Off - 1) + mask_flat_Off) % 5 + 1
        Off_tspike = tspike_flat_Off.reshape(B_Off,Tr_Off,N_Off,5)
        Off_buffer_index = buffer_flat_Off.reshape(B_Off,Tr_Off,N_Off)
        t4Off = t + cp.zeros_like(Off_tspike)
        tref4Off = Off_t_ref + cp.zeros_like(Off_tspike)
        cmpOff = t4Off <= (Off_tspike + tref4Off)
        Off_mask_ref = cp.any(cmpOff, axis=3)
        Off_V[:,:,:,-2] = cp.where(Off_mask_ref,Off_V[:,:,:,-1], Off_V[:,:,:,-2])
        Off_V[:,:,:,-1] = cp.where(Off_mask_ref, Off_V_reset,Off_V[:,:,:,-1])
        SOnOff_mask = ((SOnOff_V[:,:,:,-1] >= SOnOff_V_thresh) & (SOnOff_V[:,:,:,-2] < SOnOff_V_thresh))
        SOnOff_V[:,:,:,-2] = cp.where(SOnOff_mask,SOnOff_V[:,:,:,-1], SOnOff_V[:,:,:,-2])
        SOnOff_V[:,:,:,-1] = cp.where(SOnOff_mask,SOnOff_V_reset, SOnOff_V[:,:,:,-1])
        SOnOff_g_ad[:,:,:,-2] = cp.where(SOnOff_mask,SOnOff_g_ad[:,:,:,-1], SOnOff_g_ad[:,:,:,-2])
        SOnOff_g_ad[:,:,:,-1] = cp.where(SOnOff_mask,SOnOff_g_ad[:,:,:,-1]+SOnOff_g_inc,SOnOff_g_ad[:,:,:,-1])
        B_SOnOff, Tr_SOnOff, N_SOnOff = SOnOff_mask.shape
        b_SOnOff, tr_SOnOff, n_SOnOff = cp.where(SOnOff_mask != 0)
        flat_SOnOff = (b_SOnOff*Tr_SOnOff+tr_SOnOff) * N_SOnOff + n_SOnOff
        tspike_flat_SOnOff = SOnOff_tspike.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff * 5)
        buffer_flat_SOnOff = SOnOff_buffer_index.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff)
        row_SOnOff = ((buffer_flat_SOnOff[flat_SOnOff]-1) % 5)
        lin_SOnOff = (flat_SOnOff*5 + row_SOnOff).astype(cp.int64)
        tspike_flat_SOnOff[lin_SOnOff] = t
        mask_flat_SOnOff = (SOnOff_mask.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff)).astype(cp.int64)
        buffer_flat_SOnOff[:] = ((buffer_flat_SOnOff - 1) + mask_flat_SOnOff) % 5 + 1
        SOnOff_tspike = tspike_flat_SOnOff.reshape(B_SOnOff,Tr_SOnOff,N_SOnOff,5)
        SOnOff_buffer_index = buffer_flat_SOnOff.reshape(B_SOnOff,Tr_SOnOff,N_SOnOff)
        t4SOnOff = t + cp.zeros_like(SOnOff_tspike)
        tref4SOnOff = SOnOff_t_ref + cp.zeros_like(SOnOff_tspike)
        cmpSOnOff = t4SOnOff <= (SOnOff_tspike + tref4SOnOff)
        SOnOff_mask_ref = cp.any(cmpSOnOff, axis=3)
        SOnOff_V[:,:,:,-2] = cp.where(SOnOff_mask_ref,SOnOff_V[:,:,:,-1], SOnOff_V[:,:,:,-2])
        SOnOff_V[:,:,:,-1] = cp.where(SOnOff_mask_ref, SOnOff_V_reset,SOnOff_V[:,:,:,-1])
        ROn_mask = ((ROn_V[:,:,:,-1] >= ROn_V_thresh) & (ROn_V[:,:,:,-2] < ROn_V_thresh))
        ROn_spikes_holder[:,:,:,timestep] = ROn_mask.astype(cp.int8)
        ROn_V[:,:,:,-2] = cp.where(ROn_mask,ROn_V[:,:,:,-1], ROn_V[:,:,:,-2])
        ROn_V[:,:,:,-1] = cp.where(ROn_mask,ROn_V_reset, ROn_V[:,:,:,-1])
        ROn_g_ad[:,:,:,-2] = cp.where(ROn_mask,ROn_g_ad[:,:,:,-1], ROn_g_ad[:,:,:,-2])
        ROn_g_ad[:,:,:,-1] = cp.where(ROn_mask,ROn_g_ad[:,:,:,-1]+ROn_g_inc,ROn_g_ad[:,:,:,-1])
        B_ROn, Tr_ROn, N_ROn = ROn_mask.shape
        b_ROn, tr_ROn, n_ROn = cp.where(ROn_mask != 0)
        flat_ROn = (b_ROn*Tr_ROn+tr_ROn) * N_ROn + n_ROn
        tspike_flat_ROn = ROn_tspike.reshape(B_ROn*Tr_ROn*N_ROn * 5)
        buffer_flat_ROn = ROn_buffer_index.reshape(B_ROn*Tr_ROn*N_ROn)
        row_ROn = ((buffer_flat_ROn[flat_ROn]-1) % 5)
        lin_ROn = (flat_ROn*5 + row_ROn).astype(cp.int64)
        tspike_flat_ROn[lin_ROn] = t
        mask_flat_ROn = (ROn_mask.reshape(B_ROn*Tr_ROn*N_ROn)).astype(cp.int64)
        buffer_flat_ROn[:] = ((buffer_flat_ROn - 1) + mask_flat_ROn) % 5 + 1
        ROn_tspike = tspike_flat_ROn.reshape(B_ROn,Tr_ROn,N_ROn,5)
        ROn_buffer_index = buffer_flat_ROn.reshape(B_ROn,Tr_ROn,N_ROn)
        t4ROn = t + cp.zeros_like(ROn_tspike)
        tref4ROn = ROn_t_ref + cp.zeros_like(ROn_tspike)
        cmpROn = t4ROn <= (ROn_tspike + tref4ROn)
        ROn_mask_ref = cp.any(cmpROn, axis=3)
        ROn_V[:,:,:,-2] = cp.where(ROn_mask_ref,ROn_V[:,:,:,-1], ROn_V[:,:,:,-2])
        ROn_V[:,:,:,-1] = cp.where(ROn_mask_ref, ROn_V_reset,ROn_V[:,:,:,-1])
        tOn_ROn = t + cp.zeros_like(On_tspike)
        On_ROn_PSC_delay_cmp = On_ROn_PSC_delay + cp.zeros_like(On_tspike)
        cmpOn_ROn = tOn_ROn == (On_tspike + On_ROn_PSC_delay_cmp)
        On_ROn_mask_psc = cp.any(cmpOn_ROn, axis=3)
        On_ROn_PSC_x[:,:,:,-2] = cp.where(On_ROn_mask_psc,On_ROn_PSC_x[:,:,:,-1], On_ROn_PSC_x[:,:,:,-2])
        On_ROn_PSC_q[:,:,:,-2] = cp.where(On_ROn_mask_psc,On_ROn_PSC_q[:,:,:,-1], On_ROn_PSC_q[:,:,:,-2])
        On_ROn_PSC_F[:,:,:,-2] = cp.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1], On_ROn_PSC_F[:,:,:,-2])
        On_ROn_PSC_P[:,:,:,-2] = cp.where(On_ROn_mask_psc,On_ROn_PSC_P[:,:,:,-1], On_ROn_PSC_P[:,:,:,-2])
        On_ROn_PSC_x[:,:,:,-1] = cp.where(On_ROn_mask_psc,On_ROn_PSC_x[:,:,:,-1] + On_ROn_PSC_q[:,:,:,-1], On_ROn_PSC_x[:,:,:,-1])
        On_ROn_PSC_q[:,:,:,-1] = cp.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1] * On_ROn_PSC_P[:,:,:,-1], On_ROn_PSC_q[:,:,:,-1])
        On_ROn_PSC_F[:,:,:,-1] = cp.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1] + On_ROn_PSC_fF * (On_ROn_PSC_maxF - On_ROn_PSC_F[:,:,:,-1]), On_ROn_PSC_F[:,:,:,-1])
        On_ROn_PSC_P[:,:,:,-1] = cp.where(On_ROn_mask_psc,On_ROn_PSC_P[:,:,:,-1] * (1 - On_ROn_PSC_fP), On_ROn_PSC_P[:,:,:,-1])
        tOff_ROn = t + cp.zeros_like(Off_tspike)
        Off_ROn_PSC_delay_cmp = Off_ROn_PSC_delay + cp.zeros_like(Off_tspike)
        cmpOff_ROn = tOff_ROn == (Off_tspike + Off_ROn_PSC_delay_cmp)
        Off_ROn_mask_psc = cp.any(cmpOff_ROn, axis=3)
        Off_ROn_PSC_x[:,:,:,-2] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_x[:,:,:,-1], Off_ROn_PSC_x[:,:,:,-2])
        Off_ROn_PSC_q[:,:,:,-2] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_q[:,:,:,-1], Off_ROn_PSC_q[:,:,:,-2])
        Off_ROn_PSC_F[:,:,:,-2] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_F[:,:,:,-1], Off_ROn_PSC_F[:,:,:,-2])
        Off_ROn_PSC_P[:,:,:,-2] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_P[:,:,:,-1], Off_ROn_PSC_P[:,:,:,-2])
        Off_ROn_PSC_x[:,:,:,-1] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_x[:,:,:,-1] + Off_ROn_PSC_q[:,:,:,-1], Off_ROn_PSC_x[:,:,:,-1])
        Off_ROn_PSC_q[:,:,:,-1] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_F[:,:,:,-1] * Off_ROn_PSC_P[:,:,:,-1], Off_ROn_PSC_q[:,:,:,-1])
        Off_ROn_PSC_F[:,:,:,-1] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_F[:,:,:,-1] + Off_ROn_PSC_fF * (Off_ROn_PSC_maxF - Off_ROn_PSC_F[:,:,:,-1]), Off_ROn_PSC_F[:,:,:,-1])
        Off_ROn_PSC_P[:,:,:,-1] = cp.where(Off_ROn_mask_psc,Off_ROn_PSC_P[:,:,:,-1] * (1 - Off_ROn_PSC_fP), Off_ROn_PSC_P[:,:,:,-1])
        tOn_SOnOff = t + cp.zeros_like(On_tspike)
        On_SOnOff_PSC_delay_cmp = On_SOnOff_PSC_delay + cp.zeros_like(On_tspike)
        cmpOn_SOnOff = tOn_SOnOff == (On_tspike + On_SOnOff_PSC_delay_cmp)
        On_SOnOff_mask_psc = cp.any(cmpOn_SOnOff, axis=3)
        On_SOnOff_PSC_x[:,:,:,-2] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_x[:,:,:,-1], On_SOnOff_PSC_x[:,:,:,-2])
        On_SOnOff_PSC_q[:,:,:,-2] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_q[:,:,:,-1], On_SOnOff_PSC_q[:,:,:,-2])
        On_SOnOff_PSC_F[:,:,:,-2] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1], On_SOnOff_PSC_F[:,:,:,-2])
        On_SOnOff_PSC_P[:,:,:,-2] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_P[:,:,:,-1], On_SOnOff_PSC_P[:,:,:,-2])
        On_SOnOff_PSC_x[:,:,:,-1] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_x[:,:,:,-1] + On_SOnOff_PSC_q[:,:,:,-1], On_SOnOff_PSC_x[:,:,:,-1])
        On_SOnOff_PSC_q[:,:,:,-1] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1] * On_SOnOff_PSC_P[:,:,:,-1], On_SOnOff_PSC_q[:,:,:,-1])
        On_SOnOff_PSC_F[:,:,:,-1] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1] + On_SOnOff_PSC_fF * (On_SOnOff_PSC_maxF - On_SOnOff_PSC_F[:,:,:,-1]), On_SOnOff_PSC_F[:,:,:,-1])
        On_SOnOff_PSC_P[:,:,:,-1] = cp.where(On_SOnOff_mask_psc,On_SOnOff_PSC_P[:,:,:,-1] * (1 - On_SOnOff_PSC_fP), On_SOnOff_PSC_P[:,:,:,-1])
        tOff_SOnOff = t + cp.zeros_like(Off_tspike)
        Off_SOnOff_PSC_delay_cmp = Off_SOnOff_PSC_delay + cp.zeros_like(Off_tspike)
        cmpOff_SOnOff = tOff_SOnOff == (Off_tspike + Off_SOnOff_PSC_delay_cmp)
        Off_SOnOff_mask_psc = cp.any(cmpOff_SOnOff, axis=3)
        Off_SOnOff_PSC_x[:,:,:,-2] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_x[:,:,:,-1], Off_SOnOff_PSC_x[:,:,:,-2])
        Off_SOnOff_PSC_q[:,:,:,-2] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_q[:,:,:,-1], Off_SOnOff_PSC_q[:,:,:,-2])
        Off_SOnOff_PSC_F[:,:,:,-2] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1], Off_SOnOff_PSC_F[:,:,:,-2])
        Off_SOnOff_PSC_P[:,:,:,-2] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_P[:,:,:,-1], Off_SOnOff_PSC_P[:,:,:,-2])
        Off_SOnOff_PSC_x[:,:,:,-1] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_x[:,:,:,-1] + Off_SOnOff_PSC_q[:,:,:,-1], Off_SOnOff_PSC_x[:,:,:,-1])
        Off_SOnOff_PSC_q[:,:,:,-1] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1] * Off_SOnOff_PSC_P[:,:,:,-1], Off_SOnOff_PSC_q[:,:,:,-1])
        Off_SOnOff_PSC_F[:,:,:,-1] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1] + Off_SOnOff_PSC_fF * (Off_SOnOff_PSC_maxF - Off_SOnOff_PSC_F[:,:,:,-1]), Off_SOnOff_PSC_F[:,:,:,-1])
        Off_SOnOff_PSC_P[:,:,:,-1] = cp.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_P[:,:,:,-1] * (1 - Off_SOnOff_PSC_fP), Off_SOnOff_PSC_P[:,:,:,-1])
        tSOnOff_ROn = t + cp.zeros_like(SOnOff_tspike)
        SOnOff_ROn_PSC_delay_cmp = SOnOff_ROn_PSC_delay + cp.zeros_like(SOnOff_tspike)
        cmpSOnOff_ROn = tSOnOff_ROn == (SOnOff_tspike + SOnOff_ROn_PSC_delay_cmp)
        SOnOff_ROn_mask_psc = cp.any(cmpSOnOff_ROn, axis=3)
        SOnOff_ROn_PSC_x[:,:,:,-2] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_x[:,:,:,-1], SOnOff_ROn_PSC_x[:,:,:,-2])
        SOnOff_ROn_PSC_q[:,:,:,-2] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_q[:,:,:,-1], SOnOff_ROn_PSC_q[:,:,:,-2])
        SOnOff_ROn_PSC_F[:,:,:,-2] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1], SOnOff_ROn_PSC_F[:,:,:,-2])
        SOnOff_ROn_PSC_P[:,:,:,-2] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_P[:,:,:,-1], SOnOff_ROn_PSC_P[:,:,:,-2])
        SOnOff_ROn_PSC_x[:,:,:,-1] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_x[:,:,:,-1] + SOnOff_ROn_PSC_q[:,:,:,-1], SOnOff_ROn_PSC_x[:,:,:,-1])
        SOnOff_ROn_PSC_q[:,:,:,-1] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1] * SOnOff_ROn_PSC_P[:,:,:,-1], SOnOff_ROn_PSC_q[:,:,:,-1])
        SOnOff_ROn_PSC_F[:,:,:,-1] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1] + SOnOff_ROn_PSC_fF * (SOnOff_ROn_PSC_maxF - SOnOff_ROn_PSC_F[:,:,:,-1]), SOnOff_ROn_PSC_F[:,:,:,-1])
        SOnOff_ROn_PSC_P[:,:,:,-1] = cp.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_P[:,:,:,-1] * (1 - SOnOff_ROn_PSC_fP), SOnOff_ROn_PSC_P[:,:,:,-1])
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

        #grab loss


        if timestep % loss_bin_width == 0 and timestep != 0:
            sim_bin = cp.sum(ROn_spikes_holder[:,:,:,timestep-loss_bin_width:timestep], axis=(1,2,3))

            data_bin = cp.sum(data[:,timestep-loss_bin_width:timestep,:])

            loss_deriv = 2.0*(sim_bin - data_bin)
            grad_On_ROn += grad_On_ROn_accumulate*loss_deriv[:,None,None]
            grad_On_ROn_accumulate = 0
            grad_Off_ROn += grad_Off_ROn_accumulate*loss_deriv[:,None,None]
            grad_Off_ROn_accumulate = 0
            grad_On_SOnOff += grad_On_SOnOff_accumulate*Bk*loss_deriv[:,None,None]
            grad_On_SOnOff_accumulate = 0
            grad_Off_SOnOff += grad_Off_SOnOff_accumulate*Bk*loss_deriv[:,None,None]
            grad_Off_SOnOff_accumulate = 0
            grad_SOnOff_ROn += grad_SOnOff_ROn_accumulate*loss_deriv[:,None,None]
            grad_SOnOff_ROn_accumulate = 0
            grad_strf_gain += grad_strf_gain_accumulate*loss_deriv[:,None,None]
            grad_strf_gain_accumulate = 0
            grad_strf_latency += grad_strf_latency_accumulate*loss_deriv[:,None,None]
            grad_strf_latency_accumulate = 0
            grad_output_adaptation += grad_output_adaptation_accumulate*loss_deriv[:,None,None]
            grad_output_adaptation_accumulate = 0
    grads = cp.sum(cp.stack([grad_strf_gain,grad_strf_latency,grad_output_adaptation,grad_On_ROn,grad_Off_ROn,grad_On_SOnOff,grad_Off_SOnOff,grad_SOnOff_ROn], axis = 0), axis = 2)



    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(ROn_spikes_holder), cp.asnumpy(grads), cp.asnumpy(On_SOnOff_PSC_s_holder), cp.asnumpy(Off_SOnOff_PSC_s_holder), cp.asnumpy(losses_holder)