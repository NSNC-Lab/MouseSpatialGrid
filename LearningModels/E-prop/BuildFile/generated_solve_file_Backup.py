# pythran export solve_run(float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:]) -> Tuple[float64[:,:,:,:], float64[:,:,:]]
import numpy as np
import matplotlib.pyplot as plt
from BuildFile import calculate_loss_eprop
def solve_run(on_input,off_input,noise_token,data,p):


    #Declare Variables

    On_C = 0.1
    On_g_L = 0.005
    On_E_L = -65
    On_noise = 0
    On_t_ref = 1
    On_E_k = -80
    On_tau_ad = 5
    On_g_inc = 0.0
    On_Itonic = 0
    On_Imask = np.ones((1,1))
    On_R = 200.0
    On_tau = 20.0
    On_V_thresh = -47
    On_V_reset = -54
    is_output = 0
    is_noise = 0
    is_input = 1
    On_g_postIC = 0.17
    On_E_exc = 0
    On_netcon = np.eye(1)
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
    Off_tau_ad = 5
    Off_g_inc = 0.0
    Off_Itonic = 0
    Off_Imask = np.ones((1,1))
    Off_R = 200.0
    Off_tau = 20.0
    Off_V_thresh = -47
    Off_V_reset = -54
    is_output = 0
    is_noise = 0
    is_input = 1
    Off_g_postIC = 0.17
    Off_E_exc = 0
    Off_netcon = np.eye(1)
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
    SOnOff_tau_ad = 5
    SOnOff_g_inc = 0.0
    SOnOff_Itonic = 0
    SOnOff_Imask = np.ones((1,1))
    SOnOff_R = 100.0
    SOnOff_tau = 10.0
    SOnOff_V_thresh = -47
    SOnOff_V_reset = -52
    is_output = 0
    is_noise = 0
    is_input = 0
    SOnOff_g_postIC = 0.17
    SOnOff_E_exc = 0
    SOnOff_netcon = np.eye(1)
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
    ROn_g_inc = 0.0
    ROn_Itonic = 0
    ROn_Imask = np.ones((1,1))
    ROn_R = 200.0
    ROn_tau = 20.0
    ROn_V_thresh = -47
    ROn_V_reset = -54
    is_output = 1
    is_noise = 1
    is_input = 0
    ROn_g_postIC = 0.17
    ROn_E_exc = 0
    ROn_netcon = np.eye(1)
    ROn_nSYN = 0.011
    ROn_noise_E_exc = 0
    ROn_tauR_N = 0.7
    ROn_tauD_N = 1.5
    ROn_noise_scale = 1.9481350796278847
    final_grad_node = 1
    On_ROn_ESYN = 0
    On_ROn_tauD = 1.5
    On_ROn_tauR = 0.7
    On_ROn_PSC_delay = 3
    On_ROn_gSYN = p[0].reshape(100, 1, 1)
    On_ROn_PSC_fF = 0
    On_ROn_PSC_fP = 0.1
    On_ROn_tauF = 180
    On_ROn_tauP = 30
    On_ROn_PSC_maxF = 4
    On_ROn_netcon = np.eye(1)
    On_ROn_scale = 1.9481350796278847
    On_SOnOff_ESYN = 0
    On_SOnOff_tauD = 1
    On_SOnOff_tauR = 0.1
    On_SOnOff_PSC_delay = 1
    On_SOnOff_gSYN = 0.085
    On_SOnOff_PSC_fF = 0
    On_SOnOff_PSC_fP = 0.2
    On_SOnOff_tauF = 180
    On_SOnOff_tauP = 80
    On_SOnOff_PSC_maxF = 4
    On_SOnOff_netcon = np.eye(1)
    On_SOnOff_scale = 1.2915496650148839
    Off_SOnOff_ESYN = 0
    Off_SOnOff_tauD = 1
    Off_SOnOff_tauR = 0.1
    Off_SOnOff_PSC_delay = 1
    Off_SOnOff_gSYN = 0.045
    Off_SOnOff_PSC_fF = 0
    Off_SOnOff_PSC_fP = 0.0
    Off_SOnOff_tauF = 180
    Off_SOnOff_tauP = 80
    Off_SOnOff_PSC_maxF = 4
    Off_SOnOff_netcon = np.eye(1)
    Off_SOnOff_scale = 1.2915496650148839
    SOnOff_ROn_ESYN = -80
    SOnOff_ROn_tauD = 4.5
    SOnOff_ROn_tauR = 1
    SOnOff_ROn_PSC_delay = 0.5
    SOnOff_ROn_gSYN = 0.025
    SOnOff_ROn_PSC_fF = 0
    SOnOff_ROn_PSC_fP = 0.5
    SOnOff_ROn_tauF = 180
    SOnOff_ROn_tauP = 120
    SOnOff_ROn_PSC_maxF = 4
    SOnOff_ROn_netcon = np.eye(1)
    SOnOff_ROn_scale = 1.5368523544529802
    loss_vals = np.array([0])
    loss_bin_width = 1

    tau_IC = 1.0  # ms (try 2–5ms)
    
    # Eligibility trace parameters
    tau_eligibility = 20.0  # eligibility trace decay time constant (ms)
    tau_learning = 50.0     # learning signal smoothing time constant (ms)
    decay_e = np.exp(-0.1/tau_eligibility)
    decay_L = np.exp(-0.1/tau_learning)

    #Declare Holders
    On_IC_g = np.zeros((100,10,1,2))
    Off_IC_g = np.zeros((100,10,1,2))
    
    On_V = np.ones((100,10,1,2)) * On_E_L
    On_g_ad = np.zeros((100,10,1,2))
    On_tspike = np.ones((100,10,1,5)) * -30
    On_buffer_index = np.ones((100,10,1))
    spike_wrt_tau_ad_On = np.zeros((100,10,1))
    spike_wrt_g_inc_On = np.zeros((100,10,1))
    Off_V = np.ones((100,10,1,2)) * Off_E_L
    Off_g_ad = np.zeros((100,10,1,2))
    Off_tspike = np.ones((100,10,1,5)) * -30
    Off_buffer_index = np.ones((100,10,1))
    spike_wrt_tau_ad_Off = np.zeros((100,10,1))
    spike_wrt_g_inc_Off = np.zeros((100,10,1))
    SOnOff_V = np.ones((100,10,1,2)) * SOnOff_E_L
    SOnOff_g_ad = np.zeros((100,10,1,2))
    SOnOff_tspike = np.ones((100,10,1,5)) * -30
    SOnOff_buffer_index = np.ones((100,10,1))
    spike_wrt_tau_ad_SOnOff = np.zeros((100,10,1))
    spike_wrt_g_inc_SOnOff = np.zeros((100,10,1))
    ROn_V = np.ones((100,10,1,2)) * ROn_E_L
    ROn_g_ad = np.zeros((100,10,1,2))
    ROn_tspike = np.ones((100,10,1,5)) * -30
    ROn_buffer_index = np.ones((100,10,1))
    ROn_spikes_holder = np.zeros((100,10,1,29801), dtype=np.int8)
    On_SOnOff_PSC_s_holder = np.zeros((100,10,1,29801), dtype=np.int8)
    Off_SOnOff_PSC_s_holder = np.zeros((100,10,1,29801), dtype=np.int8)
    losses_holder = np.zeros((100,1,1,29801), dtype=np.int8)
    ROn_noise_sn = np.zeros((100,10,1,2))
    ROn_noise_xn = np.zeros((100,10,1,2))
    spike_wrt_tau_ad_ROn = np.zeros((100,10,1))
    spike_wrt_g_inc_ROn = np.zeros((100,10,1))
    On_ROn_PSC_s = np.zeros((100,10,1,2))
    On_ROn_PSC_x = np.zeros((100,10,1,2))
    On_ROn_PSC_F = np.ones((100,10,1,2))
    On_ROn_PSC_P = np.ones((100,10,1,2))
    On_ROn_PSC_q = np.ones((100,10,1,2))
    spike_wrt_gsyn_On_ROn_accumulate = np.zeros((100,10,1,loss_bin_width))
    spike_wrt_gsyn_On_ROn = np.zeros((100,10,1))
    On_SOnOff_PSC_s = np.zeros((100,10,1,2))
    On_SOnOff_PSC_x = np.zeros((100,10,1,2))
    On_SOnOff_PSC_F = np.ones((100,10,1,2))
    On_SOnOff_PSC_P = np.ones((100,10,1,2))
    On_SOnOff_PSC_q = np.ones((100,10,1,2))
    spike_wrt_gsyn_On_SOnOff_accumulate = np.zeros((100,10,1,loss_bin_width))
    spike_wrt_gsyn_On_SOnOff = np.zeros((100,10,1))
    Off_SOnOff_PSC_s = np.zeros((100,10,1,2))
    Off_SOnOff_PSC_x = np.zeros((100,10,1,2))
    Off_SOnOff_PSC_F = np.ones((100,10,1,2))
    Off_SOnOff_PSC_P = np.ones((100,10,1,2))
    Off_SOnOff_PSC_q = np.ones((100,10,1,2))
    spike_wrt_gsyn_Off_SOnOff_accumulate = np.zeros((100,10,1,loss_bin_width))
    spike_wrt_gsyn_Off_SOnOff = np.zeros((100,10,1))
    
    # Eligibility traces - decay over time, accumulate gradient contributions
    eligibility_On_SOnOff = np.zeros((100,10,1))
    eligibility_Off_SOnOff = np.zeros((100,10,1))
    eligibility_SOnOff_ROn = np.zeros((100,10,1))
    eligibility_On_ROn = np.zeros((100,10,1))
    
    # Smoothed activity for learning signal computation
    target_rate_smooth = np.zeros((100,10,1))
    actual_rate_smooth = np.zeros((100,10,1))
    
    SOnOff_ROn_PSC_s = np.zeros((100,10,1,2))
    SOnOff_ROn_PSC_x = np.zeros((100,10,1,2))
    SOnOff_ROn_PSC_F = np.ones((100,10,1,2))
    SOnOff_ROn_PSC_P = np.ones((100,10,1,2))
    SOnOff_ROn_PSC_q = np.ones((100,10,1,2))
    spike_wrt_gsyn_SOnOff_ROn_accumulate = np.zeros((100,10,1,loss_bin_width))
    spike_wrt_gsyn_SOnOff_ROn = np.zeros((100,10,1))

    # Tracking arrays for plotting
    n_timesteps = 29801
    L_t_track = np.zeros(n_timesteps)
    eligibility_On_ROn_track = np.zeros(n_timesteps)
    grad_On_ROn_cumulative = np.zeros((100, n_timesteps))  # track all batches

    for timestep,t in enumerate(np.arange(0,29801*0.1-0.1,0.1)):


        #Declare ODES

        if np.sum(on_input[:,timestep,:])>0:
            x = 1

        #On_IC_g[:,:,:,-2]  = On_IC_g[:,:,:,-1]
        #Off_IC_g[:,:,:,-2] = Off_IC_g[:,:,:,-1]

        #u_on  = on_input[:, timestep, :]  
        #u_off = off_input[:, timestep, :]

        #On_IC_g[:,:,:,-1]  = On_IC_g[:,:,:,-2]  + (0.1*(-On_IC_g[:,:,:,-2]/tau_IC)  + u_on*On_netcon)
        #Off_IC_g[:,:,:,-1] = Off_IC_g[:,:,:,-2] + (0.1*(-Off_IC_g[:,:,:,-2]/tau_IC) + u_off*Off_netcon)

        

        On_V_k1 = (((On_E_L - On_V[:,:,:,-1]) - On_R*On_g_ad[:,:,:,-1]*(On_V[:,:,:,-1]-On_E_k) - On_R*On_g_postIC*on_input[:, timestep, :]  *On_netcon*(On_V[:,:,:,-1]-On_E_exc) + On_R*On_Itonic*On_Imask) / On_tau)
        On_g_ad_k1 = -On_g_ad[:,:,:,-1] / On_tau_ad
        Off_V_k1 = (((Off_E_L - Off_V[:,:,:,-1]) - Off_R*Off_g_ad[:,:,:,-1]*(Off_V[:,:,:,-1]-Off_E_k) - Off_R*Off_g_postIC*off_input[:, timestep, :]  *Off_netcon*(Off_V[:,:,:,-1]-Off_E_exc) + Off_R*Off_Itonic*Off_Imask) / Off_tau)
        Off_g_ad_k1 = -Off_g_ad[:,:,:,-1] / Off_tau_ad
        SOnOff_V_k1 = (((SOnOff_E_L - SOnOff_V[:,:,:,-1]) - SOnOff_R*SOnOff_g_ad[:,:,:,-1]*(SOnOff_V[:,:,:,-1]-SOnOff_E_k) - SOnOff_R*(On_SOnOff_gSYN*On_SOnOff_PSC_s[:,:,:,-1]*On_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-On_SOnOff_ESYN) +Off_SOnOff_gSYN*Off_SOnOff_PSC_s[:,:,:,-1]*Off_SOnOff_netcon*(SOnOff_V[:,:,:,-1]-Off_SOnOff_ESYN) ) + SOnOff_R*SOnOff_Itonic*SOnOff_Imask) / SOnOff_tau)
        SOnOff_g_ad_k1 = -SOnOff_g_ad[:,:,:,-1] / SOnOff_tau_ad
        ROn_V_k1 = (((ROn_E_L - ROn_V[:,:,:,-1]) - ROn_R*ROn_g_ad[:,:,:,-1]*(ROn_V[:,:,:,-1]-ROn_E_k) - ROn_R*(On_ROn_gSYN*On_ROn_PSC_s[:,:,:,-1]*On_ROn_netcon*(ROn_V[:,:,:,-1]-On_ROn_ESYN) +SOnOff_ROn_gSYN*SOnOff_ROn_PSC_s[:,:,:,-1]*SOnOff_ROn_netcon*(ROn_V[:,:,:,-1]-SOnOff_ROn_ESYN) ) + ROn_R*ROn_Itonic*ROn_Imask) / ROn_tau) + ((-ROn_R * ROn_nSYN * ROn_noise_sn[:,:,:,-1]*(ROn_V[:,:,:,-1]-ROn_noise_E_exc)) / ROn_tau)
        ROn_noise_sn_k1 = (ROn_noise_scale * ROn_noise_xn[:,:,:,-1] - ROn_noise_sn[:,:,:,-1]) / ROn_tauR_N
        ROn_noise_xn_k1 = -(ROn_noise_xn[:,:,:,-1]/ROn_tauD_N) + noise_token[:,:,timestep,:]/0.1
        ROn_g_ad_k1 = -ROn_g_ad[:,:,:,-1] / ROn_tau_ad
        On_ROn_PSC_s_k1 = (On_ROn_scale*On_ROn_PSC_x[:,:,:,-1] - On_ROn_PSC_s[:,:,:,-1]) / On_ROn_tauR
        On_ROn_PSC_x_k1 = -On_ROn_PSC_x[:,:,:,-1]/On_ROn_tauD
        On_ROn_PSC_F_k1 = (1 - On_ROn_PSC_F[:,:,:,-1])/On_ROn_tauF
        On_ROn_PSC_P_k1 = (1 - On_ROn_PSC_P[:,:,:,-1])/On_ROn_tauP
        On_ROn_PSC_q_k1 = 0
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
        voltage_wrt_psc_On_ROn = (-ROn_R  * On_ROn_gSYN*On_ROn_netcon*(ROn_V[:,:,:,-2]-On_ROn_ESYN))/ROn_tau/10
        voltage_wrt_psc_On_SOnOff = (-SOnOff_R  * On_SOnOff_gSYN*On_SOnOff_netcon*(SOnOff_V[:,:,:,-2]-On_SOnOff_ESYN))/SOnOff_tau/10
        voltage_wrt_psc_Off_SOnOff = (-SOnOff_R  * Off_SOnOff_gSYN*Off_SOnOff_netcon*(SOnOff_V[:,:,:,-2]-Off_SOnOff_ESYN))/SOnOff_tau/10
        voltage_wrt_psc_SOnOff_ROn = (-ROn_R  * SOnOff_ROn_gSYN*SOnOff_ROn_netcon*(ROn_V[:,:,:,-2]-SOnOff_ROn_ESYN))/ROn_tau/10

        #psc_wrt_spiking_On_ROn = ( 0.1*(On_ROn_scale*((On_ROn_PSC_x[:,:,:,-1] + On_ROn_PSC_q[:,:,:,-1])*-2*np.tanh(t-(On_tspike[:,:,:,-1]+On_ROn_PSC_delay))*(1-(np.tanh(t-(On_tspike[:,:,:,-1]+On_ROn_PSC_delay)))**2)))  )*100 
        psc_wrt_spiking_On_ROn = 0.1*(On_ROn_scale/On_ROn_tauR)*On_ROn_PSC_q[:,:,:,-1]*-2*np.tanh(t-(On_tspike[:,:,:,-1]+On_ROn_PSC_delay))*(1-(np.tanh(t-(On_tspike[:,:,:,-1]+On_ROn_PSC_delay)))**2)
        psc_wrt_spiking_On_SOnOff = 0.1*(On_SOnOff_scale/On_SOnOff_tauR)*On_SOnOff_PSC_q[:,:,:,-1]*-2*np.tanh(t-(On_tspike[:,:,:,-1]+On_SOnOff_PSC_delay))*(1-(np.tanh(t-(On_tspike[:,:,:,-1]+On_SOnOff_PSC_delay)))**2)
        psc_wrt_spiking_Off_SOnOff = 0.1*(Off_SOnOff_scale/Off_SOnOff_tauR)*Off_SOnOff_PSC_q[:,:,:,-1]*-2*np.tanh(t-(Off_tspike[:,:,:,-1]+Off_SOnOff_PSC_delay))*(1-(np.tanh(t-(Off_tspike[:,:,:,-1]+Off_SOnOff_PSC_delay)))**2)
        psc_wrt_spiking_SOnOff_ROn = 0.1*(SOnOff_ROn_scale/SOnOff_ROn_tauR)*SOnOff_ROn_PSC_q[:,:,:,-1]*-2*np.tanh(t-(SOnOff_tspike[:,:,:,-1]+SOnOff_ROn_PSC_delay))*(1-(np.tanh(t-(SOnOff_tspike[:,:,:,-1]+SOnOff_ROn_PSC_delay)))**2)
        
        #psc_wrt_spiking_On_SOnOff = ( 0.1*(On_SOnOff_scale*((On_SOnOff_PSC_x[:,:,:,-1] + On_SOnOff_PSC_q[:,:,:,-1])*-2*np.tanh(t-(On_tspike[:,:,:,-1]+On_SOnOff_PSC_delay))*(1-(np.tanh(t-(On_tspike[:,:,:,-1]+On_SOnOff_PSC_delay)))**2)))  )*100 
        #psc_wrt_spiking_Off_SOnOff = ( 0.1*(Off_SOnOff_scale*((Off_SOnOff_PSC_x[:,:,:,-1] + Off_SOnOff_PSC_q[:,:,:,-1])*-2*np.tanh(t-(Off_tspike[:,:,:,-1]+Off_SOnOff_PSC_delay))*(1-(np.tanh(t-(Off_tspike[:,:,:,-1]+Off_SOnOff_PSC_delay)))**2)))  )*100 
        #psc_wrt_spiking_SOnOff_ROn = ( 0.1*(SOnOff_ROn_scale*((SOnOff_ROn_PSC_x[:,:,:,-1] + SOnOff_ROn_PSC_q[:,:,:,-1])*-2*np.tanh(t-(SOnOff_tspike[:,:,:,-1]+SOnOff_ROn_PSC_delay))*(1-(np.tanh(t-(SOnOff_tspike[:,:,:,-1]+SOnOff_ROn_PSC_delay)))**2)))  )*100 

        beta = 2.0   # mV-ish width
        gamma = 0.3  # gain
        tspike_wrt_odeVoltage_On = gamma * (1 - np.tanh((On_V[:,:,:,-1] - On_V_thresh)/beta)**2) / beta
        tspike_wrt_odeVoltage_Off = gamma * (1 - np.tanh((Off_V[:,:,:,-1] - Off_V_thresh)/beta)**2) / beta
        tspike_wrt_odeVoltage_SOnOff = gamma * (1 - np.tanh((SOnOff_V[:,:,:,-1] - SOnOff_V_thresh)/beta)**2) / beta
        tspike_wrt_odeVoltage_ROn = gamma * (1 - np.tanh((ROn_V[:,:,:,-1] - ROn_V_thresh)/beta)**2) / beta

        #tspike_wrt_odeVoltage_On = (  (t - On_tspike[:,:,:,-1])*np.tanh(-(On_V[:,:,:,-2]-On_V_thresh))*(1-np.tanh(On_V[:,:,:,-1]-On_V_thresh)**2)  )
        #tspike_wrt_odeVoltage_Off = (  (t - Off_tspike[:,:,:,-1])*np.tanh(-(Off_V[:,:,:,-2]-Off_V_thresh))*(1-np.tanh(Off_V[:,:,:,-1]-Off_V_thresh)**2)  )
        #tspike_wrt_odeVoltage_SOnOff = (  (t - SOnOff_tspike[:,:,:,-1])*np.tanh(-(SOnOff_V[:,:,:,-2]-SOnOff_V_thresh))*(1-np.tanh(SOnOff_V[:,:,:,-1]-SOnOff_V_thresh)**2)  )
        #tspike_wrt_odeVoltage_ROn = (  (t - ROn_tspike[:,:,:,-1])*np.tanh(-(ROn_V[:,:,:,-2]-ROn_V_thresh))*(1-np.tanh(ROn_V[:,:,:,-1]-ROn_V_thresh)**2)  )


        #Declare Conditionals

        On_mask = ((On_V[:,:,:,-1] >= On_V_thresh) & (On_V[:,:,:,-2] < On_V_thresh))
        On_V[:,:,:,-2] = np.where(On_mask,On_V[:,:,:,-1], On_V[:,:,:,-2])
        On_V[:,:,:,-1] = np.where(On_mask,On_V_reset, On_V[:,:,:,-1])
        On_g_ad[:,:,:,-2] = np.where(On_mask,On_g_ad[:,:,:,-1], On_g_ad[:,:,:,-2])
        On_g_ad[:,:,:,-1] = np.where(On_mask,On_g_ad[:,:,:,-1]+On_g_inc,On_g_ad[:,:,:,-1])
        B_On, Tr_On, N_On = On_mask.shape
        b_On, tr_On, n_On = np.where(On_mask != 0)
        flat_On = (b_On*Tr_On+tr_On) * N_On + n_On
        tspike_flat_On = On_tspike.reshape(B_On*Tr_On*N_On * 5)
        buffer_flat_On = On_buffer_index.reshape(B_On*Tr_On*N_On)
        row_On = ((buffer_flat_On[flat_On]-1) % 5)
        lin_On = (flat_On*5 + row_On).astype(np.int64)
        tspike_flat_On[lin_On] = t
        mask_flat_On = (On_mask.reshape(B_On*Tr_On*N_On)).astype(np.int64)
        buffer_flat_On[:] = ((buffer_flat_On - 1) + mask_flat_On) % 5 + 1
        On_tspike = tspike_flat_On.reshape(B_On,Tr_On,N_On,5)
        On_buffer_index = buffer_flat_On.reshape(B_On,Tr_On,N_On)
        t4On = t + np.zeros_like(On_tspike)
        tref4On = On_t_ref + np.zeros_like(On_tspike)
        cmpOn = t4On <= (On_tspike + tref4On)
        On_mask_ref = np.any(cmpOn, axis=3)
        On_V[:,:,:,-2] = np.where(On_mask_ref,On_V[:,:,:,-1], On_V[:,:,:,-2])
        On_V[:,:,:,-1] = np.where(On_mask_ref, On_V_reset,On_V[:,:,:,-1])
        Off_mask = ((Off_V[:,:,:,-1] >= Off_V_thresh) & (Off_V[:,:,:,-2] < Off_V_thresh))
        Off_V[:,:,:,-2] = np.where(Off_mask,Off_V[:,:,:,-1], Off_V[:,:,:,-2])
        Off_V[:,:,:,-1] = np.where(Off_mask,Off_V_reset, Off_V[:,:,:,-1])
        Off_g_ad[:,:,:,-2] = np.where(Off_mask,Off_g_ad[:,:,:,-1], Off_g_ad[:,:,:,-2])
        Off_g_ad[:,:,:,-1] = np.where(Off_mask,Off_g_ad[:,:,:,-1]+Off_g_inc,Off_g_ad[:,:,:,-1])
        B_Off, Tr_Off, N_Off = Off_mask.shape
        b_Off, tr_Off, n_Off = np.where(Off_mask != 0)
        flat_Off = (b_Off*Tr_Off+tr_Off) * N_Off + n_Off
        tspike_flat_Off = Off_tspike.reshape(B_Off*Tr_Off*N_Off * 5)
        buffer_flat_Off = Off_buffer_index.reshape(B_Off*Tr_Off*N_Off)
        row_Off = ((buffer_flat_Off[flat_Off]-1) % 5)
        lin_Off = (flat_Off*5 + row_Off).astype(np.int64)
        tspike_flat_Off[lin_Off] = t
        mask_flat_Off = (Off_mask.reshape(B_Off*Tr_Off*N_Off)).astype(np.int64)
        buffer_flat_Off[:] = ((buffer_flat_Off - 1) + mask_flat_Off) % 5 + 1
        Off_tspike = tspike_flat_Off.reshape(B_Off,Tr_Off,N_Off,5)
        Off_buffer_index = buffer_flat_Off.reshape(B_Off,Tr_Off,N_Off)
        t4Off = t + np.zeros_like(Off_tspike)
        tref4Off = Off_t_ref + np.zeros_like(Off_tspike)
        cmpOff = t4Off <= (Off_tspike + tref4Off)
        Off_mask_ref = np.any(cmpOff, axis=3)
        Off_V[:,:,:,-2] = np.where(Off_mask_ref,Off_V[:,:,:,-1], Off_V[:,:,:,-2])
        Off_V[:,:,:,-1] = np.where(Off_mask_ref, Off_V_reset,Off_V[:,:,:,-1])
        SOnOff_mask = ((SOnOff_V[:,:,:,-1] >= SOnOff_V_thresh) & (SOnOff_V[:,:,:,-2] < SOnOff_V_thresh))
        SOnOff_V[:,:,:,-2] = np.where(SOnOff_mask,SOnOff_V[:,:,:,-1], SOnOff_V[:,:,:,-2])
        SOnOff_V[:,:,:,-1] = np.where(SOnOff_mask,SOnOff_V_reset, SOnOff_V[:,:,:,-1])
        SOnOff_g_ad[:,:,:,-2] = np.where(SOnOff_mask,SOnOff_g_ad[:,:,:,-1], SOnOff_g_ad[:,:,:,-2])
        SOnOff_g_ad[:,:,:,-1] = np.where(SOnOff_mask,SOnOff_g_ad[:,:,:,-1]+SOnOff_g_inc,SOnOff_g_ad[:,:,:,-1])
        B_SOnOff, Tr_SOnOff, N_SOnOff = SOnOff_mask.shape
        b_SOnOff, tr_SOnOff, n_SOnOff = np.where(SOnOff_mask != 0)
        flat_SOnOff = (b_SOnOff*Tr_SOnOff+tr_SOnOff) * N_SOnOff + n_SOnOff
        tspike_flat_SOnOff = SOnOff_tspike.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff * 5)
        buffer_flat_SOnOff = SOnOff_buffer_index.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff)
        row_SOnOff = ((buffer_flat_SOnOff[flat_SOnOff]-1) % 5)
        lin_SOnOff = (flat_SOnOff*5 + row_SOnOff).astype(np.int64)
        tspike_flat_SOnOff[lin_SOnOff] = t
        mask_flat_SOnOff = (SOnOff_mask.reshape(B_SOnOff*Tr_SOnOff*N_SOnOff)).astype(np.int64)
        buffer_flat_SOnOff[:] = ((buffer_flat_SOnOff - 1) + mask_flat_SOnOff) % 5 + 1
        SOnOff_tspike = tspike_flat_SOnOff.reshape(B_SOnOff,Tr_SOnOff,N_SOnOff,5)
        SOnOff_buffer_index = buffer_flat_SOnOff.reshape(B_SOnOff,Tr_SOnOff,N_SOnOff)
        t4SOnOff = t + np.zeros_like(SOnOff_tspike)
        tref4SOnOff = SOnOff_t_ref + np.zeros_like(SOnOff_tspike)
        cmpSOnOff = t4SOnOff <= (SOnOff_tspike + tref4SOnOff)
        SOnOff_mask_ref = np.any(cmpSOnOff, axis=3)
        SOnOff_V[:,:,:,-2] = np.where(SOnOff_mask_ref,SOnOff_V[:,:,:,-1], SOnOff_V[:,:,:,-2])
        SOnOff_V[:,:,:,-1] = np.where(SOnOff_mask_ref, SOnOff_V_reset,SOnOff_V[:,:,:,-1])
        ROn_mask = ((ROn_V[:,:,:,-1] >= ROn_V_thresh) & (ROn_V[:,:,:,-2] < ROn_V_thresh))
        ROn_spikes_holder[:,:,:,timestep] = ROn_mask.astype(np.int8)
        ROn_V[:,:,:,-2] = np.where(ROn_mask,ROn_V[:,:,:,-1], ROn_V[:,:,:,-2])
        ROn_V[:,:,:,-1] = np.where(ROn_mask,ROn_V_reset, ROn_V[:,:,:,-1])
        ROn_g_ad[:,:,:,-2] = np.where(ROn_mask,ROn_g_ad[:,:,:,-1], ROn_g_ad[:,:,:,-2])
        ROn_g_ad[:,:,:,-1] = np.where(ROn_mask,ROn_g_ad[:,:,:,-1]+ROn_g_inc,ROn_g_ad[:,:,:,-1])
        B_ROn, Tr_ROn, N_ROn = ROn_mask.shape
        b_ROn, tr_ROn, n_ROn = np.where(ROn_mask != 0)
        flat_ROn = (b_ROn*Tr_ROn+tr_ROn) * N_ROn + n_ROn
        tspike_flat_ROn = ROn_tspike.reshape(B_ROn*Tr_ROn*N_ROn * 5)
        buffer_flat_ROn = ROn_buffer_index.reshape(B_ROn*Tr_ROn*N_ROn)
        row_ROn = ((buffer_flat_ROn[flat_ROn]-1) % 5)
        lin_ROn = (flat_ROn*5 + row_ROn).astype(np.int64)
        tspike_flat_ROn[lin_ROn] = t
        mask_flat_ROn = (ROn_mask.reshape(B_ROn*Tr_ROn*N_ROn)).astype(np.int64)
        buffer_flat_ROn[:] = ((buffer_flat_ROn - 1) + mask_flat_ROn) % 5 + 1
        ROn_tspike = tspike_flat_ROn.reshape(B_ROn,Tr_ROn,N_ROn,5)
        ROn_buffer_index = buffer_flat_ROn.reshape(B_ROn,Tr_ROn,N_ROn)
        t4ROn = t + np.zeros_like(ROn_tspike)
        tref4ROn = ROn_t_ref + np.zeros_like(ROn_tspike)
        cmpROn = t4ROn <= (ROn_tspike + tref4ROn)
        ROn_mask_ref = np.any(cmpROn, axis=3)
        ROn_V[:,:,:,-2] = np.where(ROn_mask_ref,ROn_V[:,:,:,-1], ROn_V[:,:,:,-2])
        ROn_V[:,:,:,-1] = np.where(ROn_mask_ref, ROn_V_reset,ROn_V[:,:,:,-1])
        tOn_ROn = t + np.zeros_like(On_tspike)
        On_ROn_PSC_delay_cmp = On_ROn_PSC_delay + np.zeros_like(On_tspike)
        cmpOn_ROn = tOn_ROn == (On_tspike + On_ROn_PSC_delay_cmp)
        On_ROn_mask_psc = np.any(cmpOn_ROn, axis=3)
        On_ROn_PSC_x[:,:,:,-2] = np.where(On_ROn_mask_psc,On_ROn_PSC_x[:,:,:,-1], On_ROn_PSC_x[:,:,:,-2])
        On_ROn_PSC_q[:,:,:,-2] = np.where(On_ROn_mask_psc,On_ROn_PSC_q[:,:,:,-1], On_ROn_PSC_q[:,:,:,-2])
        On_ROn_PSC_F[:,:,:,-2] = np.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1], On_ROn_PSC_F[:,:,:,-2])
        On_ROn_PSC_P[:,:,:,-2] = np.where(On_ROn_mask_psc,On_ROn_PSC_P[:,:,:,-1], On_ROn_PSC_P[:,:,:,-2])
        On_ROn_PSC_x[:,:,:,-1] = np.where(On_ROn_mask_psc,On_ROn_PSC_x[:,:,:,-1] + On_ROn_PSC_q[:,:,:,-1], On_ROn_PSC_x[:,:,:,-1])
        On_ROn_PSC_q[:,:,:,-1] = np.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1] * On_ROn_PSC_P[:,:,:,-1], On_ROn_PSC_q[:,:,:,-1])
        On_ROn_PSC_F[:,:,:,-1] = np.where(On_ROn_mask_psc,On_ROn_PSC_F[:,:,:,-1] + On_ROn_PSC_fF * (On_ROn_PSC_maxF - On_ROn_PSC_F[:,:,:,-1]), On_ROn_PSC_F[:,:,:,-1])
        On_ROn_PSC_P[:,:,:,-1] = np.where(On_ROn_mask_psc,On_ROn_PSC_P[:,:,:,-1] * (1 - On_ROn_PSC_fP), On_ROn_PSC_P[:,:,:,-1])
        tOn_SOnOff = t + np.zeros_like(On_tspike)
        On_SOnOff_PSC_delay_cmp = On_SOnOff_PSC_delay + np.zeros_like(On_tspike)
        cmpOn_SOnOff = tOn_SOnOff == (On_tspike + On_SOnOff_PSC_delay_cmp)
        On_SOnOff_mask_psc = np.any(cmpOn_SOnOff, axis=3)
        On_SOnOff_PSC_x[:,:,:,-2] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_x[:,:,:,-1], On_SOnOff_PSC_x[:,:,:,-2])
        On_SOnOff_PSC_q[:,:,:,-2] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_q[:,:,:,-1], On_SOnOff_PSC_q[:,:,:,-2])
        On_SOnOff_PSC_F[:,:,:,-2] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1], On_SOnOff_PSC_F[:,:,:,-2])
        On_SOnOff_PSC_P[:,:,:,-2] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_P[:,:,:,-1], On_SOnOff_PSC_P[:,:,:,-2])
        On_SOnOff_PSC_x[:,:,:,-1] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_x[:,:,:,-1] + On_SOnOff_PSC_q[:,:,:,-1], On_SOnOff_PSC_x[:,:,:,-1])
        On_SOnOff_PSC_q[:,:,:,-1] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1] * On_SOnOff_PSC_P[:,:,:,-1], On_SOnOff_PSC_q[:,:,:,-1])
        On_SOnOff_PSC_F[:,:,:,-1] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_F[:,:,:,-1] + On_SOnOff_PSC_fF * (On_SOnOff_PSC_maxF - On_SOnOff_PSC_F[:,:,:,-1]), On_SOnOff_PSC_F[:,:,:,-1])
        On_SOnOff_PSC_P[:,:,:,-1] = np.where(On_SOnOff_mask_psc,On_SOnOff_PSC_P[:,:,:,-1] * (1 - On_SOnOff_PSC_fP), On_SOnOff_PSC_P[:,:,:,-1])
        tOff_SOnOff = t + np.zeros_like(Off_tspike)
        Off_SOnOff_PSC_delay_cmp = Off_SOnOff_PSC_delay + np.zeros_like(Off_tspike)
        cmpOff_SOnOff = tOff_SOnOff == (Off_tspike + Off_SOnOff_PSC_delay_cmp)
        Off_SOnOff_mask_psc = np.any(cmpOff_SOnOff, axis=3)
        Off_SOnOff_PSC_x[:,:,:,-2] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_x[:,:,:,-1], Off_SOnOff_PSC_x[:,:,:,-2])
        Off_SOnOff_PSC_q[:,:,:,-2] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_q[:,:,:,-1], Off_SOnOff_PSC_q[:,:,:,-2])
        Off_SOnOff_PSC_F[:,:,:,-2] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1], Off_SOnOff_PSC_F[:,:,:,-2])
        Off_SOnOff_PSC_P[:,:,:,-2] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_P[:,:,:,-1], Off_SOnOff_PSC_P[:,:,:,-2])
        Off_SOnOff_PSC_x[:,:,:,-1] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_x[:,:,:,-1] + Off_SOnOff_PSC_q[:,:,:,-1], Off_SOnOff_PSC_x[:,:,:,-1])
        Off_SOnOff_PSC_q[:,:,:,-1] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1] * Off_SOnOff_PSC_P[:,:,:,-1], Off_SOnOff_PSC_q[:,:,:,-1])
        Off_SOnOff_PSC_F[:,:,:,-1] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_F[:,:,:,-1] + Off_SOnOff_PSC_fF * (Off_SOnOff_PSC_maxF - Off_SOnOff_PSC_F[:,:,:,-1]), Off_SOnOff_PSC_F[:,:,:,-1])
        Off_SOnOff_PSC_P[:,:,:,-1] = np.where(Off_SOnOff_mask_psc,Off_SOnOff_PSC_P[:,:,:,-1] * (1 - Off_SOnOff_PSC_fP), Off_SOnOff_PSC_P[:,:,:,-1])
        tSOnOff_ROn = t + np.zeros_like(SOnOff_tspike)
        SOnOff_ROn_PSC_delay_cmp = SOnOff_ROn_PSC_delay + np.zeros_like(SOnOff_tspike)
        cmpSOnOff_ROn = tSOnOff_ROn == (SOnOff_tspike + SOnOff_ROn_PSC_delay_cmp)
        SOnOff_ROn_mask_psc = np.any(cmpSOnOff_ROn, axis=3)
        SOnOff_ROn_PSC_x[:,:,:,-2] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_x[:,:,:,-1], SOnOff_ROn_PSC_x[:,:,:,-2])
        SOnOff_ROn_PSC_q[:,:,:,-2] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_q[:,:,:,-1], SOnOff_ROn_PSC_q[:,:,:,-2])
        SOnOff_ROn_PSC_F[:,:,:,-2] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1], SOnOff_ROn_PSC_F[:,:,:,-2])
        SOnOff_ROn_PSC_P[:,:,:,-2] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_P[:,:,:,-1], SOnOff_ROn_PSC_P[:,:,:,-2])
        SOnOff_ROn_PSC_x[:,:,:,-1] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_x[:,:,:,-1] + SOnOff_ROn_PSC_q[:,:,:,-1], SOnOff_ROn_PSC_x[:,:,:,-1])
        SOnOff_ROn_PSC_q[:,:,:,-1] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1] * SOnOff_ROn_PSC_P[:,:,:,-1], SOnOff_ROn_PSC_q[:,:,:,-1])
        SOnOff_ROn_PSC_F[:,:,:,-1] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_F[:,:,:,-1] + SOnOff_ROn_PSC_fF * (SOnOff_ROn_PSC_maxF - SOnOff_ROn_PSC_F[:,:,:,-1]), SOnOff_ROn_PSC_F[:,:,:,-1])
        SOnOff_ROn_PSC_P[:,:,:,-1] = np.where(SOnOff_ROn_mask_psc,SOnOff_ROn_PSC_P[:,:,:,-1] * (1 - SOnOff_ROn_PSC_fP), SOnOff_ROn_PSC_P[:,:,:,-1])
        voltage_wrt_gsyn_On_ROn = 0.1*(-ROn_R * On_ROn_PSC_s[:,:,:,-2]*On_ROn_netcon*(ROn_V[:,:,:,-2]-On_ROn_ESYN))/ROn_tau





        spiking_wrt_voltage_ROn = tspike_wrt_odeVoltage_ROn



        #Parameter Updates


        spike_wrt_gsyn_On_ROn_accumulate[:,:,:,timestep%loss_bin_width] = spiking_wrt_voltage_ROn*voltage_wrt_gsyn_On_ROn


        #grab loss


        #if timestep % loss_bin_width == (loss_bin_width-1) and timestep != 0:
        #    loss_vals = calculate_loss_eprop.calculate(ROn_spikes_holder,timestep,loss_bin_width)
        #    losses_holder[:,:,:,timestep] = loss_vals[:,None,None]
        #    # On_ROn and SOnOff_ROn still use bin-boundary updates (direct connections to output)
        #    spike_wrt_gsyn_On_ROn += np.sum(spike_wrt_gsyn_On_ROn_accumulate,axis=-1)*loss_vals[:,None,None]
        #    spike_wrt_gsyn_SOnOff_ROn += np.sum(spike_wrt_gsyn_SOnOff_ROn_accumulate,axis=-1)*loss_vals[:,None,None]
            # Note: On_SOnOff and Off_SOnOff now use per-timestep eligibility trace updates above
        
        # Track for plotting
        # L_t equivalent: instantaneous error (output - target)
        L_t_instant = ROn_mask.astype(np.float64) - data[:,timestep,:].reshape(1,10,1)
        L_t_track[timestep] = L_t_instant.mean()
        # Eligibility: ψ × ∂V/∂g
        eligibility_On_ROn_track = spiking_wrt_voltage_ROn * voltage_wrt_gsyn_On_ROn
        # Cumulative gradient
        
        spike_wrt_gsyn_On_ROn += L_t_instant*eligibility_On_ROn_track
        
        grad_On_ROn_cumulative[:, timestep] = spike_wrt_gsyn_On_ROn[:, :, 0].mean(axis=1)  # avg over trials

    grads = np.sum(np.stack([spike_wrt_gsyn_On_ROn], axis = 0), axis = 2)

    # # ============== PLOTTING ==============
    # time_ms = np.arange(n_timesteps) * 0.1
    
    # fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    # fig.suptitle('Chain Rule Gradient Analysis: On→ROn', fontsize=14)
    
    # # Plot 1: L_t over time
    # ax = axes[0]
    # ax.plot(time_ms, L_t_track, 'r-', alpha=0.7, linewidth=0.5)
    # ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    # ax.set_ylabel('L_t')
    # ax.set_title('Instantaneous Error L_t = (output - target)')
    # ax.set_xlim([0, 3000])
    
    # # Plot 2: Eligibility over time
    # ax = axes[1]
    # ax.plot(time_ms, eligibility_On_ROn_track, 'g-', alpha=0.7, linewidth=0.5)
    # ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    # ax.set_ylabel('Eligibility')
    # ax.set_title('Eligibility = ψ × ∂V/∂g (On→ROn)')
    # ax.set_xlim([0, 3000])
    
    # # Plot 3: Cumulative gradient for all batches
    # ax = axes[2]
    # for batch in range(100):
    #     ax.plot(time_ms, grad_On_ROn_cumulative[batch, :], alpha=0.3, linewidth=0.5)
    # ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    # ax.set_xlabel('Time (ms)')
    # ax.set_ylabel('Cumulative Gradient')
    # ax.set_title('Cumulative Gradient On→ROn (all 100 batches)')
    # ax.set_xlim([0, 3000])
    
    # plt.tight_layout()
    # plt.savefig('chainrule_grad_analysis.png', dpi=150)
    # #plt.show()
    
    # # ============== RASTER PLOTS ==============
    # fig_raster, axes_raster = plt.subplots(2, 1, figsize=(14, 8))
    # fig_raster.suptitle('Spike Raster: Target vs Simulation Output', fontsize=14)
    
    # # Plot 1: Target raster (from data)
    # ax = axes_raster[0]
    # for trial in range(10):
    #     target_spikes = np.where(data[trial, :, 0] > 0.5)[0] * 0.1
    #     if len(target_spikes) > 0:
    #         ax.scatter(target_spikes, np.ones_like(target_spikes) * trial, 
    #                    c='green', s=1, marker='|')
    # ax.set_ylabel('Trial')
    # ax.set_title('TARGET Raster (from data)')
    # ax.set_xlim([0, 3000])
    # ax.set_ylim([-0.5, 9.5])
    # ax.set_yticks(range(10))
    
    # # Plot 2: Simulated ROn raster (use batch 0 for visualization)
    # ax = axes_raster[1]
    # for trial in range(10):
    #     output_spikes = np.where(ROn_spikes_holder[0, trial, 0, :] > 0.5)[0] * 0.1
    #     if len(output_spikes) > 0:
    #         ax.scatter(output_spikes, np.ones_like(output_spikes) * trial, 
    #                    c='blue', s=1, marker='|')
    # ax.set_xlabel('Time (ms)')
    # ax.set_ylabel('Trial')
    # ax.set_title('SIMULATED ROn Raster (output, batch 0)')
    # ax.set_xlim([0, 3000])
    # ax.set_ylim([-0.5, 9.5])
    # ax.set_yticks(range(10))
    
    # plt.tight_layout()
    # plt.savefig('chainrule_raster_comparison.png', dpi=150)
    # #plt.show()

    return ROn_spikes_holder, grads, On_SOnOff_PSC_s_holder, Off_SOnOff_PSC_s_holder, losses_holder