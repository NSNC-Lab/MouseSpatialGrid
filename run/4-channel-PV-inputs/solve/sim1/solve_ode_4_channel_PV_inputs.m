function [T,On_V,On_g_ad,Off_V,Off_g_ad,ROn_V,ROn_g_ad,ROff_V,ROff_g_ad,SOnOff_V,SOnOff_g_ad,TD_V,TD_g_ad,X_V,X_g_ad,C_V,C_g_ad,ROn_On_PSC3_s,ROn_On_PSC3_x,ROn_On_PSC3_F,ROn_On_PSC3_P,ROn_On_PSC3_q,SOnOff_On_PSC_s,SOnOff_On_PSC_x,SOnOff_On_PSC_F,SOnOff_On_PSC_P,SOnOff_On_PSC_q,ROn_SOnOff_PSC3_s,ROn_SOnOff_PSC3_x,ROn_SOnOff_PSC3_F,ROn_SOnOff_PSC3_P,ROn_SOnOff_PSC3_q,ROff_SOnOff_PSC_s,ROff_SOnOff_PSC_x,ROff_SOnOff_PSC_F,ROff_SOnOff_PSC_P,ROff_SOnOff_PSC_q,ROff_Off_PSC_s,ROff_Off_PSC_x,ROff_Off_PSC_F,ROff_Off_PSC_P,ROff_Off_PSC_q,SOnOff_Off_PSC_s,SOnOff_Off_PSC_x,SOnOff_Off_PSC_F,SOnOff_Off_PSC_P,SOnOff_Off_PSC_q,ROn_ROn_iNoise_V3_sn,ROn_ROn_iNoise_V3_xn,X_ROn_PSC_s,X_ROn_PSC_x,X_ROn_PSC_F,X_ROn_PSC_P,X_ROn_PSC_q,ROn_X_PSC3_s,ROn_X_PSC3_x,ROn_X_PSC3_F,ROn_X_PSC3_P,ROn_X_PSC3_q,ROn_TD_PSC_s,ROn_TD_PSC_x,ROn_TD_PSC_F,ROn_TD_PSC_P,ROn_TD_PSC_q,ROff_TD_PSC_s,ROff_TD_PSC_x,ROff_TD_PSC_F,ROff_TD_PSC_P,ROff_TD_PSC_q,X_TD_PSC_s,X_TD_PSC_x,X_TD_PSC_F,X_TD_PSC_P,X_TD_PSC_q,C_ROn_PSC3_s,C_ROn_PSC3_x,C_ROn_PSC3_F,C_ROn_PSC3_P,C_ROn_PSC3_q,On_V_spikes,Off_V_spikes,ROn_V_spikes,ROff_V_spikes,SOnOff_V_spikes,TD_V_spikes,X_V_spikes,C_V_spikes,On_On_IC_iIC,Off_Off_IC_iIC,ROn_On_PSC3_syn,SOnOff_On_PSC_syn,ROn_SOnOff_PSC3_syn,ROff_SOnOff_PSC_syn,ROff_Off_PSC_syn,SOnOff_Off_PSC_syn,X_ROn_PSC_syn,ROn_X_PSC3_syn,ROn_TD_PSC_syn,ROff_TD_PSC_syn,X_TD_PSC_syn,C_ROn_PSC3_syn,On_R,On_tau,On_Imask,Off_R,Off_tau,Off_Imask,ROn_R,ROn_tau,ROn_Imask,ROff_R,ROff_tau,ROff_Imask,SOnOff_R,SOnOff_tau,SOnOff_Imask,TD_R,TD_tau,TD_Imask,TD_Icur,X_R,X_tau,X_Imask,C_R,C_tau,C_Imask,On_On_IC_netcon,On_On_IC_input,Off_Off_IC_netcon,Off_Off_IC_input,ROn_On_PSC3_netcon,ROn_On_PSC3_scale,SOnOff_On_PSC_netcon,SOnOff_On_PSC_scale,ROn_SOnOff_PSC3_netcon,ROn_SOnOff_PSC3_scale,ROff_SOnOff_PSC_netcon,ROff_SOnOff_PSC_scale,ROff_Off_PSC_netcon,ROff_Off_PSC_scale,SOnOff_Off_PSC_netcon,SOnOff_Off_PSC_scale,ROn_ROn_iNoise_V3_netcon,ROn_ROn_iNoise_V3_token,ROn_ROn_iNoise_V3_scale,X_ROn_PSC_netcon,X_ROn_PSC_scale,ROn_X_PSC3_netcon,ROn_X_PSC3_scale,ROn_TD_PSC_netcon,ROn_TD_PSC_scale,ROff_TD_PSC_netcon,ROff_TD_PSC_scale,X_TD_PSC_netcon,X_TD_PSC_scale,C_ROn_PSC3_netcon,C_ROn_PSC3_scale]=solve_ode(tspan, downsample_factor, random_seed, solver, disk_flag, dt, datafile, mex_flag, verbose_flag, On_C, On_g_L, On_E_L, On_noise, On_t_ref, On_E_k, On_tau_ad, On_g_inc, On_Itonic, On_V_thresh, On_V_reset, On_Npop, Off_C, Off_g_L, Off_E_L, Off_noise, Off_t_ref, Off_E_k, Off_tau_ad, Off_g_inc, Off_Itonic, Off_V_thresh, Off_V_reset, Off_Npop, ROn_C, ROn_g_L, ROn_E_L, ROn_noise, ROn_t_ref, ROn_E_k, ROn_tau_ad, ROn_g_inc, ROn_Itonic, ROn_V_thresh, ROn_V_reset, ROn_Npop, ROff_C, ROff_g_L, ROff_E_L, ROff_noise, ROff_t_ref, ROff_E_k, ROff_tau_ad, ROff_g_inc, ROff_Itonic, ROff_V_thresh, ROff_V_reset, ROff_Npop, SOnOff_C, SOnOff_g_L, SOnOff_E_L, SOnOff_noise, SOnOff_t_ref, SOnOff_E_k, SOnOff_tau_ad, SOnOff_g_inc, SOnOff_Itonic, SOnOff_V_thresh, SOnOff_V_reset, SOnOff_Npop, TD_C, TD_g_L, TD_E_L, TD_noise, TD_t_ref, TD_E_k, TD_tau_ad, TD_g_inc, TD_V_thresh, TD_V_reset, TD_Itonic, TD_numLocs, TD_Npop, X_C, X_g_L, X_E_L, X_noise, X_t_ref, X_E_k, X_tau_ad, X_g_inc, X_Itonic, X_V_thresh, X_V_reset, X_Npop, C_C, C_g_L, C_E_L, C_noise, C_t_ref, C_E_k, C_tau_ad, C_g_inc, C_Itonic, C_V_thresh, C_V_reset, C_Npop, On_On_IC_trial, On_On_IC_locNum, On_On_IC_label, On_On_IC_t_ref, On_On_IC_t_ref_rel, On_On_IC_rec, On_On_IC_g_postIC, On_On_IC_E_exc, Off_Off_IC_trial, Off_Off_IC_locNum, Off_Off_IC_label, Off_Off_IC_t_ref, Off_Off_IC_t_ref_rel, Off_Off_IC_rec, Off_Off_IC_g_postIC, Off_Off_IC_E_exc, ROn_On_PSC3_ESYN, ROn_On_PSC3_tauD, ROn_On_PSC3_tauR, ROn_On_PSC3_delay, ROn_On_PSC3_gSYN, ROn_On_PSC3_fF, ROn_On_PSC3_fP, ROn_On_PSC3_tauF, ROn_On_PSC3_tauP, ROn_On_PSC3_maxF, SOnOff_On_PSC_ESYN, SOnOff_On_PSC_tauD, SOnOff_On_PSC_tauR, SOnOff_On_PSC_delay, SOnOff_On_PSC_gSYN, SOnOff_On_PSC_fF, SOnOff_On_PSC_fP, SOnOff_On_PSC_tauF, SOnOff_On_PSC_tauP, SOnOff_On_PSC_maxF, ROn_SOnOff_PSC3_ESYN, ROn_SOnOff_PSC3_tauD, ROn_SOnOff_PSC3_tauR, ROn_SOnOff_PSC3_delay, ROn_SOnOff_PSC3_gSYN, ROn_SOnOff_PSC3_fF, ROn_SOnOff_PSC3_fP, ROn_SOnOff_PSC3_tauF, ROn_SOnOff_PSC3_tauP, ROn_SOnOff_PSC3_maxF, ROff_SOnOff_PSC_ESYN, ROff_SOnOff_PSC_tauD, ROff_SOnOff_PSC_tauR, ROff_SOnOff_PSC_delay, ROff_SOnOff_PSC_gSYN, ROff_SOnOff_PSC_fF, ROff_SOnOff_PSC_fP, ROff_SOnOff_PSC_tauF, ROff_SOnOff_PSC_tauP, ROff_SOnOff_PSC_maxF, ROff_Off_PSC_ESYN, ROff_Off_PSC_tauD, ROff_Off_PSC_tauR, ROff_Off_PSC_delay, ROff_Off_PSC_gSYN, ROff_Off_PSC_fF, ROff_Off_PSC_fP, ROff_Off_PSC_tauF, ROff_Off_PSC_tauP, ROff_Off_PSC_maxF, SOnOff_Off_PSC_ESYN, SOnOff_Off_PSC_tauD, SOnOff_Off_PSC_tauR, SOnOff_Off_PSC_delay, SOnOff_Off_PSC_gSYN, SOnOff_Off_PSC_fF, SOnOff_Off_PSC_fP, SOnOff_Off_PSC_tauF, SOnOff_Off_PSC_tauP, SOnOff_Off_PSC_maxF, ROn_ROn_iNoise_V3_FR, ROn_ROn_iNoise_V3_sigma, ROn_ROn_iNoise_V3_dt, ROn_ROn_iNoise_V3_nSYN, ROn_ROn_iNoise_V3_simlen, ROn_ROn_iNoise_V3_tauD_N, ROn_ROn_iNoise_V3_tauR_N, ROn_ROn_iNoise_V3_E_exc, X_ROn_PSC_ESYN, X_ROn_PSC_tauD, X_ROn_PSC_tauR, X_ROn_PSC_delay, X_ROn_PSC_gSYN, X_ROn_PSC_fF, X_ROn_PSC_fP, X_ROn_PSC_tauF, X_ROn_PSC_tauP, X_ROn_PSC_maxF, ROn_X_PSC3_ESYN, ROn_X_PSC3_tauD, ROn_X_PSC3_tauR, ROn_X_PSC3_delay, ROn_X_PSC3_gSYN, ROn_X_PSC3_fF, ROn_X_PSC3_fP, ROn_X_PSC3_tauF, ROn_X_PSC3_tauP, ROn_X_PSC3_maxF, ROn_TD_PSC_ESYN, ROn_TD_PSC_tauD, ROn_TD_PSC_tauR, ROn_TD_PSC_delay, ROn_TD_PSC_gSYN, ROn_TD_PSC_fF, ROn_TD_PSC_fP, ROn_TD_PSC_tauF, ROn_TD_PSC_tauP, ROn_TD_PSC_maxF, ROff_TD_PSC_ESYN, ROff_TD_PSC_tauD, ROff_TD_PSC_tauR, ROff_TD_PSC_delay, ROff_TD_PSC_gSYN, ROff_TD_PSC_fF, ROff_TD_PSC_fP, ROff_TD_PSC_tauF, ROff_TD_PSC_tauP, ROff_TD_PSC_maxF, X_TD_PSC_ESYN, X_TD_PSC_tauD, X_TD_PSC_tauR, X_TD_PSC_delay, X_TD_PSC_gSYN, X_TD_PSC_fF, X_TD_PSC_fP, X_TD_PSC_tauF, X_TD_PSC_tauP, X_TD_PSC_maxF, C_ROn_PSC3_ESYN, C_ROn_PSC3_tauD, C_ROn_PSC3_tauR, C_ROn_PSC3_delay, C_ROn_PSC3_gSYN, C_ROn_PSC3_fF, C_ROn_PSC3_fP, C_ROn_PSC3_tauF, C_ROn_PSC3_tauP, C_ROn_PSC3_maxF, ROn_X_PSC3_netcon, ROn_SOnOff_PSC3_netcon, C_ROn_PSC3_netcon)

% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
storeFields = fieldnames(p);
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);

% seed the random number generator
rng(p.random_seed);

% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
On_R = 1/On_g_L;
On_tau = On_C*On_R;
On_Imask =  ones(1,On_Npop);
Off_R = 1/Off_g_L;
Off_tau = Off_C*Off_R;
Off_Imask =  ones(1,Off_Npop);
ROn_R = 1/ROn_g_L;
ROn_tau = ROn_C*ROn_R;
ROn_Imask =  ones(1,ROn_Npop);
ROff_R = 1/ROff_g_L;
ROff_tau = ROff_C*ROff_R;
ROff_Imask =  ones(1,ROff_Npop);
SOnOff_R = 1/SOnOff_g_L;
SOnOff_tau = SOnOff_C*SOnOff_R;
SOnOff_Imask =  ones(1,SOnOff_Npop);
TD_R = 1/TD_g_L;
TD_tau = TD_C*TD_R;
TD_Imask =  ones(1,TD_Npop);
TD_Icur =  TD_Itonic*buildTonicCurrent(T,TD_Npop,dt,TD_numLocs);
X_R = 1/X_g_L;
X_tau = X_C*X_R;
X_Imask =  ones(1,X_Npop);
C_R = 1/C_g_L;
C_tau = C_C*C_R;
C_Imask =  ones(1,C_Npop);
On_On_IC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
On_On_IC_input =  genPoissonInputs(On_On_IC_trial,On_On_IC_locNum,On_On_IC_label,On_On_IC_t_ref,On_On_IC_t_ref_rel,On_On_IC_rec);
Off_Off_IC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
Off_Off_IC_input =  genPoissonInputs(Off_Off_IC_trial,Off_Off_IC_locNum,Off_Off_IC_label,Off_Off_IC_t_ref,Off_Off_IC_t_ref_rel,Off_Off_IC_rec);
ROn_On_PSC3_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
ROn_On_PSC3_scale = (ROn_On_PSC3_tauD/ROn_On_PSC3_tauR)^(ROn_On_PSC3_tauR/(ROn_On_PSC3_tauD-ROn_On_PSC3_tauR));
SOnOff_On_PSC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
SOnOff_On_PSC_scale = (SOnOff_On_PSC_tauD/SOnOff_On_PSC_tauR)^(SOnOff_On_PSC_tauR/(SOnOff_On_PSC_tauD-SOnOff_On_PSC_tauR));
ROn_SOnOff_PSC3_netcon = ROn_SOnOff_PSC3_netcon;
ROn_SOnOff_PSC3_scale = (ROn_SOnOff_PSC3_tauD/ROn_SOnOff_PSC3_tauR)^(ROn_SOnOff_PSC3_tauR/(ROn_SOnOff_PSC3_tauD-ROn_SOnOff_PSC3_tauR));
ROff_SOnOff_PSC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
ROff_SOnOff_PSC_scale = (ROff_SOnOff_PSC_tauD/ROff_SOnOff_PSC_tauR)^(ROff_SOnOff_PSC_tauR/(ROff_SOnOff_PSC_tauD-ROff_SOnOff_PSC_tauR));
ROff_Off_PSC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
ROff_Off_PSC_scale = (ROff_Off_PSC_tauD/ROff_Off_PSC_tauR)^(ROff_Off_PSC_tauR/(ROff_Off_PSC_tauD-ROff_Off_PSC_tauR));
SOnOff_Off_PSC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
SOnOff_Off_PSC_scale = (SOnOff_Off_PSC_tauD/SOnOff_Off_PSC_tauR)^(SOnOff_Off_PSC_tauR/(SOnOff_Off_PSC_tauD-SOnOff_Off_PSC_tauR));
ROn_ROn_iNoise_V3_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
ROn_ROn_iNoise_V3_token = genPoissonTimes(ROn_Npop,ROn_ROn_iNoise_V3_dt,ROn_ROn_iNoise_V3_FR,ROn_ROn_iNoise_V3_sigma,ROn_ROn_iNoise_V3_simlen);
ROn_ROn_iNoise_V3_scale =  (ROn_ROn_iNoise_V3_tauD_N/ROn_ROn_iNoise_V3_tauR_N)^(ROn_ROn_iNoise_V3_tauR_N/(ROn_ROn_iNoise_V3_tauD_N-ROn_ROn_iNoise_V3_tauR_N));
X_ROn_PSC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
X_ROn_PSC_scale = (X_ROn_PSC_tauD/X_ROn_PSC_tauR)^(X_ROn_PSC_tauR/(X_ROn_PSC_tauD-X_ROn_PSC_tauR));
ROn_X_PSC3_netcon = ROn_X_PSC3_netcon;
ROn_X_PSC3_scale = (ROn_X_PSC3_tauD/ROn_X_PSC3_tauR)^(ROn_X_PSC3_tauR/(ROn_X_PSC3_tauD-ROn_X_PSC3_tauR));
ROn_TD_PSC_netcon = [+0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00];
ROn_TD_PSC_scale = (ROn_TD_PSC_tauD/ROn_TD_PSC_tauR)^(ROn_TD_PSC_tauR/(ROn_TD_PSC_tauD-ROn_TD_PSC_tauR));
ROff_TD_PSC_netcon = [+0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00];
ROff_TD_PSC_scale = (ROff_TD_PSC_tauD/ROff_TD_PSC_tauR)^(ROff_TD_PSC_tauR/(ROff_TD_PSC_tauD-ROff_TD_PSC_tauR));
X_TD_PSC_netcon = [+0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00];
X_TD_PSC_scale = (X_TD_PSC_tauD/X_TD_PSC_tauR)^(X_TD_PSC_tauR/(X_TD_PSC_tauD-X_TD_PSC_tauR));
C_ROn_PSC3_netcon = C_ROn_PSC3_netcon;
C_ROn_PSC3_scale = (C_ROn_PSC3_tauD/C_ROn_PSC3_tauR)^(C_ROn_PSC3_tauR/(C_ROn_PSC3_tauD-C_ROn_PSC3_tauR));

% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
t=0; k=1;

% STATE_VARIABLES:
On_V = zeros(2,On_Npop);
On_V(1,:) =  On_E_L*ones(1,On_Npop);
On_g_ad = zeros(2,On_Npop);
On_g_ad(1,:) =  zeros(1,On_Npop);
Off_V = zeros(2,Off_Npop);
Off_V(1,:) =  Off_E_L*ones(1,Off_Npop);
Off_g_ad = zeros(2,Off_Npop);
Off_g_ad(1,:) =  zeros(1,Off_Npop);
ROn_V = zeros(2,ROn_Npop);
ROn_V(1,:) =  ROn_E_L*ones(1,ROn_Npop);
ROn_g_ad = zeros(2,ROn_Npop);
ROn_g_ad(1,:) =  zeros(1,ROn_Npop);
ROff_V = zeros(2,ROff_Npop);
ROff_V(1,:) =  ROff_E_L*ones(1,ROff_Npop);
ROff_g_ad = zeros(2,ROff_Npop);
ROff_g_ad(1,:) =  zeros(1,ROff_Npop);
SOnOff_V = zeros(2,SOnOff_Npop);
SOnOff_V(1,:) =  SOnOff_E_L*ones(1,SOnOff_Npop);
SOnOff_g_ad = zeros(2,SOnOff_Npop);
SOnOff_g_ad(1,:) =  zeros(1,SOnOff_Npop);
TD_V = zeros(2,TD_Npop);
TD_V(1,:) =  TD_E_L*ones(1,TD_Npop);
TD_g_ad = zeros(2,TD_Npop);
TD_g_ad(1,:) =  zeros(1,TD_Npop);
X_V = zeros(2,X_Npop);
X_V(1,:) =  X_E_L*ones(1,X_Npop);
X_g_ad = zeros(2,X_Npop);
X_g_ad(1,:) =  zeros(1,X_Npop);
C_V = zeros(2,C_Npop);
C_V(1,:) =  C_E_L*ones(1,C_Npop);
C_g_ad = zeros(2,C_Npop);
C_g_ad(1,:) =  zeros(1,C_Npop);
ROn_On_PSC3_s = zeros(2,On_Npop);
ROn_On_PSC3_s(1,:) =  zeros(1,On_Npop);
ROn_On_PSC3_x = zeros(2,On_Npop);
ROn_On_PSC3_x(1,:) =  zeros(1,On_Npop);
ROn_On_PSC3_F = zeros(2,On_Npop);
ROn_On_PSC3_F(1,:) =  ones(1,On_Npop);
ROn_On_PSC3_P = zeros(2,On_Npop);
ROn_On_PSC3_P(1,:) =  ones(1,On_Npop);
ROn_On_PSC3_q = zeros(2,On_Npop);
ROn_On_PSC3_q(1,:) =  ones(1,On_Npop);
SOnOff_On_PSC_s = zeros(2,On_Npop);
SOnOff_On_PSC_s(1,:) =  zeros(1,On_Npop);
SOnOff_On_PSC_x = zeros(2,On_Npop);
SOnOff_On_PSC_x(1,:) =  zeros(1,On_Npop);
SOnOff_On_PSC_F = zeros(2,On_Npop);
SOnOff_On_PSC_F(1,:) =  ones(1,On_Npop);
SOnOff_On_PSC_P = zeros(2,On_Npop);
SOnOff_On_PSC_P(1,:) =  ones(1,On_Npop);
SOnOff_On_PSC_q = zeros(2,On_Npop);
SOnOff_On_PSC_q(1,:) =  ones(1,On_Npop);
ROn_SOnOff_PSC3_s = zeros(2,SOnOff_Npop);
ROn_SOnOff_PSC3_s(1,:) =  zeros(1,SOnOff_Npop);
ROn_SOnOff_PSC3_x = zeros(2,SOnOff_Npop);
ROn_SOnOff_PSC3_x(1,:) =  zeros(1,SOnOff_Npop);
ROn_SOnOff_PSC3_F = zeros(2,SOnOff_Npop);
ROn_SOnOff_PSC3_F(1,:) =  ones(1,SOnOff_Npop);
ROn_SOnOff_PSC3_P = zeros(2,SOnOff_Npop);
ROn_SOnOff_PSC3_P(1,:) =  ones(1,SOnOff_Npop);
ROn_SOnOff_PSC3_q = zeros(2,SOnOff_Npop);
ROn_SOnOff_PSC3_q(1,:) =  ones(1,SOnOff_Npop);
ROff_SOnOff_PSC_s = zeros(2,SOnOff_Npop);
ROff_SOnOff_PSC_s(1,:) =  zeros(1,SOnOff_Npop);
ROff_SOnOff_PSC_x = zeros(2,SOnOff_Npop);
ROff_SOnOff_PSC_x(1,:) =  zeros(1,SOnOff_Npop);
ROff_SOnOff_PSC_F = zeros(2,SOnOff_Npop);
ROff_SOnOff_PSC_F(1,:) =  ones(1,SOnOff_Npop);
ROff_SOnOff_PSC_P = zeros(2,SOnOff_Npop);
ROff_SOnOff_PSC_P(1,:) =  ones(1,SOnOff_Npop);
ROff_SOnOff_PSC_q = zeros(2,SOnOff_Npop);
ROff_SOnOff_PSC_q(1,:) =  ones(1,SOnOff_Npop);
ROff_Off_PSC_s = zeros(2,Off_Npop);
ROff_Off_PSC_s(1,:) =  zeros(1,Off_Npop);
ROff_Off_PSC_x = zeros(2,Off_Npop);
ROff_Off_PSC_x(1,:) =  zeros(1,Off_Npop);
ROff_Off_PSC_F = zeros(2,Off_Npop);
ROff_Off_PSC_F(1,:) =  ones(1,Off_Npop);
ROff_Off_PSC_P = zeros(2,Off_Npop);
ROff_Off_PSC_P(1,:) =  ones(1,Off_Npop);
ROff_Off_PSC_q = zeros(2,Off_Npop);
ROff_Off_PSC_q(1,:) =  ones(1,Off_Npop);
SOnOff_Off_PSC_s = zeros(2,Off_Npop);
SOnOff_Off_PSC_s(1,:) =  zeros(1,Off_Npop);
SOnOff_Off_PSC_x = zeros(2,Off_Npop);
SOnOff_Off_PSC_x(1,:) =  zeros(1,Off_Npop);
SOnOff_Off_PSC_F = zeros(2,Off_Npop);
SOnOff_Off_PSC_F(1,:) =  ones(1,Off_Npop);
SOnOff_Off_PSC_P = zeros(2,Off_Npop);
SOnOff_Off_PSC_P(1,:) =  ones(1,Off_Npop);
SOnOff_Off_PSC_q = zeros(2,Off_Npop);
SOnOff_Off_PSC_q(1,:) =  ones(1,Off_Npop);
ROn_ROn_iNoise_V3_sn = zeros(2,ROn_Npop);
ROn_ROn_iNoise_V3_sn(1,:) =  0 * ones(1,ROn_Npop);
ROn_ROn_iNoise_V3_xn = zeros(2,ROn_Npop);
ROn_ROn_iNoise_V3_xn(1,:) =  0 * ones(1,ROn_Npop);
X_ROn_PSC_s = zeros(2,ROn_Npop);
X_ROn_PSC_s(1,:) =  zeros(1,ROn_Npop);
X_ROn_PSC_x = zeros(2,ROn_Npop);
X_ROn_PSC_x(1,:) =  zeros(1,ROn_Npop);
X_ROn_PSC_F = zeros(2,ROn_Npop);
X_ROn_PSC_F(1,:) =  ones(1,ROn_Npop);
X_ROn_PSC_P = zeros(2,ROn_Npop);
X_ROn_PSC_P(1,:) =  ones(1,ROn_Npop);
X_ROn_PSC_q = zeros(2,ROn_Npop);
X_ROn_PSC_q(1,:) =  ones(1,ROn_Npop);
ROn_X_PSC3_s = zeros(2,X_Npop);
ROn_X_PSC3_s(1,:) =  zeros(1,X_Npop);
ROn_X_PSC3_x = zeros(2,X_Npop);
ROn_X_PSC3_x(1,:) =  zeros(1,X_Npop);
ROn_X_PSC3_F = zeros(2,X_Npop);
ROn_X_PSC3_F(1,:) =  ones(1,X_Npop);
ROn_X_PSC3_P = zeros(2,X_Npop);
ROn_X_PSC3_P(1,:) =  ones(1,X_Npop);
ROn_X_PSC3_q = zeros(2,X_Npop);
ROn_X_PSC3_q(1,:) =  ones(1,X_Npop);
ROn_TD_PSC_s = zeros(2,TD_Npop);
ROn_TD_PSC_s(1,:) =  zeros(1,TD_Npop);
ROn_TD_PSC_x = zeros(2,TD_Npop);
ROn_TD_PSC_x(1,:) =  zeros(1,TD_Npop);
ROn_TD_PSC_F = zeros(2,TD_Npop);
ROn_TD_PSC_F(1,:) =  ones(1,TD_Npop);
ROn_TD_PSC_P = zeros(2,TD_Npop);
ROn_TD_PSC_P(1,:) =  ones(1,TD_Npop);
ROn_TD_PSC_q = zeros(2,TD_Npop);
ROn_TD_PSC_q(1,:) =  ones(1,TD_Npop);
ROff_TD_PSC_s = zeros(2,TD_Npop);
ROff_TD_PSC_s(1,:) =  zeros(1,TD_Npop);
ROff_TD_PSC_x = zeros(2,TD_Npop);
ROff_TD_PSC_x(1,:) =  zeros(1,TD_Npop);
ROff_TD_PSC_F = zeros(2,TD_Npop);
ROff_TD_PSC_F(1,:) =  ones(1,TD_Npop);
ROff_TD_PSC_P = zeros(2,TD_Npop);
ROff_TD_PSC_P(1,:) =  ones(1,TD_Npop);
ROff_TD_PSC_q = zeros(2,TD_Npop);
ROff_TD_PSC_q(1,:) =  ones(1,TD_Npop);
X_TD_PSC_s = zeros(2,TD_Npop);
X_TD_PSC_s(1,:) =  zeros(1,TD_Npop);
X_TD_PSC_x = zeros(2,TD_Npop);
X_TD_PSC_x(1,:) =  zeros(1,TD_Npop);
X_TD_PSC_F = zeros(2,TD_Npop);
X_TD_PSC_F(1,:) =  ones(1,TD_Npop);
X_TD_PSC_P = zeros(2,TD_Npop);
X_TD_PSC_P(1,:) =  ones(1,TD_Npop);
X_TD_PSC_q = zeros(2,TD_Npop);
X_TD_PSC_q(1,:) =  ones(1,TD_Npop);
C_ROn_PSC3_s = zeros(2,ROn_Npop);
C_ROn_PSC3_s(1,:) =  zeros(1,ROn_Npop);
C_ROn_PSC3_x = zeros(2,ROn_Npop);
C_ROn_PSC3_x(1,:) =  zeros(1,ROn_Npop);
C_ROn_PSC3_F = zeros(2,ROn_Npop);
C_ROn_PSC3_F(1,:) =  ones(1,ROn_Npop);
C_ROn_PSC3_P = zeros(2,ROn_Npop);
C_ROn_PSC3_P(1,:) =  ones(1,ROn_Npop);
C_ROn_PSC3_q = zeros(2,ROn_Npop);
C_ROn_PSC3_q(1,:) =  ones(1,ROn_Npop);

% MONITORS:
On_tspike = -1e32*ones(5,On_Npop);
On_buffer_index = ones(1,On_Npop);
On_V_spikes = zeros(nsamp,On_Npop);
Off_tspike = -1e32*ones(5,Off_Npop);
Off_buffer_index = ones(1,Off_Npop);
Off_V_spikes = zeros(nsamp,Off_Npop);
ROn_tspike = -1e32*ones(5,ROn_Npop);
ROn_buffer_index = ones(1,ROn_Npop);
ROn_V_spikes = zeros(nsamp,ROn_Npop);
ROff_tspike = -1e32*ones(5,ROff_Npop);
ROff_buffer_index = ones(1,ROff_Npop);
ROff_V_spikes = zeros(nsamp,ROff_Npop);
SOnOff_tspike = -1e32*ones(5,SOnOff_Npop);
SOnOff_buffer_index = ones(1,SOnOff_Npop);
SOnOff_V_spikes = zeros(nsamp,SOnOff_Npop);
TD_tspike = -1e32*ones(5,TD_Npop);
TD_buffer_index = ones(1,TD_Npop);
TD_V_spikes = zeros(nsamp,TD_Npop);
X_tspike = -1e32*ones(5,X_Npop);
X_buffer_index = ones(1,X_Npop);
X_V_spikes = zeros(nsamp,X_Npop);
C_tspike = -1e32*ones(5,C_Npop);
C_buffer_index = ones(1,C_Npop);
C_V_spikes = zeros(nsamp,C_Npop);
On_On_IC_iIC = zeros(nsamp,On_Npop);
On_On_IC_iIC(1,:)=On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(1,:)-On_On_IC_E_exc);
Off_Off_IC_iIC = zeros(nsamp,Off_Npop);
Off_Off_IC_iIC(1,:)=Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(1,:)-Off_Off_IC_E_exc);
ROn_On_PSC3_syn = zeros(nsamp,ROn_Npop);
ROn_On_PSC3_syn(1,:)=((ROn_On_PSC3_s(1,:)*(ROn_On_PSC3_netcon.*ROn_On_PSC3_gSYN)).*(ROn_V(1,:)-ROn_On_PSC3_ESYN));
SOnOff_On_PSC_syn = zeros(nsamp,SOnOff_Npop);
SOnOff_On_PSC_syn(1,:)=SOnOff_On_PSC_gSYN.*(SOnOff_On_PSC_s(1,:)*SOnOff_On_PSC_netcon).*(SOnOff_V(1,:)-SOnOff_On_PSC_ESYN);
ROn_SOnOff_PSC3_syn = zeros(nsamp,ROn_Npop);
ROn_SOnOff_PSC3_syn(1,:)=((ROn_SOnOff_PSC3_s(1,:)*(ROn_SOnOff_PSC3_netcon.*ROn_SOnOff_PSC3_gSYN)).*(ROn_V(1,:)-ROn_SOnOff_PSC3_ESYN));
ROff_SOnOff_PSC_syn = zeros(nsamp,ROff_Npop);
ROff_SOnOff_PSC_syn(1,:)=ROff_SOnOff_PSC_gSYN.*(ROff_SOnOff_PSC_s(1,:)*ROff_SOnOff_PSC_netcon).*(ROff_V(1,:)-ROff_SOnOff_PSC_ESYN);
ROff_Off_PSC_syn = zeros(nsamp,ROff_Npop);
ROff_Off_PSC_syn(1,:)=ROff_Off_PSC_gSYN.*(ROff_Off_PSC_s(1,:)*ROff_Off_PSC_netcon).*(ROff_V(1,:)-ROff_Off_PSC_ESYN);
SOnOff_Off_PSC_syn = zeros(nsamp,SOnOff_Npop);
SOnOff_Off_PSC_syn(1,:)=SOnOff_Off_PSC_gSYN.*(SOnOff_Off_PSC_s(1,:)*SOnOff_Off_PSC_netcon).*(SOnOff_V(1,:)-SOnOff_Off_PSC_ESYN);
X_ROn_PSC_syn = zeros(nsamp,X_Npop);
X_ROn_PSC_syn(1,:)=X_ROn_PSC_gSYN.*(X_ROn_PSC_s(1,:)*X_ROn_PSC_netcon).*(X_V(1,:)-X_ROn_PSC_ESYN);
ROn_X_PSC3_syn = zeros(nsamp,ROn_Npop);
ROn_X_PSC3_syn(1,:)=((ROn_X_PSC3_s(1,:)*(ROn_X_PSC3_netcon.*ROn_X_PSC3_gSYN)).*(ROn_V(1,:)-ROn_X_PSC3_ESYN));
ROn_TD_PSC_syn = zeros(nsamp,ROn_Npop);
ROn_TD_PSC_syn(1,:)=ROn_TD_PSC_gSYN.*(ROn_TD_PSC_s(1,:)*ROn_TD_PSC_netcon).*(ROn_V(1,:)-ROn_TD_PSC_ESYN);
ROff_TD_PSC_syn = zeros(nsamp,ROff_Npop);
ROff_TD_PSC_syn(1,:)=ROff_TD_PSC_gSYN.*(ROff_TD_PSC_s(1,:)*ROff_TD_PSC_netcon).*(ROff_V(1,:)-ROff_TD_PSC_ESYN);
X_TD_PSC_syn = zeros(nsamp,X_Npop);
X_TD_PSC_syn(1,:)=X_TD_PSC_gSYN.*(X_TD_PSC_s(1,:)*X_TD_PSC_netcon).*(X_V(1,:)-X_TD_PSC_ESYN);
C_ROn_PSC3_syn = zeros(nsamp,C_Npop);
C_ROn_PSC3_syn(1,:)=((C_ROn_PSC3_s(1,:)*(C_ROn_PSC3_netcon.*C_ROn_PSC3_gSYN)).*(C_V(1,:)-C_ROn_PSC3_ESYN));

% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  On_V_k1 = ( (On_E_L-On_V(1,:)) - On_R*On_g_ad(1,:).*(On_V(1,:)-On_E_k) - On_R*((((On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(1,:)-On_On_IC_E_exc))))) + On_R*On_Itonic.*On_Imask + On_R*On_noise.*randn(1,On_Npop) ) / On_tau;
  On_g_ad_k1 = -On_g_ad(1,:) / On_tau_ad;
  Off_V_k1 = ( (Off_E_L-Off_V(1,:)) - Off_R*Off_g_ad(1,:).*(Off_V(1,:)-Off_E_k) - Off_R*((((Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(1,:)-Off_Off_IC_E_exc))))) + Off_R*Off_Itonic.*Off_Imask + Off_R*Off_noise.*randn(1,Off_Npop) ) / Off_tau;
  Off_g_ad_k1 = -Off_g_ad(1,:) / Off_tau_ad;
  ROn_V_k1 = ( (ROn_E_L-ROn_V(1,:)) - ROn_R*ROn_g_ad(1,:).*(ROn_V(1,:)-ROn_E_k) - ROn_R*((((((ROn_On_PSC3_s(1,:)*(ROn_On_PSC3_netcon.*ROn_On_PSC3_gSYN)).*(ROn_V(1,:)-ROn_On_PSC3_ESYN)))))+((((((ROn_SOnOff_PSC3_s(1,:)*(ROn_SOnOff_PSC3_netcon.*ROn_SOnOff_PSC3_gSYN)).*(ROn_V(1,:)-ROn_SOnOff_PSC3_ESYN)))))+((((ROn_ROn_iNoise_V3_nSYN.*(ROn_ROn_iNoise_V3_sn(1,:)*ROn_ROn_iNoise_V3_netcon).*(ROn_V(1,:)-ROn_ROn_iNoise_V3_E_exc))))+((((((ROn_X_PSC3_s(1,:)*(ROn_X_PSC3_netcon.*ROn_X_PSC3_gSYN)).*(ROn_V(1,:)-ROn_X_PSC3_ESYN)))))+((((ROn_TD_PSC_gSYN.*(ROn_TD_PSC_s(1,:)*ROn_TD_PSC_netcon).*(ROn_V(1,:)-ROn_TD_PSC_ESYN))))))))) + ROn_R*ROn_Itonic.*ROn_Imask + ROn_R*ROn_noise.*randn(1,ROn_Npop) ) / ROn_tau;
  ROn_g_ad_k1 = -ROn_g_ad(1,:) / ROn_tau_ad;
  ROff_V_k1 = ( (ROff_E_L-ROff_V(1,:)) - ROff_R*ROff_g_ad(1,:).*(ROff_V(1,:)-ROff_E_k) - ROff_R*((((ROff_SOnOff_PSC_gSYN.*(ROff_SOnOff_PSC_s(1,:)*ROff_SOnOff_PSC_netcon).*(ROff_V(1,:)-ROff_SOnOff_PSC_ESYN))))+((((ROff_Off_PSC_gSYN.*(ROff_Off_PSC_s(1,:)*ROff_Off_PSC_netcon).*(ROff_V(1,:)-ROff_Off_PSC_ESYN))))+((((ROff_TD_PSC_gSYN.*(ROff_TD_PSC_s(1,:)*ROff_TD_PSC_netcon).*(ROff_V(1,:)-ROff_TD_PSC_ESYN))))))) + ROff_R*ROff_Itonic.*ROff_Imask + ROff_R*ROff_noise.*randn(1,ROff_Npop) ) / ROff_tau;
  ROff_g_ad_k1 = -ROff_g_ad(1,:) / ROff_tau_ad;
  SOnOff_V_k1 = ( (SOnOff_E_L-SOnOff_V(1,:)) - SOnOff_R*SOnOff_g_ad(1,:).*(SOnOff_V(1,:)-SOnOff_E_k) - SOnOff_R*((((SOnOff_On_PSC_gSYN.*(SOnOff_On_PSC_s(1,:)*SOnOff_On_PSC_netcon).*(SOnOff_V(1,:)-SOnOff_On_PSC_ESYN))))+((((SOnOff_Off_PSC_gSYN.*(SOnOff_Off_PSC_s(1,:)*SOnOff_Off_PSC_netcon).*(SOnOff_V(1,:)-SOnOff_Off_PSC_ESYN)))))) + SOnOff_R*SOnOff_Itonic.*SOnOff_Imask + SOnOff_R*SOnOff_noise.*randn(1,SOnOff_Npop) ) / SOnOff_tau;
  SOnOff_g_ad_k1 = -SOnOff_g_ad(1,:) / SOnOff_tau_ad;
  TD_V_k1 = ( (TD_E_L-TD_V(1,:)) - TD_R*TD_g_ad(1,:).*(TD_V(1,:)-TD_E_k) + TD_R*TD_Icur(k,:).*TD_Imask + TD_R*TD_noise.*randn(1,TD_Npop) ) / TD_tau;
  TD_g_ad_k1 = -TD_g_ad(1,:) / TD_tau_ad;
  X_V_k1 = ( (X_E_L-X_V(1,:)) - X_R*X_g_ad(1,:).*(X_V(1,:)-X_E_k) - X_R*((((X_ROn_PSC_gSYN.*(X_ROn_PSC_s(1,:)*X_ROn_PSC_netcon).*(X_V(1,:)-X_ROn_PSC_ESYN))))+((((X_TD_PSC_gSYN.*(X_TD_PSC_s(1,:)*X_TD_PSC_netcon).*(X_V(1,:)-X_TD_PSC_ESYN)))))) + X_R*X_Itonic.*X_Imask + X_R*X_noise.*randn(1,X_Npop) ) / X_tau;
  X_g_ad_k1 = -X_g_ad(1,:) / X_tau_ad;
  C_V_k1 = ( (C_E_L-C_V(1,:)) - C_R*C_g_ad(1,:).*(C_V(1,:)-C_E_k) - C_R*((((((C_ROn_PSC3_s(1,:)*(C_ROn_PSC3_netcon.*C_ROn_PSC3_gSYN)).*(C_V(1,:)-C_ROn_PSC3_ESYN)))))) + C_R*C_Itonic.*C_Imask + C_R*C_noise.*randn(1,C_Npop) ) / C_tau;
  C_g_ad_k1 = -C_g_ad(1,:) / C_tau_ad;
  ROn_On_PSC3_s_k1 = ( ROn_On_PSC3_scale * ROn_On_PSC3_x(1,:) - ROn_On_PSC3_s(1,:) )/ROn_On_PSC3_tauR;
  ROn_On_PSC3_x_k1 = -ROn_On_PSC3_x(1,:)/ROn_On_PSC3_tauD;
  ROn_On_PSC3_F_k1 = (1 - ROn_On_PSC3_F(1,:))/ROn_On_PSC3_tauF;
  ROn_On_PSC3_P_k1 = (1 - ROn_On_PSC3_P(1,:))/ROn_On_PSC3_tauP;
  ROn_On_PSC3_q_k1 = 0;
  SOnOff_On_PSC_s_k1 = ( SOnOff_On_PSC_scale * SOnOff_On_PSC_x(1,:) - SOnOff_On_PSC_s(1,:) )/SOnOff_On_PSC_tauR;
  SOnOff_On_PSC_x_k1 = -SOnOff_On_PSC_x(1,:)/SOnOff_On_PSC_tauD;
  SOnOff_On_PSC_F_k1 = (1 - SOnOff_On_PSC_F(1,:))/SOnOff_On_PSC_tauF;
  SOnOff_On_PSC_P_k1 = (1 - SOnOff_On_PSC_P(1,:))/SOnOff_On_PSC_tauP;
  SOnOff_On_PSC_q_k1 = 0;
  ROn_SOnOff_PSC3_s_k1 = ( ROn_SOnOff_PSC3_scale * ROn_SOnOff_PSC3_x(1,:) - ROn_SOnOff_PSC3_s(1,:) )/ROn_SOnOff_PSC3_tauR;
  ROn_SOnOff_PSC3_x_k1 = -ROn_SOnOff_PSC3_x(1,:)/ROn_SOnOff_PSC3_tauD;
  ROn_SOnOff_PSC3_F_k1 = (1 - ROn_SOnOff_PSC3_F(1,:))/ROn_SOnOff_PSC3_tauF;
  ROn_SOnOff_PSC3_P_k1 = (1 - ROn_SOnOff_PSC3_P(1,:))/ROn_SOnOff_PSC3_tauP;
  ROn_SOnOff_PSC3_q_k1 = 0;
  ROff_SOnOff_PSC_s_k1 = ( ROff_SOnOff_PSC_scale * ROff_SOnOff_PSC_x(1,:) - ROff_SOnOff_PSC_s(1,:) )/ROff_SOnOff_PSC_tauR;
  ROff_SOnOff_PSC_x_k1 = -ROff_SOnOff_PSC_x(1,:)/ROff_SOnOff_PSC_tauD;
  ROff_SOnOff_PSC_F_k1 = (1 - ROff_SOnOff_PSC_F(1,:))/ROff_SOnOff_PSC_tauF;
  ROff_SOnOff_PSC_P_k1 = (1 - ROff_SOnOff_PSC_P(1,:))/ROff_SOnOff_PSC_tauP;
  ROff_SOnOff_PSC_q_k1 = 0;
  ROff_Off_PSC_s_k1 = ( ROff_Off_PSC_scale * ROff_Off_PSC_x(1,:) - ROff_Off_PSC_s(1,:) )/ROff_Off_PSC_tauR;
  ROff_Off_PSC_x_k1 = -ROff_Off_PSC_x(1,:)/ROff_Off_PSC_tauD;
  ROff_Off_PSC_F_k1 = (1 - ROff_Off_PSC_F(1,:))/ROff_Off_PSC_tauF;
  ROff_Off_PSC_P_k1 = (1 - ROff_Off_PSC_P(1,:))/ROff_Off_PSC_tauP;
  ROff_Off_PSC_q_k1 = 0;
  SOnOff_Off_PSC_s_k1 = ( SOnOff_Off_PSC_scale * SOnOff_Off_PSC_x(1,:) - SOnOff_Off_PSC_s(1,:) )/SOnOff_Off_PSC_tauR;
  SOnOff_Off_PSC_x_k1 = -SOnOff_Off_PSC_x(1,:)/SOnOff_Off_PSC_tauD;
  SOnOff_Off_PSC_F_k1 = (1 - SOnOff_Off_PSC_F(1,:))/SOnOff_Off_PSC_tauF;
  SOnOff_Off_PSC_P_k1 = (1 - SOnOff_Off_PSC_P(1,:))/SOnOff_Off_PSC_tauP;
  SOnOff_Off_PSC_q_k1 = 0;
  ROn_ROn_iNoise_V3_sn_k1 = ( ROn_ROn_iNoise_V3_scale * ROn_ROn_iNoise_V3_xn(1,:) - ROn_ROn_iNoise_V3_sn(1,:) )/ROn_ROn_iNoise_V3_tauR_N;
  ROn_ROn_iNoise_V3_xn_k1 = -ROn_ROn_iNoise_V3_xn(1,:)/ROn_ROn_iNoise_V3_tauD_N + ROn_ROn_iNoise_V3_token(k,:)/ROn_ROn_iNoise_V3_dt;
  X_ROn_PSC_s_k1 = ( X_ROn_PSC_scale * X_ROn_PSC_x(1,:) - X_ROn_PSC_s(1,:) )/X_ROn_PSC_tauR;
  X_ROn_PSC_x_k1 = -X_ROn_PSC_x(1,:)/X_ROn_PSC_tauD;
  X_ROn_PSC_F_k1 = (1 - X_ROn_PSC_F(1,:))/X_ROn_PSC_tauF;
  X_ROn_PSC_P_k1 = (1 - X_ROn_PSC_P(1,:))/X_ROn_PSC_tauP;
  X_ROn_PSC_q_k1 = 0;
  ROn_X_PSC3_s_k1 = ( ROn_X_PSC3_scale * ROn_X_PSC3_x(1,:) - ROn_X_PSC3_s(1,:) )/ROn_X_PSC3_tauR;
  ROn_X_PSC3_x_k1 = -ROn_X_PSC3_x(1,:)/ROn_X_PSC3_tauD;
  ROn_X_PSC3_F_k1 = (1 - ROn_X_PSC3_F(1,:))/ROn_X_PSC3_tauF;
  ROn_X_PSC3_P_k1 = (1 - ROn_X_PSC3_P(1,:))/ROn_X_PSC3_tauP;
  ROn_X_PSC3_q_k1 = 0;
  ROn_TD_PSC_s_k1 = ( ROn_TD_PSC_scale * ROn_TD_PSC_x(1,:) - ROn_TD_PSC_s(1,:) )/ROn_TD_PSC_tauR;
  ROn_TD_PSC_x_k1 = -ROn_TD_PSC_x(1,:)/ROn_TD_PSC_tauD;
  ROn_TD_PSC_F_k1 = (1 - ROn_TD_PSC_F(1,:))/ROn_TD_PSC_tauF;
  ROn_TD_PSC_P_k1 = (1 - ROn_TD_PSC_P(1,:))/ROn_TD_PSC_tauP;
  ROn_TD_PSC_q_k1 = 0;
  ROff_TD_PSC_s_k1 = ( ROff_TD_PSC_scale * ROff_TD_PSC_x(1,:) - ROff_TD_PSC_s(1,:) )/ROff_TD_PSC_tauR;
  ROff_TD_PSC_x_k1 = -ROff_TD_PSC_x(1,:)/ROff_TD_PSC_tauD;
  ROff_TD_PSC_F_k1 = (1 - ROff_TD_PSC_F(1,:))/ROff_TD_PSC_tauF;
  ROff_TD_PSC_P_k1 = (1 - ROff_TD_PSC_P(1,:))/ROff_TD_PSC_tauP;
  ROff_TD_PSC_q_k1 = 0;
  X_TD_PSC_s_k1 = ( X_TD_PSC_scale * X_TD_PSC_x(1,:) - X_TD_PSC_s(1,:) )/X_TD_PSC_tauR;
  X_TD_PSC_x_k1 = -X_TD_PSC_x(1,:)/X_TD_PSC_tauD;
  X_TD_PSC_F_k1 = (1 - X_TD_PSC_F(1,:))/X_TD_PSC_tauF;
  X_TD_PSC_P_k1 = (1 - X_TD_PSC_P(1,:))/X_TD_PSC_tauP;
  X_TD_PSC_q_k1 = 0;
  C_ROn_PSC3_s_k1 = ( C_ROn_PSC3_scale * C_ROn_PSC3_x(1,:) - C_ROn_PSC3_s(1,:) )/C_ROn_PSC3_tauR;
  C_ROn_PSC3_x_k1 = -C_ROn_PSC3_x(1,:)/C_ROn_PSC3_tauD;
  C_ROn_PSC3_F_k1 = (1 - C_ROn_PSC3_F(1,:))/C_ROn_PSC3_tauF;
  C_ROn_PSC3_P_k1 = (1 - C_ROn_PSC3_P(1,:))/C_ROn_PSC3_tauP;
  C_ROn_PSC3_q_k1 = 0;

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  On_V(2,:) = On_V(1,:)+dt*On_V_k1;
  On_g_ad(2,:) = On_g_ad(1,:)+dt*On_g_ad_k1;
  Off_V(2,:) = Off_V(1,:)+dt*Off_V_k1;
  Off_g_ad(2,:) = Off_g_ad(1,:)+dt*Off_g_ad_k1;
  ROn_V(2,:) = ROn_V(1,:)+dt*ROn_V_k1;
  ROn_g_ad(2,:) = ROn_g_ad(1,:)+dt*ROn_g_ad_k1;
  ROff_V(2,:) = ROff_V(1,:)+dt*ROff_V_k1;
  ROff_g_ad(2,:) = ROff_g_ad(1,:)+dt*ROff_g_ad_k1;
  SOnOff_V(2,:) = SOnOff_V(1,:)+dt*SOnOff_V_k1;
  SOnOff_g_ad(2,:) = SOnOff_g_ad(1,:)+dt*SOnOff_g_ad_k1;
  TD_V(2,:) = TD_V(1,:)+dt*TD_V_k1;
  TD_g_ad(2,:) = TD_g_ad(1,:)+dt*TD_g_ad_k1;
  X_V(2,:) = X_V(1,:)+dt*X_V_k1;
  X_g_ad(2,:) = X_g_ad(1,:)+dt*X_g_ad_k1;
  C_V(2,:) = C_V(1,:)+dt*C_V_k1;
  C_g_ad(2,:) = C_g_ad(1,:)+dt*C_g_ad_k1;
  ROn_On_PSC3_s(2,:) = ROn_On_PSC3_s(1,:)+dt*ROn_On_PSC3_s_k1;
  ROn_On_PSC3_x(2,:) = ROn_On_PSC3_x(1,:)+dt*ROn_On_PSC3_x_k1;
  ROn_On_PSC3_F(2,:) = ROn_On_PSC3_F(1,:)+dt*ROn_On_PSC3_F_k1;
  ROn_On_PSC3_P(2,:) = ROn_On_PSC3_P(1,:)+dt*ROn_On_PSC3_P_k1;
  ROn_On_PSC3_q(2,:) = ROn_On_PSC3_q(1,:)+dt*ROn_On_PSC3_q_k1;
  SOnOff_On_PSC_s(2,:) = SOnOff_On_PSC_s(1,:)+dt*SOnOff_On_PSC_s_k1;
  SOnOff_On_PSC_x(2,:) = SOnOff_On_PSC_x(1,:)+dt*SOnOff_On_PSC_x_k1;
  SOnOff_On_PSC_F(2,:) = SOnOff_On_PSC_F(1,:)+dt*SOnOff_On_PSC_F_k1;
  SOnOff_On_PSC_P(2,:) = SOnOff_On_PSC_P(1,:)+dt*SOnOff_On_PSC_P_k1;
  SOnOff_On_PSC_q(2,:) = SOnOff_On_PSC_q(1,:)+dt*SOnOff_On_PSC_q_k1;
  ROn_SOnOff_PSC3_s(2,:) = ROn_SOnOff_PSC3_s(1,:)+dt*ROn_SOnOff_PSC3_s_k1;
  ROn_SOnOff_PSC3_x(2,:) = ROn_SOnOff_PSC3_x(1,:)+dt*ROn_SOnOff_PSC3_x_k1;
  ROn_SOnOff_PSC3_F(2,:) = ROn_SOnOff_PSC3_F(1,:)+dt*ROn_SOnOff_PSC3_F_k1;
  ROn_SOnOff_PSC3_P(2,:) = ROn_SOnOff_PSC3_P(1,:)+dt*ROn_SOnOff_PSC3_P_k1;
  ROn_SOnOff_PSC3_q(2,:) = ROn_SOnOff_PSC3_q(1,:)+dt*ROn_SOnOff_PSC3_q_k1;
  ROff_SOnOff_PSC_s(2,:) = ROff_SOnOff_PSC_s(1,:)+dt*ROff_SOnOff_PSC_s_k1;
  ROff_SOnOff_PSC_x(2,:) = ROff_SOnOff_PSC_x(1,:)+dt*ROff_SOnOff_PSC_x_k1;
  ROff_SOnOff_PSC_F(2,:) = ROff_SOnOff_PSC_F(1,:)+dt*ROff_SOnOff_PSC_F_k1;
  ROff_SOnOff_PSC_P(2,:) = ROff_SOnOff_PSC_P(1,:)+dt*ROff_SOnOff_PSC_P_k1;
  ROff_SOnOff_PSC_q(2,:) = ROff_SOnOff_PSC_q(1,:)+dt*ROff_SOnOff_PSC_q_k1;
  ROff_Off_PSC_s(2,:) = ROff_Off_PSC_s(1,:)+dt*ROff_Off_PSC_s_k1;
  ROff_Off_PSC_x(2,:) = ROff_Off_PSC_x(1,:)+dt*ROff_Off_PSC_x_k1;
  ROff_Off_PSC_F(2,:) = ROff_Off_PSC_F(1,:)+dt*ROff_Off_PSC_F_k1;
  ROff_Off_PSC_P(2,:) = ROff_Off_PSC_P(1,:)+dt*ROff_Off_PSC_P_k1;
  ROff_Off_PSC_q(2,:) = ROff_Off_PSC_q(1,:)+dt*ROff_Off_PSC_q_k1;
  SOnOff_Off_PSC_s(2,:) = SOnOff_Off_PSC_s(1,:)+dt*SOnOff_Off_PSC_s_k1;
  SOnOff_Off_PSC_x(2,:) = SOnOff_Off_PSC_x(1,:)+dt*SOnOff_Off_PSC_x_k1;
  SOnOff_Off_PSC_F(2,:) = SOnOff_Off_PSC_F(1,:)+dt*SOnOff_Off_PSC_F_k1;
  SOnOff_Off_PSC_P(2,:) = SOnOff_Off_PSC_P(1,:)+dt*SOnOff_Off_PSC_P_k1;
  SOnOff_Off_PSC_q(2,:) = SOnOff_Off_PSC_q(1,:)+dt*SOnOff_Off_PSC_q_k1;
  ROn_ROn_iNoise_V3_sn(2,:) = ROn_ROn_iNoise_V3_sn(1,:)+dt*ROn_ROn_iNoise_V3_sn_k1;
  ROn_ROn_iNoise_V3_xn(2,:) = ROn_ROn_iNoise_V3_xn(1,:)+dt*ROn_ROn_iNoise_V3_xn_k1;
  X_ROn_PSC_s(2,:) = X_ROn_PSC_s(1,:)+dt*X_ROn_PSC_s_k1;
  X_ROn_PSC_x(2,:) = X_ROn_PSC_x(1,:)+dt*X_ROn_PSC_x_k1;
  X_ROn_PSC_F(2,:) = X_ROn_PSC_F(1,:)+dt*X_ROn_PSC_F_k1;
  X_ROn_PSC_P(2,:) = X_ROn_PSC_P(1,:)+dt*X_ROn_PSC_P_k1;
  X_ROn_PSC_q(2,:) = X_ROn_PSC_q(1,:)+dt*X_ROn_PSC_q_k1;
  ROn_X_PSC3_s(2,:) = ROn_X_PSC3_s(1,:)+dt*ROn_X_PSC3_s_k1;
  ROn_X_PSC3_x(2,:) = ROn_X_PSC3_x(1,:)+dt*ROn_X_PSC3_x_k1;
  ROn_X_PSC3_F(2,:) = ROn_X_PSC3_F(1,:)+dt*ROn_X_PSC3_F_k1;
  ROn_X_PSC3_P(2,:) = ROn_X_PSC3_P(1,:)+dt*ROn_X_PSC3_P_k1;
  ROn_X_PSC3_q(2,:) = ROn_X_PSC3_q(1,:)+dt*ROn_X_PSC3_q_k1;
  ROn_TD_PSC_s(2,:) = ROn_TD_PSC_s(1,:)+dt*ROn_TD_PSC_s_k1;
  ROn_TD_PSC_x(2,:) = ROn_TD_PSC_x(1,:)+dt*ROn_TD_PSC_x_k1;
  ROn_TD_PSC_F(2,:) = ROn_TD_PSC_F(1,:)+dt*ROn_TD_PSC_F_k1;
  ROn_TD_PSC_P(2,:) = ROn_TD_PSC_P(1,:)+dt*ROn_TD_PSC_P_k1;
  ROn_TD_PSC_q(2,:) = ROn_TD_PSC_q(1,:)+dt*ROn_TD_PSC_q_k1;
  ROff_TD_PSC_s(2,:) = ROff_TD_PSC_s(1,:)+dt*ROff_TD_PSC_s_k1;
  ROff_TD_PSC_x(2,:) = ROff_TD_PSC_x(1,:)+dt*ROff_TD_PSC_x_k1;
  ROff_TD_PSC_F(2,:) = ROff_TD_PSC_F(1,:)+dt*ROff_TD_PSC_F_k1;
  ROff_TD_PSC_P(2,:) = ROff_TD_PSC_P(1,:)+dt*ROff_TD_PSC_P_k1;
  ROff_TD_PSC_q(2,:) = ROff_TD_PSC_q(1,:)+dt*ROff_TD_PSC_q_k1;
  X_TD_PSC_s(2,:) = X_TD_PSC_s(1,:)+dt*X_TD_PSC_s_k1;
  X_TD_PSC_x(2,:) = X_TD_PSC_x(1,:)+dt*X_TD_PSC_x_k1;
  X_TD_PSC_F(2,:) = X_TD_PSC_F(1,:)+dt*X_TD_PSC_F_k1;
  X_TD_PSC_P(2,:) = X_TD_PSC_P(1,:)+dt*X_TD_PSC_P_k1;
  X_TD_PSC_q(2,:) = X_TD_PSC_q(1,:)+dt*X_TD_PSC_q_k1;
  C_ROn_PSC3_s(2,:) = C_ROn_PSC3_s(1,:)+dt*C_ROn_PSC3_s_k1;
  C_ROn_PSC3_x(2,:) = C_ROn_PSC3_x(1,:)+dt*C_ROn_PSC3_x_k1;
  C_ROn_PSC3_F(2,:) = C_ROn_PSC3_F(1,:)+dt*C_ROn_PSC3_F_k1;
  C_ROn_PSC3_P(2,:) = C_ROn_PSC3_P(1,:)+dt*C_ROn_PSC3_P_k1;
  C_ROn_PSC3_q(2,:) = C_ROn_PSC3_q(1,:)+dt*C_ROn_PSC3_q_k1;

  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=any(C_V(2,:)>=C_V_thresh&C_V(1,:)<C_V_thresh);
  conditional_indx=(C_V(2,:)>=C_V_thresh&C_V(1,:)<C_V_thresh);
  if conditional_test, C_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); C_tspike(C_buffer_index(i),i)=t; C_buffer_index(i)=mod(-1+(C_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(X_V(2,:)>=X_V_thresh&X_V(1,:)<X_V_thresh);
  conditional_indx=(X_V(2,:)>=X_V_thresh&X_V(1,:)<X_V_thresh);
  if conditional_test, X_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); X_tspike(X_buffer_index(i),i)=t; X_buffer_index(i)=mod(-1+(X_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(TD_V(2,:)>=TD_V_thresh&TD_V(1,:)<TD_V_thresh);
  conditional_indx=(TD_V(2,:)>=TD_V_thresh&TD_V(1,:)<TD_V_thresh);
  if conditional_test, TD_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); TD_tspike(TD_buffer_index(i),i)=t; TD_buffer_index(i)=mod(-1+(TD_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(SOnOff_V(2,:)>=SOnOff_V_thresh&SOnOff_V(1,:)<SOnOff_V_thresh);
  conditional_indx=(SOnOff_V(2,:)>=SOnOff_V_thresh&SOnOff_V(1,:)<SOnOff_V_thresh);
  if conditional_test, SOnOff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); SOnOff_tspike(SOnOff_buffer_index(i),i)=t; SOnOff_buffer_index(i)=mod(-1+(SOnOff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(ROff_V(2,:)>=ROff_V_thresh&ROff_V(1,:)<ROff_V_thresh);
  conditional_indx=(ROff_V(2,:)>=ROff_V_thresh&ROff_V(1,:)<ROff_V_thresh);
  if conditional_test, ROff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); ROff_tspike(ROff_buffer_index(i),i)=t; ROff_buffer_index(i)=mod(-1+(ROff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(ROn_V(2,:)>=ROn_V_thresh&ROn_V(1,:)<ROn_V_thresh);
  conditional_indx=(ROn_V(2,:)>=ROn_V_thresh&ROn_V(1,:)<ROn_V_thresh);
  if conditional_test, ROn_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); ROn_tspike(ROn_buffer_index(i),i)=t; ROn_buffer_index(i)=mod(-1+(ROn_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(Off_V(2,:)>=Off_V_thresh&Off_V(1,:)<Off_V_thresh);
  conditional_indx=(Off_V(2,:)>=Off_V_thresh&Off_V(1,:)<Off_V_thresh);
  if conditional_test, Off_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); Off_tspike(Off_buffer_index(i),i)=t; Off_buffer_index(i)=mod(-1+(Off_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(On_V(2,:)>=On_V_thresh&On_V(1,:)<On_V_thresh);
  conditional_indx=(On_V(2,:)>=On_V_thresh&On_V(1,:)<On_V_thresh);
  if conditional_test, On_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); On_tspike(On_buffer_index(i),i)=t; On_buffer_index(i)=mod(-1+(On_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(On_V(2,:) > On_V_thresh);
  conditional_indx=(On_V(2,:) > On_V_thresh);
  if conditional_test, On_V(2,conditional_indx) = On_V_reset; On_g_ad(2,conditional_indx) = On_g_ad(2,conditional_indx) + On_g_inc; end
  conditional_test=any(any(t<=On_tspike+On_t_ref,1));
  conditional_indx=(any(t<=On_tspike+On_t_ref,1));
  if conditional_test, On_V(2,conditional_indx) = On_V_reset; end
  conditional_test=any(Off_V(2,:) > Off_V_thresh);
  conditional_indx=(Off_V(2,:) > Off_V_thresh);
  if conditional_test, Off_V(2,conditional_indx) = Off_V_reset; Off_g_ad(2,conditional_indx) = Off_g_ad(2,conditional_indx) + Off_g_inc; end
  conditional_test=any(any(t<=Off_tspike+Off_t_ref,1));
  conditional_indx=(any(t<=Off_tspike+Off_t_ref,1));
  if conditional_test, Off_V(2,conditional_indx) = Off_V_reset; end
  conditional_test=any(ROn_V(2,:) > ROn_V_thresh);
  conditional_indx=(ROn_V(2,:) > ROn_V_thresh);
  if conditional_test, ROn_V(2,conditional_indx) = ROn_V_reset; ROn_g_ad(2,conditional_indx) = ROn_g_ad(2,conditional_indx) + ROn_g_inc; end
  conditional_test=any(any(t<=ROn_tspike+ROn_t_ref,1));
  conditional_indx=(any(t<=ROn_tspike+ROn_t_ref,1));
  if conditional_test, ROn_V(2,conditional_indx) = ROn_V_reset; end
  conditional_test=any(ROff_V(2,:) > ROff_V_thresh);
  conditional_indx=(ROff_V(2,:) > ROff_V_thresh);
  if conditional_test, ROff_V(2,conditional_indx) = ROff_V_reset; ROff_g_ad(2,conditional_indx) = ROff_g_ad(2,conditional_indx) + ROff_g_inc; end
  conditional_test=any(any(t<=ROff_tspike+ROff_t_ref,1));
  conditional_indx=(any(t<=ROff_tspike+ROff_t_ref,1));
  if conditional_test, ROff_V(2,conditional_indx) = ROff_V_reset; end
  conditional_test=any(SOnOff_V(2,:) > SOnOff_V_thresh);
  conditional_indx=(SOnOff_V(2,:) > SOnOff_V_thresh);
  if conditional_test, SOnOff_V(2,conditional_indx) = SOnOff_V_reset; SOnOff_g_ad(2,conditional_indx) = SOnOff_g_ad(2,conditional_indx) + SOnOff_g_inc; end
  conditional_test=any(any(t<=SOnOff_tspike+SOnOff_t_ref,1));
  conditional_indx=(any(t<=SOnOff_tspike+SOnOff_t_ref,1));
  if conditional_test, SOnOff_V(2,conditional_indx) = SOnOff_V_reset; end
  conditional_test=any(TD_V(2,:) > TD_V_thresh);
  conditional_indx=(TD_V(2,:) > TD_V_thresh);
  if conditional_test, TD_V(2,conditional_indx) = TD_V_reset; TD_g_ad(2,conditional_indx) = TD_g_ad(2,conditional_indx) + TD_g_inc; end
  conditional_test=any(any(t<=TD_tspike+TD_t_ref,1));
  conditional_indx=(any(t<=TD_tspike+TD_t_ref,1));
  if conditional_test, TD_V(2,conditional_indx) = TD_V_reset; end
  conditional_test=any(X_V(2,:) > X_V_thresh);
  conditional_indx=(X_V(2,:) > X_V_thresh);
  if conditional_test, X_V(2,conditional_indx) = X_V_reset; X_g_ad(2,conditional_indx) = X_g_ad(2,conditional_indx) + X_g_inc; end
  conditional_test=any(any(t<=X_tspike+X_t_ref,1));
  conditional_indx=(any(t<=X_tspike+X_t_ref,1));
  if conditional_test, X_V(2,conditional_indx) = X_V_reset; end
  conditional_test=any(C_V(2,:) > C_V_thresh);
  conditional_indx=(C_V(2,:) > C_V_thresh);
  if conditional_test, C_V(2,conditional_indx) = C_V_reset; C_g_ad(2,conditional_indx) = C_g_ad(2,conditional_indx) + C_g_inc; end
  conditional_test=any(any(t<=C_tspike+C_t_ref,1));
  conditional_indx=(any(t<=C_tspike+C_t_ref,1));
  if conditional_test, C_V(2,conditional_indx) = C_V_reset; end
  conditional_test=any(any(t == On_tspike+ROn_On_PSC3_delay,1));
  conditional_indx=(any(t == On_tspike+ROn_On_PSC3_delay,1));
  if conditional_test, ROn_On_PSC3_x(2,conditional_indx) = ROn_On_PSC3_x(2,conditional_indx) + ROn_On_PSC3_q(2,conditional_indx);ROn_On_PSC3_q(2,conditional_indx) = ROn_On_PSC3_F(2,conditional_indx).*ROn_On_PSC3_P(2,conditional_indx);ROn_On_PSC3_F(2,conditional_indx) = ROn_On_PSC3_F(2,conditional_indx) + ROn_On_PSC3_fF*(ROn_On_PSC3_maxF-ROn_On_PSC3_F(2,conditional_indx)); ROn_On_PSC3_P(2,conditional_indx) = ROn_On_PSC3_P(2,conditional_indx)*(1 - ROn_On_PSC3_fP); end
  conditional_test=any(any(t == On_tspike+SOnOff_On_PSC_delay,1));
  conditional_indx=(any(t == On_tspike+SOnOff_On_PSC_delay,1));
  if conditional_test, SOnOff_On_PSC_x(2,conditional_indx) = SOnOff_On_PSC_x(2,conditional_indx) + SOnOff_On_PSC_q(2,conditional_indx);SOnOff_On_PSC_q(2,conditional_indx) = SOnOff_On_PSC_F(2,conditional_indx).*SOnOff_On_PSC_P(2,conditional_indx);SOnOff_On_PSC_F(2,conditional_indx) = SOnOff_On_PSC_F(2,conditional_indx) + SOnOff_On_PSC_fF*(SOnOff_On_PSC_maxF-SOnOff_On_PSC_F(2,conditional_indx)); SOnOff_On_PSC_P(2,conditional_indx) = SOnOff_On_PSC_P(2,conditional_indx)*(1 - SOnOff_On_PSC_fP); end
  conditional_test=any(any(t == SOnOff_tspike+ROn_SOnOff_PSC3_delay,1));
  conditional_indx=(any(t == SOnOff_tspike+ROn_SOnOff_PSC3_delay,1));
  if conditional_test, ROn_SOnOff_PSC3_x(2,conditional_indx) = ROn_SOnOff_PSC3_x(2,conditional_indx) + ROn_SOnOff_PSC3_q(2,conditional_indx);ROn_SOnOff_PSC3_q(2,conditional_indx) = ROn_SOnOff_PSC3_F(2,conditional_indx).*ROn_SOnOff_PSC3_P(2,conditional_indx);ROn_SOnOff_PSC3_F(2,conditional_indx) = ROn_SOnOff_PSC3_F(2,conditional_indx) + ROn_SOnOff_PSC3_fF*(ROn_SOnOff_PSC3_maxF-ROn_SOnOff_PSC3_F(2,conditional_indx)); ROn_SOnOff_PSC3_P(2,conditional_indx) = ROn_SOnOff_PSC3_P(2,conditional_indx)*(1 - ROn_SOnOff_PSC3_fP); end
  conditional_test=any(any(t == SOnOff_tspike+ROff_SOnOff_PSC_delay,1));
  conditional_indx=(any(t == SOnOff_tspike+ROff_SOnOff_PSC_delay,1));
  if conditional_test, ROff_SOnOff_PSC_x(2,conditional_indx) = ROff_SOnOff_PSC_x(2,conditional_indx) + ROff_SOnOff_PSC_q(2,conditional_indx);ROff_SOnOff_PSC_q(2,conditional_indx) = ROff_SOnOff_PSC_F(2,conditional_indx).*ROff_SOnOff_PSC_P(2,conditional_indx);ROff_SOnOff_PSC_F(2,conditional_indx) = ROff_SOnOff_PSC_F(2,conditional_indx) + ROff_SOnOff_PSC_fF*(ROff_SOnOff_PSC_maxF-ROff_SOnOff_PSC_F(2,conditional_indx)); ROff_SOnOff_PSC_P(2,conditional_indx) = ROff_SOnOff_PSC_P(2,conditional_indx)*(1 - ROff_SOnOff_PSC_fP); end
  conditional_test=any(any(t == Off_tspike+ROff_Off_PSC_delay,1));
  conditional_indx=(any(t == Off_tspike+ROff_Off_PSC_delay,1));
  if conditional_test, ROff_Off_PSC_x(2,conditional_indx) = ROff_Off_PSC_x(2,conditional_indx) + ROff_Off_PSC_q(2,conditional_indx);ROff_Off_PSC_q(2,conditional_indx) = ROff_Off_PSC_F(2,conditional_indx).*ROff_Off_PSC_P(2,conditional_indx);ROff_Off_PSC_F(2,conditional_indx) = ROff_Off_PSC_F(2,conditional_indx) + ROff_Off_PSC_fF*(ROff_Off_PSC_maxF-ROff_Off_PSC_F(2,conditional_indx)); ROff_Off_PSC_P(2,conditional_indx) = ROff_Off_PSC_P(2,conditional_indx)*(1 - ROff_Off_PSC_fP); end
  conditional_test=any(any(t == Off_tspike+SOnOff_Off_PSC_delay,1));
  conditional_indx=(any(t == Off_tspike+SOnOff_Off_PSC_delay,1));
  if conditional_test, SOnOff_Off_PSC_x(2,conditional_indx) = SOnOff_Off_PSC_x(2,conditional_indx) + SOnOff_Off_PSC_q(2,conditional_indx);SOnOff_Off_PSC_q(2,conditional_indx) = SOnOff_Off_PSC_F(2,conditional_indx).*SOnOff_Off_PSC_P(2,conditional_indx);SOnOff_Off_PSC_F(2,conditional_indx) = SOnOff_Off_PSC_F(2,conditional_indx) + SOnOff_Off_PSC_fF*(SOnOff_Off_PSC_maxF-SOnOff_Off_PSC_F(2,conditional_indx)); SOnOff_Off_PSC_P(2,conditional_indx) = SOnOff_Off_PSC_P(2,conditional_indx)*(1 - SOnOff_Off_PSC_fP); end
  conditional_test=any(any(t == ROn_tspike+X_ROn_PSC_delay,1));
  conditional_indx=(any(t == ROn_tspike+X_ROn_PSC_delay,1));
  if conditional_test, X_ROn_PSC_x(2,conditional_indx) = X_ROn_PSC_x(2,conditional_indx) + X_ROn_PSC_q(2,conditional_indx);X_ROn_PSC_q(2,conditional_indx) = X_ROn_PSC_F(2,conditional_indx).*X_ROn_PSC_P(2,conditional_indx);X_ROn_PSC_F(2,conditional_indx) = X_ROn_PSC_F(2,conditional_indx) + X_ROn_PSC_fF*(X_ROn_PSC_maxF-X_ROn_PSC_F(2,conditional_indx)); X_ROn_PSC_P(2,conditional_indx) = X_ROn_PSC_P(2,conditional_indx)*(1 - X_ROn_PSC_fP); end
  conditional_test=any(any(t == X_tspike+ROn_X_PSC3_delay,1));
  conditional_indx=(any(t == X_tspike+ROn_X_PSC3_delay,1));
  if conditional_test, ROn_X_PSC3_x(2,conditional_indx) = ROn_X_PSC3_x(2,conditional_indx) + ROn_X_PSC3_q(2,conditional_indx);ROn_X_PSC3_q(2,conditional_indx) = ROn_X_PSC3_F(2,conditional_indx).*ROn_X_PSC3_P(2,conditional_indx);ROn_X_PSC3_F(2,conditional_indx) = ROn_X_PSC3_F(2,conditional_indx) + ROn_X_PSC3_fF*(ROn_X_PSC3_maxF-ROn_X_PSC3_F(2,conditional_indx)); ROn_X_PSC3_P(2,conditional_indx) = ROn_X_PSC3_P(2,conditional_indx)*(1 - ROn_X_PSC3_fP); end
  conditional_test=any(any(t == TD_tspike+ROn_TD_PSC_delay,1));
  conditional_indx=(any(t == TD_tspike+ROn_TD_PSC_delay,1));
  if conditional_test, ROn_TD_PSC_x(2,conditional_indx) = ROn_TD_PSC_x(2,conditional_indx) + ROn_TD_PSC_q(2,conditional_indx);ROn_TD_PSC_q(2,conditional_indx) = ROn_TD_PSC_F(2,conditional_indx).*ROn_TD_PSC_P(2,conditional_indx);ROn_TD_PSC_F(2,conditional_indx) = ROn_TD_PSC_F(2,conditional_indx) + ROn_TD_PSC_fF*(ROn_TD_PSC_maxF-ROn_TD_PSC_F(2,conditional_indx)); ROn_TD_PSC_P(2,conditional_indx) = ROn_TD_PSC_P(2,conditional_indx)*(1 - ROn_TD_PSC_fP); end
  conditional_test=any(any(t == TD_tspike+ROff_TD_PSC_delay,1));
  conditional_indx=(any(t == TD_tspike+ROff_TD_PSC_delay,1));
  if conditional_test, ROff_TD_PSC_x(2,conditional_indx) = ROff_TD_PSC_x(2,conditional_indx) + ROff_TD_PSC_q(2,conditional_indx);ROff_TD_PSC_q(2,conditional_indx) = ROff_TD_PSC_F(2,conditional_indx).*ROff_TD_PSC_P(2,conditional_indx);ROff_TD_PSC_F(2,conditional_indx) = ROff_TD_PSC_F(2,conditional_indx) + ROff_TD_PSC_fF*(ROff_TD_PSC_maxF-ROff_TD_PSC_F(2,conditional_indx)); ROff_TD_PSC_P(2,conditional_indx) = ROff_TD_PSC_P(2,conditional_indx)*(1 - ROff_TD_PSC_fP); end
  conditional_test=any(any(t == TD_tspike+X_TD_PSC_delay,1));
  conditional_indx=(any(t == TD_tspike+X_TD_PSC_delay,1));
  if conditional_test, X_TD_PSC_x(2,conditional_indx) = X_TD_PSC_x(2,conditional_indx) + X_TD_PSC_q(2,conditional_indx);X_TD_PSC_q(2,conditional_indx) = X_TD_PSC_F(2,conditional_indx).*X_TD_PSC_P(2,conditional_indx);X_TD_PSC_F(2,conditional_indx) = X_TD_PSC_F(2,conditional_indx) + X_TD_PSC_fF*(X_TD_PSC_maxF-X_TD_PSC_F(2,conditional_indx)); X_TD_PSC_P(2,conditional_indx) = X_TD_PSC_P(2,conditional_indx)*(1 - X_TD_PSC_fP); end
  conditional_test=any(any(t == ROn_tspike+C_ROn_PSC3_delay,1));
  conditional_indx=(any(t == ROn_tspike+C_ROn_PSC3_delay,1));
  if conditional_test, C_ROn_PSC3_x(2,conditional_indx) = C_ROn_PSC3_x(2,conditional_indx) + C_ROn_PSC3_q(2,conditional_indx);C_ROn_PSC3_q(2,conditional_indx) = C_ROn_PSC3_F(2,conditional_indx).*C_ROn_PSC3_P(2,conditional_indx);C_ROn_PSC3_F(2,conditional_indx) = C_ROn_PSC3_F(2,conditional_indx) + C_ROn_PSC3_fF*(C_ROn_PSC3_maxF-C_ROn_PSC3_F(2,conditional_indx)); C_ROn_PSC3_P(2,conditional_indx) = C_ROn_PSC3_P(2,conditional_indx)*(1 - C_ROn_PSC3_fP); end

  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  On_On_IC_iIC(n,:)=On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(2,:)-On_On_IC_E_exc);
  Off_Off_IC_iIC(n,:)=Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(2,:)-Off_Off_IC_E_exc);
  ROn_On_PSC3_syn(n,:)=((ROn_On_PSC3_s(2,:)*(ROn_On_PSC3_netcon.*ROn_On_PSC3_gSYN)).*(ROn_V(2,:)-ROn_On_PSC3_ESYN));
  SOnOff_On_PSC_syn(n,:)=SOnOff_On_PSC_gSYN.*(SOnOff_On_PSC_s(2,:)*SOnOff_On_PSC_netcon).*(SOnOff_V(2,:)-SOnOff_On_PSC_ESYN);
  ROn_SOnOff_PSC3_syn(n,:)=((ROn_SOnOff_PSC3_s(2,:)*(ROn_SOnOff_PSC3_netcon.*ROn_SOnOff_PSC3_gSYN)).*(ROn_V(2,:)-ROn_SOnOff_PSC3_ESYN));
  ROff_SOnOff_PSC_syn(n,:)=ROff_SOnOff_PSC_gSYN.*(ROff_SOnOff_PSC_s(2,:)*ROff_SOnOff_PSC_netcon).*(ROff_V(2,:)-ROff_SOnOff_PSC_ESYN);
  ROff_Off_PSC_syn(n,:)=ROff_Off_PSC_gSYN.*(ROff_Off_PSC_s(2,:)*ROff_Off_PSC_netcon).*(ROff_V(2,:)-ROff_Off_PSC_ESYN);
  SOnOff_Off_PSC_syn(n,:)=SOnOff_Off_PSC_gSYN.*(SOnOff_Off_PSC_s(2,:)*SOnOff_Off_PSC_netcon).*(SOnOff_V(2,:)-SOnOff_Off_PSC_ESYN);
  X_ROn_PSC_syn(n,:)=X_ROn_PSC_gSYN.*(X_ROn_PSC_s(2,:)*X_ROn_PSC_netcon).*(X_V(2,:)-X_ROn_PSC_ESYN);
  ROn_X_PSC3_syn(n,:)=((ROn_X_PSC3_s(2,:)*(ROn_X_PSC3_netcon.*ROn_X_PSC3_gSYN)).*(ROn_V(2,:)-ROn_X_PSC3_ESYN));
  ROn_TD_PSC_syn(n,:)=ROn_TD_PSC_gSYN.*(ROn_TD_PSC_s(2,:)*ROn_TD_PSC_netcon).*(ROn_V(2,:)-ROn_TD_PSC_ESYN);
  ROff_TD_PSC_syn(n,:)=ROff_TD_PSC_gSYN.*(ROff_TD_PSC_s(2,:)*ROff_TD_PSC_netcon).*(ROff_V(2,:)-ROff_TD_PSC_ESYN);
  X_TD_PSC_syn(n,:)=X_TD_PSC_gSYN.*(X_TD_PSC_s(2,:)*X_TD_PSC_netcon).*(X_V(2,:)-X_TD_PSC_ESYN);
  C_ROn_PSC3_syn(n,:)=((C_ROn_PSC3_s(2,:)*(C_ROn_PSC3_netcon.*C_ROn_PSC3_gSYN)).*(C_V(2,:)-C_ROn_PSC3_ESYN));

  % ------------------------------------------------------------
  % Replace n=1 for state variables IB:
  % ------------------------------------------------------------
On_V(1,:)=On_V(2,:);
On_g_ad(1,:)=On_g_ad(2,:);
Off_V(1,:)=Off_V(2,:);
Off_g_ad(1,:)=Off_g_ad(2,:);
ROn_V(1,:)=ROn_V(2,:);
ROn_g_ad(1,:)=ROn_g_ad(2,:);
ROff_V(1,:)=ROff_V(2,:);
ROff_g_ad(1,:)=ROff_g_ad(2,:);
SOnOff_V(1,:)=SOnOff_V(2,:);
SOnOff_g_ad(1,:)=SOnOff_g_ad(2,:);
TD_V(1,:)=TD_V(2,:);
TD_g_ad(1,:)=TD_g_ad(2,:);
X_V(1,:)=X_V(2,:);
X_g_ad(1,:)=X_g_ad(2,:);
C_V(1,:)=C_V(2,:);
C_g_ad(1,:)=C_g_ad(2,:);
ROn_On_PSC3_s(1,:)=ROn_On_PSC3_s(2,:);
ROn_On_PSC3_x(1,:)=ROn_On_PSC3_x(2,:);
ROn_On_PSC3_F(1,:)=ROn_On_PSC3_F(2,:);
ROn_On_PSC3_P(1,:)=ROn_On_PSC3_P(2,:);
ROn_On_PSC3_q(1,:)=ROn_On_PSC3_q(2,:);
SOnOff_On_PSC_s(1,:)=SOnOff_On_PSC_s(2,:);
SOnOff_On_PSC_x(1,:)=SOnOff_On_PSC_x(2,:);
SOnOff_On_PSC_F(1,:)=SOnOff_On_PSC_F(2,:);
SOnOff_On_PSC_P(1,:)=SOnOff_On_PSC_P(2,:);
SOnOff_On_PSC_q(1,:)=SOnOff_On_PSC_q(2,:);
ROn_SOnOff_PSC3_s(1,:)=ROn_SOnOff_PSC3_s(2,:);
ROn_SOnOff_PSC3_x(1,:)=ROn_SOnOff_PSC3_x(2,:);
ROn_SOnOff_PSC3_F(1,:)=ROn_SOnOff_PSC3_F(2,:);
ROn_SOnOff_PSC3_P(1,:)=ROn_SOnOff_PSC3_P(2,:);
ROn_SOnOff_PSC3_q(1,:)=ROn_SOnOff_PSC3_q(2,:);
ROff_SOnOff_PSC_s(1,:)=ROff_SOnOff_PSC_s(2,:);
ROff_SOnOff_PSC_x(1,:)=ROff_SOnOff_PSC_x(2,:);
ROff_SOnOff_PSC_F(1,:)=ROff_SOnOff_PSC_F(2,:);
ROff_SOnOff_PSC_P(1,:)=ROff_SOnOff_PSC_P(2,:);
ROff_SOnOff_PSC_q(1,:)=ROff_SOnOff_PSC_q(2,:);
ROff_Off_PSC_s(1,:)=ROff_Off_PSC_s(2,:);
ROff_Off_PSC_x(1,:)=ROff_Off_PSC_x(2,:);
ROff_Off_PSC_F(1,:)=ROff_Off_PSC_F(2,:);
ROff_Off_PSC_P(1,:)=ROff_Off_PSC_P(2,:);
ROff_Off_PSC_q(1,:)=ROff_Off_PSC_q(2,:);
SOnOff_Off_PSC_s(1,:)=SOnOff_Off_PSC_s(2,:);
SOnOff_Off_PSC_x(1,:)=SOnOff_Off_PSC_x(2,:);
SOnOff_Off_PSC_F(1,:)=SOnOff_Off_PSC_F(2,:);
SOnOff_Off_PSC_P(1,:)=SOnOff_Off_PSC_P(2,:);
SOnOff_Off_PSC_q(1,:)=SOnOff_Off_PSC_q(2,:);
ROn_ROn_iNoise_V3_sn(1,:)=ROn_ROn_iNoise_V3_sn(2,:);
ROn_ROn_iNoise_V3_xn(1,:)=ROn_ROn_iNoise_V3_xn(2,:);
X_ROn_PSC_s(1,:)=X_ROn_PSC_s(2,:);
X_ROn_PSC_x(1,:)=X_ROn_PSC_x(2,:);
X_ROn_PSC_F(1,:)=X_ROn_PSC_F(2,:);
X_ROn_PSC_P(1,:)=X_ROn_PSC_P(2,:);
X_ROn_PSC_q(1,:)=X_ROn_PSC_q(2,:);
ROn_X_PSC3_s(1,:)=ROn_X_PSC3_s(2,:);
ROn_X_PSC3_x(1,:)=ROn_X_PSC3_x(2,:);
ROn_X_PSC3_F(1,:)=ROn_X_PSC3_F(2,:);
ROn_X_PSC3_P(1,:)=ROn_X_PSC3_P(2,:);
ROn_X_PSC3_q(1,:)=ROn_X_PSC3_q(2,:);
ROn_TD_PSC_s(1,:)=ROn_TD_PSC_s(2,:);
ROn_TD_PSC_x(1,:)=ROn_TD_PSC_x(2,:);
ROn_TD_PSC_F(1,:)=ROn_TD_PSC_F(2,:);
ROn_TD_PSC_P(1,:)=ROn_TD_PSC_P(2,:);
ROn_TD_PSC_q(1,:)=ROn_TD_PSC_q(2,:);
ROff_TD_PSC_s(1,:)=ROff_TD_PSC_s(2,:);
ROff_TD_PSC_x(1,:)=ROff_TD_PSC_x(2,:);
ROff_TD_PSC_F(1,:)=ROff_TD_PSC_F(2,:);
ROff_TD_PSC_P(1,:)=ROff_TD_PSC_P(2,:);
ROff_TD_PSC_q(1,:)=ROff_TD_PSC_q(2,:);
X_TD_PSC_s(1,:)=X_TD_PSC_s(2,:);
X_TD_PSC_x(1,:)=X_TD_PSC_x(2,:);
X_TD_PSC_F(1,:)=X_TD_PSC_F(2,:);
X_TD_PSC_P(1,:)=X_TD_PSC_P(2,:);
X_TD_PSC_q(1,:)=X_TD_PSC_q(2,:);
C_ROn_PSC3_s(1,:)=C_ROn_PSC3_s(2,:);
C_ROn_PSC3_x(1,:)=C_ROn_PSC3_x(2,:);
C_ROn_PSC3_F(1,:)=C_ROn_PSC3_F(2,:);
C_ROn_PSC3_P(1,:)=C_ROn_PSC3_P(2,:);
C_ROn_PSC3_q(1,:)=C_ROn_PSC3_q(2,:);
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end
