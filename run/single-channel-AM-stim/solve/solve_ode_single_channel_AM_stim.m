function [T,On_V,On_g_ad,Off_V,Off_g_ad,R1On_V,R1On_g_ad,R1Off_V,R1Off_g_ad,S1OnOff_V,S1OnOff_g_ad,R2On_V,R2On_g_ad,R2Off_V,R2Off_g_ad,S2OnOff_V,S2OnOff_g_ad,R1On_On_PSC_s,R1On_On_PSC_x,R1On_On_PSC_F,R1On_On_PSC_P,R1On_On_PSC_q,S1OnOff_On_PSC_s,S1OnOff_On_PSC_x,S1OnOff_On_PSC_F,S1OnOff_On_PSC_P,S1OnOff_On_PSC_q,R1On_S1OnOff_PSC_s,R1On_S1OnOff_PSC_x,R1On_S1OnOff_PSC_F,R1On_S1OnOff_PSC_P,R1On_S1OnOff_PSC_q,R1Off_S1OnOff_PSC_s,R1Off_S1OnOff_PSC_x,R1Off_S1OnOff_PSC_F,R1Off_S1OnOff_PSC_P,R1Off_S1OnOff_PSC_q,R1Off_Off_PSC_s,R1Off_Off_PSC_x,R1Off_Off_PSC_F,R1Off_Off_PSC_P,R1Off_Off_PSC_q,S1OnOff_Off_PSC_s,S1OnOff_Off_PSC_x,S1OnOff_Off_PSC_F,S1OnOff_Off_PSC_P,S1OnOff_Off_PSC_q,R2On_R1On_PSC_s,R2On_R1On_PSC_x,R2On_R1On_PSC_F,R2On_R1On_PSC_P,R2On_R1On_PSC_q,S2OnOff_R1On_PSC_s,S2OnOff_R1On_PSC_x,S2OnOff_R1On_PSC_F,S2OnOff_R1On_PSC_P,S2OnOff_R1On_PSC_q,R2On_S2OnOff_PSC_s,R2On_S2OnOff_PSC_x,R2On_S2OnOff_PSC_F,R2On_S2OnOff_PSC_P,R2On_S2OnOff_PSC_q,R2Off_S2OnOff_PSC_s,R2Off_S2OnOff_PSC_x,R2Off_S2OnOff_PSC_F,R2Off_S2OnOff_PSC_P,R2Off_S2OnOff_PSC_q,R2Off_R1Off_PSC_s,R2Off_R1Off_PSC_x,R2Off_R1Off_PSC_F,R2Off_R1Off_PSC_P,R2Off_R1Off_PSC_q,S2OnOff_R1Off_PSC_s,S2OnOff_R1Off_PSC_x,S2OnOff_R1Off_PSC_F,S2OnOff_R1Off_PSC_P,S2OnOff_R1Off_PSC_q,R2On_R2On_iNoise_V3_sn,R2On_R2On_iNoise_V3_xn,On_V_spikes,Off_V_spikes,R1On_V_spikes,R1Off_V_spikes,S1OnOff_V_spikes,R2On_V_spikes,R2Off_V_spikes,S2OnOff_V_spikes,On_On_IC_iIC,Off_Off_IC_iIC,R1On_On_PSC_syn,S1OnOff_On_PSC_syn,R1On_S1OnOff_PSC_syn,R1Off_S1OnOff_PSC_syn,R1Off_Off_PSC_syn,S1OnOff_Off_PSC_syn,R2On_R1On_PSC_syn,S2OnOff_R1On_PSC_syn,R2On_S2OnOff_PSC_syn,R2Off_S2OnOff_PSC_syn,R2Off_R1Off_PSC_syn,S2OnOff_R1Off_PSC_syn,On_R,On_tau,On_Imask,Off_R,Off_tau,Off_Imask,R1On_R,R1On_tau,R1On_Imask,R1Off_R,R1Off_tau,R1Off_Imask,S1OnOff_R,S1OnOff_tau,S1OnOff_Imask,R2On_R,R2On_tau,R2On_Imask,R2Off_R,R2Off_tau,R2Off_Imask,S2OnOff_R,S2OnOff_tau,S2OnOff_Imask,On_On_IC_netcon,On_On_IC_input,Off_Off_IC_netcon,Off_Off_IC_input,R1On_On_PSC_netcon,R1On_On_PSC_scale,S1OnOff_On_PSC_netcon,S1OnOff_On_PSC_scale,R1On_S1OnOff_PSC_netcon,R1On_S1OnOff_PSC_scale,R1Off_S1OnOff_PSC_netcon,R1Off_S1OnOff_PSC_scale,R1Off_Off_PSC_netcon,R1Off_Off_PSC_scale,S1OnOff_Off_PSC_netcon,S1OnOff_Off_PSC_scale,R2On_R1On_PSC_netcon,R2On_R1On_PSC_scale,S2OnOff_R1On_PSC_netcon,S2OnOff_R1On_PSC_scale,R2On_S2OnOff_PSC_netcon,R2On_S2OnOff_PSC_scale,R2Off_S2OnOff_PSC_netcon,R2Off_S2OnOff_PSC_scale,R2Off_R1Off_PSC_netcon,R2Off_R1Off_PSC_scale,S2OnOff_R1Off_PSC_netcon,S2OnOff_R1Off_PSC_scale,R2On_R2On_iNoise_V3_netcon,R2On_R2On_iNoise_V3_token,R2On_R2On_iNoise_V3_scale]=solve_ode

% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);

% the 'shuffle' random seed has been set in advance.

% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
On_R = 1/p.On_g_L;
On_tau = p.On_C*On_R;
On_Imask =  ones(1,p.On_Npop);
Off_R = 1/p.Off_g_L;
Off_tau = p.Off_C*Off_R;
Off_Imask =  ones(1,p.Off_Npop);
R1On_R = 1/p.R1On_g_L;
R1On_tau = p.R1On_C*R1On_R;
R1On_Imask =  ones(1,p.R1On_Npop);
R1Off_R = 1/p.R1Off_g_L;
R1Off_tau = p.R1Off_C*R1Off_R;
R1Off_Imask =  ones(1,p.R1Off_Npop);
S1OnOff_R = 1/p.S1OnOff_g_L;
S1OnOff_tau = p.S1OnOff_C*S1OnOff_R;
S1OnOff_Imask =  ones(1,p.S1OnOff_Npop);
R2On_R = 1/p.R2On_g_L;
R2On_tau = p.R2On_C*R2On_R;
R2On_Imask =  ones(1,p.R2On_Npop);
R2Off_R = 1/p.R2Off_g_L;
R2Off_tau = p.R2Off_C*R2Off_R;
R2Off_Imask =  ones(1,p.R2Off_Npop);
S2OnOff_R = 1/p.S2OnOff_g_L;
S2OnOff_tau = p.S2OnOff_C*S2OnOff_R;
S2OnOff_Imask =  ones(1,p.S2OnOff_Npop);
On_On_IC_netcon = +1.000000000000000e+00;
On_On_IC_input =  genPoissonInputs(p.On_On_IC_trial,p.On_On_IC_locNum,p.On_On_IC_label,p.On_On_IC_t_ref,p.On_On_IC_t_ref_rel,p.On_On_IC_rec);
Off_Off_IC_netcon = +1.000000000000000e+00;
Off_Off_IC_input =  genPoissonInputs(p.Off_Off_IC_trial,p.Off_Off_IC_locNum,p.Off_Off_IC_label,p.Off_Off_IC_t_ref,p.Off_Off_IC_t_ref_rel,p.Off_Off_IC_rec);
R1On_On_PSC_netcon = eye(p.On_Npop,p.R1On_Npop);
R1On_On_PSC_scale = (p.R1On_On_PSC_tauD/p.R1On_On_PSC_tauR)^(p.R1On_On_PSC_tauR/(p.R1On_On_PSC_tauD-p.R1On_On_PSC_tauR));
S1OnOff_On_PSC_netcon = eye(p.On_Npop,p.S1OnOff_Npop);
S1OnOff_On_PSC_scale = (p.S1OnOff_On_PSC_tauD/p.S1OnOff_On_PSC_tauR)^(p.S1OnOff_On_PSC_tauR/(p.S1OnOff_On_PSC_tauD-p.S1OnOff_On_PSC_tauR));
R1On_S1OnOff_PSC_netcon = eye(p.S1OnOff_Npop,p.R1On_Npop);
R1On_S1OnOff_PSC_scale = (p.R1On_S1OnOff_PSC_tauD/p.R1On_S1OnOff_PSC_tauR)^(p.R1On_S1OnOff_PSC_tauR/(p.R1On_S1OnOff_PSC_tauD-p.R1On_S1OnOff_PSC_tauR));
R1Off_S1OnOff_PSC_netcon = eye(p.S1OnOff_Npop,p.R1Off_Npop);
R1Off_S1OnOff_PSC_scale = (p.R1Off_S1OnOff_PSC_tauD/p.R1Off_S1OnOff_PSC_tauR)^(p.R1Off_S1OnOff_PSC_tauR/(p.R1Off_S1OnOff_PSC_tauD-p.R1Off_S1OnOff_PSC_tauR));
R1Off_Off_PSC_netcon = eye(p.Off_Npop,p.R1Off_Npop);
R1Off_Off_PSC_scale = (p.R1Off_Off_PSC_tauD/p.R1Off_Off_PSC_tauR)^(p.R1Off_Off_PSC_tauR/(p.R1Off_Off_PSC_tauD-p.R1Off_Off_PSC_tauR));
S1OnOff_Off_PSC_netcon = eye(p.Off_Npop,p.S1OnOff_Npop);
S1OnOff_Off_PSC_scale = (p.S1OnOff_Off_PSC_tauD/p.S1OnOff_Off_PSC_tauR)^(p.S1OnOff_Off_PSC_tauR/(p.S1OnOff_Off_PSC_tauD-p.S1OnOff_Off_PSC_tauR));
R2On_R1On_PSC_netcon = eye(p.R1On_Npop,p.R2On_Npop);
R2On_R1On_PSC_scale = (p.R2On_R1On_PSC_tauD/p.R2On_R1On_PSC_tauR)^(p.R2On_R1On_PSC_tauR/(p.R2On_R1On_PSC_tauD-p.R2On_R1On_PSC_tauR));
S2OnOff_R1On_PSC_netcon = eye(p.R1On_Npop,p.S2OnOff_Npop);
S2OnOff_R1On_PSC_scale = (p.S2OnOff_R1On_PSC_tauD/p.S2OnOff_R1On_PSC_tauR)^(p.S2OnOff_R1On_PSC_tauR/(p.S2OnOff_R1On_PSC_tauD-p.S2OnOff_R1On_PSC_tauR));
R2On_S2OnOff_PSC_netcon = eye(p.S2OnOff_Npop,p.R2On_Npop);
R2On_S2OnOff_PSC_scale = (p.R2On_S2OnOff_PSC_tauD/p.R2On_S2OnOff_PSC_tauR)^(p.R2On_S2OnOff_PSC_tauR/(p.R2On_S2OnOff_PSC_tauD-p.R2On_S2OnOff_PSC_tauR));
R2Off_S2OnOff_PSC_netcon = eye(p.S2OnOff_Npop,p.R2Off_Npop);
R2Off_S2OnOff_PSC_scale = (p.R2Off_S2OnOff_PSC_tauD/p.R2Off_S2OnOff_PSC_tauR)^(p.R2Off_S2OnOff_PSC_tauR/(p.R2Off_S2OnOff_PSC_tauD-p.R2Off_S2OnOff_PSC_tauR));
R2Off_R1Off_PSC_netcon = eye(p.R1Off_Npop,p.R2Off_Npop);
R2Off_R1Off_PSC_scale = (p.R2Off_R1Off_PSC_tauD/p.R2Off_R1Off_PSC_tauR)^(p.R2Off_R1Off_PSC_tauR/(p.R2Off_R1Off_PSC_tauD-p.R2Off_R1Off_PSC_tauR));
S2OnOff_R1Off_PSC_netcon = eye(p.R1Off_Npop,p.S2OnOff_Npop);
S2OnOff_R1Off_PSC_scale = (p.S2OnOff_R1Off_PSC_tauD/p.S2OnOff_R1Off_PSC_tauR)^(p.S2OnOff_R1Off_PSC_tauR/(p.S2OnOff_R1Off_PSC_tauD-p.S2OnOff_R1Off_PSC_tauR));
R2On_R2On_iNoise_V3_netcon =  eye(p.R2On_Npop,p.R2On_Npop);
R2On_R2On_iNoise_V3_token = genPoissonTimes(p.R2On_Npop,p.R2On_R2On_iNoise_V3_dt,p.R2On_R2On_iNoise_V3_FR,p.R2On_R2On_iNoise_V3_sigma,p.R2On_R2On_iNoise_V3_simlen);
R2On_R2On_iNoise_V3_scale =  (p.R2On_R2On_iNoise_V3_tauD_N/p.R2On_R2On_iNoise_V3_tauR_N)^(p.R2On_R2On_iNoise_V3_tauR_N/(p.R2On_R2On_iNoise_V3_tauD_N-p.R2On_R2On_iNoise_V3_tauR_N));

% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
t=0; k=1;

% STATE_VARIABLES:
On_V = zeros(nsamp,p.On_Npop);
On_V(1,:) =  p.On_E_L*ones(1,p.On_Npop);
On_g_ad = zeros(nsamp,p.On_Npop);
On_g_ad(1,:) =  zeros(1,p.On_Npop);
Off_V = zeros(nsamp,p.Off_Npop);
Off_V(1,:) =  p.Off_E_L*ones(1,p.Off_Npop);
Off_g_ad = zeros(nsamp,p.Off_Npop);
Off_g_ad(1,:) =  zeros(1,p.Off_Npop);
R1On_V = zeros(nsamp,p.R1On_Npop);
R1On_V(1,:) =  p.R1On_E_L*ones(1,p.R1On_Npop);
R1On_g_ad = zeros(nsamp,p.R1On_Npop);
R1On_g_ad(1,:) =  zeros(1,p.R1On_Npop);
R1Off_V = zeros(nsamp,p.R1Off_Npop);
R1Off_V(1,:) =  p.R1Off_E_L*ones(1,p.R1Off_Npop);
R1Off_g_ad = zeros(nsamp,p.R1Off_Npop);
R1Off_g_ad(1,:) =  zeros(1,p.R1Off_Npop);
S1OnOff_V = zeros(nsamp,p.S1OnOff_Npop);
S1OnOff_V(1,:) =  p.S1OnOff_E_L*ones(1,p.S1OnOff_Npop);
S1OnOff_g_ad = zeros(nsamp,p.S1OnOff_Npop);
S1OnOff_g_ad(1,:) =  zeros(1,p.S1OnOff_Npop);
R2On_V = zeros(nsamp,p.R2On_Npop);
R2On_V(1,:) =  p.R2On_E_L*ones(1,p.R2On_Npop);
R2On_g_ad = zeros(nsamp,p.R2On_Npop);
R2On_g_ad(1,:) =  zeros(1,p.R2On_Npop);
R2Off_V = zeros(nsamp,p.R2Off_Npop);
R2Off_V(1,:) =  p.R2Off_E_L*ones(1,p.R2Off_Npop);
R2Off_g_ad = zeros(nsamp,p.R2Off_Npop);
R2Off_g_ad(1,:) =  zeros(1,p.R2Off_Npop);
S2OnOff_V = zeros(nsamp,p.S2OnOff_Npop);
S2OnOff_V(1,:) =  p.S2OnOff_E_L*ones(1,p.S2OnOff_Npop);
S2OnOff_g_ad = zeros(nsamp,p.S2OnOff_Npop);
S2OnOff_g_ad(1,:) =  zeros(1,p.S2OnOff_Npop);
R1On_On_PSC_s = zeros(nsamp,p.On_Npop);
R1On_On_PSC_s(1,:) =  zeros(1,p.On_Npop);
R1On_On_PSC_x = zeros(nsamp,p.On_Npop);
R1On_On_PSC_x(1,:) =  zeros(1,p.On_Npop);
R1On_On_PSC_F = zeros(nsamp,p.On_Npop);
R1On_On_PSC_F(1,:) =  ones(1,p.On_Npop);
R1On_On_PSC_P = zeros(nsamp,p.On_Npop);
R1On_On_PSC_P(1,:) =  ones(1,p.On_Npop);
R1On_On_PSC_q = zeros(nsamp,p.On_Npop);
R1On_On_PSC_q(1,:) =  ones(1,p.On_Npop);
S1OnOff_On_PSC_s = zeros(nsamp,p.On_Npop);
S1OnOff_On_PSC_s(1,:) =  zeros(1,p.On_Npop);
S1OnOff_On_PSC_x = zeros(nsamp,p.On_Npop);
S1OnOff_On_PSC_x(1,:) =  zeros(1,p.On_Npop);
S1OnOff_On_PSC_F = zeros(nsamp,p.On_Npop);
S1OnOff_On_PSC_F(1,:) =  ones(1,p.On_Npop);
S1OnOff_On_PSC_P = zeros(nsamp,p.On_Npop);
S1OnOff_On_PSC_P(1,:) =  ones(1,p.On_Npop);
S1OnOff_On_PSC_q = zeros(nsamp,p.On_Npop);
S1OnOff_On_PSC_q(1,:) =  ones(1,p.On_Npop);
R1On_S1OnOff_PSC_s = zeros(nsamp,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_s(1,:) =  zeros(1,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_x = zeros(nsamp,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_x(1,:) =  zeros(1,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_F = zeros(nsamp,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_F(1,:) =  ones(1,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_P = zeros(nsamp,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_P(1,:) =  ones(1,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_q = zeros(nsamp,p.S1OnOff_Npop);
R1On_S1OnOff_PSC_q(1,:) =  ones(1,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_s = zeros(nsamp,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_s(1,:) =  zeros(1,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_x = zeros(nsamp,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_x(1,:) =  zeros(1,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_F = zeros(nsamp,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_F(1,:) =  ones(1,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_P = zeros(nsamp,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_P(1,:) =  ones(1,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_q = zeros(nsamp,p.S1OnOff_Npop);
R1Off_S1OnOff_PSC_q(1,:) =  ones(1,p.S1OnOff_Npop);
R1Off_Off_PSC_s = zeros(nsamp,p.Off_Npop);
R1Off_Off_PSC_s(1,:) =  zeros(1,p.Off_Npop);
R1Off_Off_PSC_x = zeros(nsamp,p.Off_Npop);
R1Off_Off_PSC_x(1,:) =  zeros(1,p.Off_Npop);
R1Off_Off_PSC_F = zeros(nsamp,p.Off_Npop);
R1Off_Off_PSC_F(1,:) =  ones(1,p.Off_Npop);
R1Off_Off_PSC_P = zeros(nsamp,p.Off_Npop);
R1Off_Off_PSC_P(1,:) =  ones(1,p.Off_Npop);
R1Off_Off_PSC_q = zeros(nsamp,p.Off_Npop);
R1Off_Off_PSC_q(1,:) =  ones(1,p.Off_Npop);
S1OnOff_Off_PSC_s = zeros(nsamp,p.Off_Npop);
S1OnOff_Off_PSC_s(1,:) =  zeros(1,p.Off_Npop);
S1OnOff_Off_PSC_x = zeros(nsamp,p.Off_Npop);
S1OnOff_Off_PSC_x(1,:) =  zeros(1,p.Off_Npop);
S1OnOff_Off_PSC_F = zeros(nsamp,p.Off_Npop);
S1OnOff_Off_PSC_F(1,:) =  ones(1,p.Off_Npop);
S1OnOff_Off_PSC_P = zeros(nsamp,p.Off_Npop);
S1OnOff_Off_PSC_P(1,:) =  ones(1,p.Off_Npop);
S1OnOff_Off_PSC_q = zeros(nsamp,p.Off_Npop);
S1OnOff_Off_PSC_q(1,:) =  ones(1,p.Off_Npop);
R2On_R1On_PSC_s = zeros(nsamp,p.R1On_Npop);
R2On_R1On_PSC_s(1,:) =  zeros(1,p.R1On_Npop);
R2On_R1On_PSC_x = zeros(nsamp,p.R1On_Npop);
R2On_R1On_PSC_x(1,:) =  zeros(1,p.R1On_Npop);
R2On_R1On_PSC_F = zeros(nsamp,p.R1On_Npop);
R2On_R1On_PSC_F(1,:) =  ones(1,p.R1On_Npop);
R2On_R1On_PSC_P = zeros(nsamp,p.R1On_Npop);
R2On_R1On_PSC_P(1,:) =  ones(1,p.R1On_Npop);
R2On_R1On_PSC_q = zeros(nsamp,p.R1On_Npop);
R2On_R1On_PSC_q(1,:) =  ones(1,p.R1On_Npop);
S2OnOff_R1On_PSC_s = zeros(nsamp,p.R1On_Npop);
S2OnOff_R1On_PSC_s(1,:) =  zeros(1,p.R1On_Npop);
S2OnOff_R1On_PSC_x = zeros(nsamp,p.R1On_Npop);
S2OnOff_R1On_PSC_x(1,:) =  zeros(1,p.R1On_Npop);
S2OnOff_R1On_PSC_F = zeros(nsamp,p.R1On_Npop);
S2OnOff_R1On_PSC_F(1,:) =  ones(1,p.R1On_Npop);
S2OnOff_R1On_PSC_P = zeros(nsamp,p.R1On_Npop);
S2OnOff_R1On_PSC_P(1,:) =  ones(1,p.R1On_Npop);
S2OnOff_R1On_PSC_q = zeros(nsamp,p.R1On_Npop);
S2OnOff_R1On_PSC_q(1,:) =  ones(1,p.R1On_Npop);
R2On_S2OnOff_PSC_s = zeros(nsamp,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_s(1,:) =  zeros(1,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_x = zeros(nsamp,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_x(1,:) =  zeros(1,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_F = zeros(nsamp,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_F(1,:) =  ones(1,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_P = zeros(nsamp,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_P(1,:) =  ones(1,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_q = zeros(nsamp,p.S2OnOff_Npop);
R2On_S2OnOff_PSC_q(1,:) =  ones(1,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_s = zeros(nsamp,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_s(1,:) =  zeros(1,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_x = zeros(nsamp,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_x(1,:) =  zeros(1,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_F = zeros(nsamp,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_F(1,:) =  ones(1,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_P = zeros(nsamp,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_P(1,:) =  ones(1,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_q = zeros(nsamp,p.S2OnOff_Npop);
R2Off_S2OnOff_PSC_q(1,:) =  ones(1,p.S2OnOff_Npop);
R2Off_R1Off_PSC_s = zeros(nsamp,p.R1Off_Npop);
R2Off_R1Off_PSC_s(1,:) =  zeros(1,p.R1Off_Npop);
R2Off_R1Off_PSC_x = zeros(nsamp,p.R1Off_Npop);
R2Off_R1Off_PSC_x(1,:) =  zeros(1,p.R1Off_Npop);
R2Off_R1Off_PSC_F = zeros(nsamp,p.R1Off_Npop);
R2Off_R1Off_PSC_F(1,:) =  ones(1,p.R1Off_Npop);
R2Off_R1Off_PSC_P = zeros(nsamp,p.R1Off_Npop);
R2Off_R1Off_PSC_P(1,:) =  ones(1,p.R1Off_Npop);
R2Off_R1Off_PSC_q = zeros(nsamp,p.R1Off_Npop);
R2Off_R1Off_PSC_q(1,:) =  ones(1,p.R1Off_Npop);
S2OnOff_R1Off_PSC_s = zeros(nsamp,p.R1Off_Npop);
S2OnOff_R1Off_PSC_s(1,:) =  zeros(1,p.R1Off_Npop);
S2OnOff_R1Off_PSC_x = zeros(nsamp,p.R1Off_Npop);
S2OnOff_R1Off_PSC_x(1,:) =  zeros(1,p.R1Off_Npop);
S2OnOff_R1Off_PSC_F = zeros(nsamp,p.R1Off_Npop);
S2OnOff_R1Off_PSC_F(1,:) =  ones(1,p.R1Off_Npop);
S2OnOff_R1Off_PSC_P = zeros(nsamp,p.R1Off_Npop);
S2OnOff_R1Off_PSC_P(1,:) =  ones(1,p.R1Off_Npop);
S2OnOff_R1Off_PSC_q = zeros(nsamp,p.R1Off_Npop);
S2OnOff_R1Off_PSC_q(1,:) =  ones(1,p.R1Off_Npop);
R2On_R2On_iNoise_V3_sn = zeros(nsamp,p.R2On_Npop);
R2On_R2On_iNoise_V3_sn(1,:) =  0 * ones(1,p.R2On_Npop);
R2On_R2On_iNoise_V3_xn = zeros(nsamp,p.R2On_Npop);
R2On_R2On_iNoise_V3_xn(1,:) =  0 * ones(1,p.R2On_Npop);

% MONITORS:
On_tspike = -1e32*ones(5,p.On_Npop);
On_buffer_index = ones(1,p.On_Npop);
On_V_spikes = zeros(nsamp,p.On_Npop);
Off_tspike = -1e32*ones(5,p.Off_Npop);
Off_buffer_index = ones(1,p.Off_Npop);
Off_V_spikes = zeros(nsamp,p.Off_Npop);
R1On_tspike = -1e32*ones(5,p.R1On_Npop);
R1On_buffer_index = ones(1,p.R1On_Npop);
R1On_V_spikes = zeros(nsamp,p.R1On_Npop);
R1Off_tspike = -1e32*ones(5,p.R1Off_Npop);
R1Off_buffer_index = ones(1,p.R1Off_Npop);
R1Off_V_spikes = zeros(nsamp,p.R1Off_Npop);
S1OnOff_tspike = -1e32*ones(5,p.S1OnOff_Npop);
S1OnOff_buffer_index = ones(1,p.S1OnOff_Npop);
S1OnOff_V_spikes = zeros(nsamp,p.S1OnOff_Npop);
R2On_tspike = -1e32*ones(5,p.R2On_Npop);
R2On_buffer_index = ones(1,p.R2On_Npop);
R2On_V_spikes = zeros(nsamp,p.R2On_Npop);
R2Off_tspike = -1e32*ones(5,p.R2Off_Npop);
R2Off_buffer_index = ones(1,p.R2Off_Npop);
R2Off_V_spikes = zeros(nsamp,p.R2Off_Npop);
S2OnOff_tspike = -1e32*ones(5,p.S2OnOff_Npop);
S2OnOff_buffer_index = ones(1,p.S2OnOff_Npop);
S2OnOff_V_spikes = zeros(nsamp,p.S2OnOff_Npop);
On_On_IC_iIC = zeros(nsamp,p.On_Npop);
On_On_IC_iIC(1,:)=p.On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(1,:)-p.On_On_IC_E_exc);
Off_Off_IC_iIC = zeros(nsamp,p.Off_Npop);
Off_Off_IC_iIC(1,:)=p.Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(1,:)-p.Off_Off_IC_E_exc);
R1On_On_PSC_syn = zeros(nsamp,p.R1On_Npop);
R1On_On_PSC_syn(1,:)=p.R1On_On_PSC_gSYN.*(R1On_On_PSC_s(1,:)*R1On_On_PSC_netcon).*(R1On_V(1,:)-p.R1On_On_PSC_ESYN);
S1OnOff_On_PSC_syn = zeros(nsamp,p.S1OnOff_Npop);
S1OnOff_On_PSC_syn(1,:)=p.S1OnOff_On_PSC_gSYN.*(S1OnOff_On_PSC_s(1,:)*S1OnOff_On_PSC_netcon).*(S1OnOff_V(1,:)-p.S1OnOff_On_PSC_ESYN);
R1On_S1OnOff_PSC_syn = zeros(nsamp,p.R1On_Npop);
R1On_S1OnOff_PSC_syn(1,:)=p.R1On_S1OnOff_PSC_gSYN.*(R1On_S1OnOff_PSC_s(1,:)*R1On_S1OnOff_PSC_netcon).*(R1On_V(1,:)-p.R1On_S1OnOff_PSC_ESYN);
R1Off_S1OnOff_PSC_syn = zeros(nsamp,p.R1Off_Npop);
R1Off_S1OnOff_PSC_syn(1,:)=p.R1Off_S1OnOff_PSC_gSYN.*(R1Off_S1OnOff_PSC_s(1,:)*R1Off_S1OnOff_PSC_netcon).*(R1Off_V(1,:)-p.R1Off_S1OnOff_PSC_ESYN);
R1Off_Off_PSC_syn = zeros(nsamp,p.R1Off_Npop);
R1Off_Off_PSC_syn(1,:)=p.R1Off_Off_PSC_gSYN.*(R1Off_Off_PSC_s(1,:)*R1Off_Off_PSC_netcon).*(R1Off_V(1,:)-p.R1Off_Off_PSC_ESYN);
S1OnOff_Off_PSC_syn = zeros(nsamp,p.S1OnOff_Npop);
S1OnOff_Off_PSC_syn(1,:)=p.S1OnOff_Off_PSC_gSYN.*(S1OnOff_Off_PSC_s(1,:)*S1OnOff_Off_PSC_netcon).*(S1OnOff_V(1,:)-p.S1OnOff_Off_PSC_ESYN);
R2On_R1On_PSC_syn = zeros(nsamp,p.R2On_Npop);
R2On_R1On_PSC_syn(1,:)=p.R2On_R1On_PSC_gSYN.*(R2On_R1On_PSC_s(1,:)*R2On_R1On_PSC_netcon).*(R2On_V(1,:)-p.R2On_R1On_PSC_ESYN);
S2OnOff_R1On_PSC_syn = zeros(nsamp,p.S2OnOff_Npop);
S2OnOff_R1On_PSC_syn(1,:)=p.S2OnOff_R1On_PSC_gSYN.*(S2OnOff_R1On_PSC_s(1,:)*S2OnOff_R1On_PSC_netcon).*(S2OnOff_V(1,:)-p.S2OnOff_R1On_PSC_ESYN);
R2On_S2OnOff_PSC_syn = zeros(nsamp,p.R2On_Npop);
R2On_S2OnOff_PSC_syn(1,:)=p.R2On_S2OnOff_PSC_gSYN.*(R2On_S2OnOff_PSC_s(1,:)*R2On_S2OnOff_PSC_netcon).*(R2On_V(1,:)-p.R2On_S2OnOff_PSC_ESYN);
R2Off_S2OnOff_PSC_syn = zeros(nsamp,p.R2Off_Npop);
R2Off_S2OnOff_PSC_syn(1,:)=p.R2Off_S2OnOff_PSC_gSYN.*(R2Off_S2OnOff_PSC_s(1,:)*R2Off_S2OnOff_PSC_netcon).*(R2Off_V(1,:)-p.R2Off_S2OnOff_PSC_ESYN);
R2Off_R1Off_PSC_syn = zeros(nsamp,p.R2Off_Npop);
R2Off_R1Off_PSC_syn(1,:)=p.R2Off_R1Off_PSC_gSYN.*(R2Off_R1Off_PSC_s(1,:)*R2Off_R1Off_PSC_netcon).*(R2Off_V(1,:)-p.R2Off_R1Off_PSC_ESYN);
S2OnOff_R1Off_PSC_syn = zeros(nsamp,p.S2OnOff_Npop);
S2OnOff_R1Off_PSC_syn(1,:)=p.S2OnOff_R1Off_PSC_gSYN.*(S2OnOff_R1Off_PSC_s(1,:)*S2OnOff_R1Off_PSC_netcon).*(S2OnOff_V(1,:)-p.S2OnOff_R1Off_PSC_ESYN);

% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  On_V_k1 = ( (p.On_E_L-On_V(n-1,:)) - On_R*On_g_ad(n-1,:).*(On_V(n-1,:)-p.On_E_k) - On_R*((((p.On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(n-1,:)-p.On_On_IC_E_exc))))) + On_R*p.On_Itonic.*On_Imask + On_R*p.On_noise.*randn(1,p.On_Npop) ) / On_tau;
  On_g_ad_k1 = -On_g_ad(n-1,:) / p.On_tau_ad;
  Off_V_k1 = ( (p.Off_E_L-Off_V(n-1,:)) - Off_R*Off_g_ad(n-1,:).*(Off_V(n-1,:)-p.Off_E_k) - Off_R*((((p.Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(n-1,:)-p.Off_Off_IC_E_exc))))) + Off_R*p.Off_Itonic.*Off_Imask + Off_R*p.Off_noise.*randn(1,p.Off_Npop) ) / Off_tau;
  Off_g_ad_k1 = -Off_g_ad(n-1,:) / p.Off_tau_ad;
  R1On_V_k1 = ( (p.R1On_E_L-R1On_V(n-1,:)) - R1On_R*R1On_g_ad(n-1,:).*(R1On_V(n-1,:)-p.R1On_E_k) - R1On_R*((((p.R1On_On_PSC_gSYN.*(R1On_On_PSC_s(n-1,:)*R1On_On_PSC_netcon).*(R1On_V(n-1,:)-p.R1On_On_PSC_ESYN))))+((((p.R1On_S1OnOff_PSC_gSYN.*(R1On_S1OnOff_PSC_s(n-1,:)*R1On_S1OnOff_PSC_netcon).*(R1On_V(n-1,:)-p.R1On_S1OnOff_PSC_ESYN)))))) + R1On_R*p.R1On_Itonic.*R1On_Imask + R1On_R*p.R1On_noise.*randn(1,p.R1On_Npop) ) / R1On_tau;
  R1On_g_ad_k1 = -R1On_g_ad(n-1,:) / p.R1On_tau_ad;
  R1Off_V_k1 = ( (p.R1Off_E_L-R1Off_V(n-1,:)) - R1Off_R*R1Off_g_ad(n-1,:).*(R1Off_V(n-1,:)-p.R1Off_E_k) - R1Off_R*((((p.R1Off_S1OnOff_PSC_gSYN.*(R1Off_S1OnOff_PSC_s(n-1,:)*R1Off_S1OnOff_PSC_netcon).*(R1Off_V(n-1,:)-p.R1Off_S1OnOff_PSC_ESYN))))+((((p.R1Off_Off_PSC_gSYN.*(R1Off_Off_PSC_s(n-1,:)*R1Off_Off_PSC_netcon).*(R1Off_V(n-1,:)-p.R1Off_Off_PSC_ESYN)))))) + R1Off_R*p.R1Off_Itonic.*R1Off_Imask + R1Off_R*p.R1Off_noise.*randn(1,p.R1Off_Npop) ) / R1Off_tau;
  R1Off_g_ad_k1 = -R1Off_g_ad(n-1,:) / p.R1Off_tau_ad;
  S1OnOff_V_k1 = ( (p.S1OnOff_E_L-S1OnOff_V(n-1,:)) - S1OnOff_R*S1OnOff_g_ad(n-1,:).*(S1OnOff_V(n-1,:)-p.S1OnOff_E_k) - S1OnOff_R*((((p.S1OnOff_On_PSC_gSYN.*(S1OnOff_On_PSC_s(n-1,:)*S1OnOff_On_PSC_netcon).*(S1OnOff_V(n-1,:)-p.S1OnOff_On_PSC_ESYN))))+((((p.S1OnOff_Off_PSC_gSYN.*(S1OnOff_Off_PSC_s(n-1,:)*S1OnOff_Off_PSC_netcon).*(S1OnOff_V(n-1,:)-p.S1OnOff_Off_PSC_ESYN)))))) + S1OnOff_R*p.S1OnOff_Itonic.*S1OnOff_Imask + S1OnOff_R*p.S1OnOff_noise.*randn(1,p.S1OnOff_Npop) ) / S1OnOff_tau;
  S1OnOff_g_ad_k1 = -S1OnOff_g_ad(n-1,:) / p.S1OnOff_tau_ad;
  R2On_V_k1 = ( (p.R2On_E_L-R2On_V(n-1,:)) - R2On_R*R2On_g_ad(n-1,:).*(R2On_V(n-1,:)-p.R2On_E_k) - R2On_R*((((p.R2On_R1On_PSC_gSYN.*(R2On_R1On_PSC_s(n-1,:)*R2On_R1On_PSC_netcon).*(R2On_V(n-1,:)-p.R2On_R1On_PSC_ESYN))))+((((p.R2On_S2OnOff_PSC_gSYN.*(R2On_S2OnOff_PSC_s(n-1,:)*R2On_S2OnOff_PSC_netcon).*(R2On_V(n-1,:)-p.R2On_S2OnOff_PSC_ESYN))))+((((p.R2On_R2On_iNoise_V3_nSYN.*(R2On_R2On_iNoise_V3_sn(n-1,:)*R2On_R2On_iNoise_V3_netcon).*(R2On_V(n-1,:)-p.R2On_R2On_iNoise_V3_E_exc))))))) + R2On_R*p.R2On_Itonic.*R2On_Imask + R2On_R*p.R2On_noise.*randn(1,p.R2On_Npop) ) / R2On_tau;
  R2On_g_ad_k1 = -R2On_g_ad(n-1,:) / p.R2On_tau_ad;
  R2Off_V_k1 = ( (p.R2Off_E_L-R2Off_V(n-1,:)) - R2Off_R*R2Off_g_ad(n-1,:).*(R2Off_V(n-1,:)-p.R2Off_E_k) - R2Off_R*((((p.R2Off_S2OnOff_PSC_gSYN.*(R2Off_S2OnOff_PSC_s(n-1,:)*R2Off_S2OnOff_PSC_netcon).*(R2Off_V(n-1,:)-p.R2Off_S2OnOff_PSC_ESYN))))+((((p.R2Off_R1Off_PSC_gSYN.*(R2Off_R1Off_PSC_s(n-1,:)*R2Off_R1Off_PSC_netcon).*(R2Off_V(n-1,:)-p.R2Off_R1Off_PSC_ESYN)))))) + R2Off_R*p.R2Off_Itonic.*R2Off_Imask + R2Off_R*p.R2Off_noise.*randn(1,p.R2Off_Npop) ) / R2Off_tau;
  R2Off_g_ad_k1 = -R2Off_g_ad(n-1,:) / p.R2Off_tau_ad;
  S2OnOff_V_k1 = ( (p.S2OnOff_E_L-S2OnOff_V(n-1,:)) - S2OnOff_R*S2OnOff_g_ad(n-1,:).*(S2OnOff_V(n-1,:)-p.S2OnOff_E_k) - S2OnOff_R*((((p.S2OnOff_R1On_PSC_gSYN.*(S2OnOff_R1On_PSC_s(n-1,:)*S2OnOff_R1On_PSC_netcon).*(S2OnOff_V(n-1,:)-p.S2OnOff_R1On_PSC_ESYN))))+((((p.S2OnOff_R1Off_PSC_gSYN.*(S2OnOff_R1Off_PSC_s(n-1,:)*S2OnOff_R1Off_PSC_netcon).*(S2OnOff_V(n-1,:)-p.S2OnOff_R1Off_PSC_ESYN)))))) + S2OnOff_R*p.S2OnOff_Itonic.*S2OnOff_Imask + S2OnOff_R*p.S2OnOff_noise.*randn(1,p.S2OnOff_Npop) ) / S2OnOff_tau;
  S2OnOff_g_ad_k1 = -S2OnOff_g_ad(n-1,:) / p.S2OnOff_tau_ad;
  R1On_On_PSC_s_k1 = ( R1On_On_PSC_scale * R1On_On_PSC_x(n-1,:) - R1On_On_PSC_s(n-1,:) )/p.R1On_On_PSC_tauR;
  R1On_On_PSC_x_k1 = -R1On_On_PSC_x(n-1,:)/p.R1On_On_PSC_tauD;
  R1On_On_PSC_F_k1 = (1 - R1On_On_PSC_F(n-1,:))/p.R1On_On_PSC_tauF;
  R1On_On_PSC_P_k1 = (1 - R1On_On_PSC_P(n-1,:))/p.R1On_On_PSC_tauP;
  R1On_On_PSC_q_k1 = 0;
  S1OnOff_On_PSC_s_k1 = ( S1OnOff_On_PSC_scale * S1OnOff_On_PSC_x(n-1,:) - S1OnOff_On_PSC_s(n-1,:) )/p.S1OnOff_On_PSC_tauR;
  S1OnOff_On_PSC_x_k1 = -S1OnOff_On_PSC_x(n-1,:)/p.S1OnOff_On_PSC_tauD;
  S1OnOff_On_PSC_F_k1 = (1 - S1OnOff_On_PSC_F(n-1,:))/p.S1OnOff_On_PSC_tauF;
  S1OnOff_On_PSC_P_k1 = (1 - S1OnOff_On_PSC_P(n-1,:))/p.S1OnOff_On_PSC_tauP;
  S1OnOff_On_PSC_q_k1 = 0;
  R1On_S1OnOff_PSC_s_k1 = ( R1On_S1OnOff_PSC_scale * R1On_S1OnOff_PSC_x(n-1,:) - R1On_S1OnOff_PSC_s(n-1,:) )/p.R1On_S1OnOff_PSC_tauR;
  R1On_S1OnOff_PSC_x_k1 = -R1On_S1OnOff_PSC_x(n-1,:)/p.R1On_S1OnOff_PSC_tauD;
  R1On_S1OnOff_PSC_F_k1 = (1 - R1On_S1OnOff_PSC_F(n-1,:))/p.R1On_S1OnOff_PSC_tauF;
  R1On_S1OnOff_PSC_P_k1 = (1 - R1On_S1OnOff_PSC_P(n-1,:))/p.R1On_S1OnOff_PSC_tauP;
  R1On_S1OnOff_PSC_q_k1 = 0;
  R1Off_S1OnOff_PSC_s_k1 = ( R1Off_S1OnOff_PSC_scale * R1Off_S1OnOff_PSC_x(n-1,:) - R1Off_S1OnOff_PSC_s(n-1,:) )/p.R1Off_S1OnOff_PSC_tauR;
  R1Off_S1OnOff_PSC_x_k1 = -R1Off_S1OnOff_PSC_x(n-1,:)/p.R1Off_S1OnOff_PSC_tauD;
  R1Off_S1OnOff_PSC_F_k1 = (1 - R1Off_S1OnOff_PSC_F(n-1,:))/p.R1Off_S1OnOff_PSC_tauF;
  R1Off_S1OnOff_PSC_P_k1 = (1 - R1Off_S1OnOff_PSC_P(n-1,:))/p.R1Off_S1OnOff_PSC_tauP;
  R1Off_S1OnOff_PSC_q_k1 = 0;
  R1Off_Off_PSC_s_k1 = ( R1Off_Off_PSC_scale * R1Off_Off_PSC_x(n-1,:) - R1Off_Off_PSC_s(n-1,:) )/p.R1Off_Off_PSC_tauR;
  R1Off_Off_PSC_x_k1 = -R1Off_Off_PSC_x(n-1,:)/p.R1Off_Off_PSC_tauD;
  R1Off_Off_PSC_F_k1 = (1 - R1Off_Off_PSC_F(n-1,:))/p.R1Off_Off_PSC_tauF;
  R1Off_Off_PSC_P_k1 = (1 - R1Off_Off_PSC_P(n-1,:))/p.R1Off_Off_PSC_tauP;
  R1Off_Off_PSC_q_k1 = 0;
  S1OnOff_Off_PSC_s_k1 = ( S1OnOff_Off_PSC_scale * S1OnOff_Off_PSC_x(n-1,:) - S1OnOff_Off_PSC_s(n-1,:) )/p.S1OnOff_Off_PSC_tauR;
  S1OnOff_Off_PSC_x_k1 = -S1OnOff_Off_PSC_x(n-1,:)/p.S1OnOff_Off_PSC_tauD;
  S1OnOff_Off_PSC_F_k1 = (1 - S1OnOff_Off_PSC_F(n-1,:))/p.S1OnOff_Off_PSC_tauF;
  S1OnOff_Off_PSC_P_k1 = (1 - S1OnOff_Off_PSC_P(n-1,:))/p.S1OnOff_Off_PSC_tauP;
  S1OnOff_Off_PSC_q_k1 = 0;
  R2On_R1On_PSC_s_k1 = ( R2On_R1On_PSC_scale * R2On_R1On_PSC_x(n-1,:) - R2On_R1On_PSC_s(n-1,:) )/p.R2On_R1On_PSC_tauR;
  R2On_R1On_PSC_x_k1 = -R2On_R1On_PSC_x(n-1,:)/p.R2On_R1On_PSC_tauD;
  R2On_R1On_PSC_F_k1 = (1 - R2On_R1On_PSC_F(n-1,:))/p.R2On_R1On_PSC_tauF;
  R2On_R1On_PSC_P_k1 = (1 - R2On_R1On_PSC_P(n-1,:))/p.R2On_R1On_PSC_tauP;
  R2On_R1On_PSC_q_k1 = 0;
  S2OnOff_R1On_PSC_s_k1 = ( S2OnOff_R1On_PSC_scale * S2OnOff_R1On_PSC_x(n-1,:) - S2OnOff_R1On_PSC_s(n-1,:) )/p.S2OnOff_R1On_PSC_tauR;
  S2OnOff_R1On_PSC_x_k1 = -S2OnOff_R1On_PSC_x(n-1,:)/p.S2OnOff_R1On_PSC_tauD;
  S2OnOff_R1On_PSC_F_k1 = (1 - S2OnOff_R1On_PSC_F(n-1,:))/p.S2OnOff_R1On_PSC_tauF;
  S2OnOff_R1On_PSC_P_k1 = (1 - S2OnOff_R1On_PSC_P(n-1,:))/p.S2OnOff_R1On_PSC_tauP;
  S2OnOff_R1On_PSC_q_k1 = 0;
  R2On_S2OnOff_PSC_s_k1 = ( R2On_S2OnOff_PSC_scale * R2On_S2OnOff_PSC_x(n-1,:) - R2On_S2OnOff_PSC_s(n-1,:) )/p.R2On_S2OnOff_PSC_tauR;
  R2On_S2OnOff_PSC_x_k1 = -R2On_S2OnOff_PSC_x(n-1,:)/p.R2On_S2OnOff_PSC_tauD;
  R2On_S2OnOff_PSC_F_k1 = (1 - R2On_S2OnOff_PSC_F(n-1,:))/p.R2On_S2OnOff_PSC_tauF;
  R2On_S2OnOff_PSC_P_k1 = (1 - R2On_S2OnOff_PSC_P(n-1,:))/p.R2On_S2OnOff_PSC_tauP;
  R2On_S2OnOff_PSC_q_k1 = 0;
  R2Off_S2OnOff_PSC_s_k1 = ( R2Off_S2OnOff_PSC_scale * R2Off_S2OnOff_PSC_x(n-1,:) - R2Off_S2OnOff_PSC_s(n-1,:) )/p.R2Off_S2OnOff_PSC_tauR;
  R2Off_S2OnOff_PSC_x_k1 = -R2Off_S2OnOff_PSC_x(n-1,:)/p.R2Off_S2OnOff_PSC_tauD;
  R2Off_S2OnOff_PSC_F_k1 = (1 - R2Off_S2OnOff_PSC_F(n-1,:))/p.R2Off_S2OnOff_PSC_tauF;
  R2Off_S2OnOff_PSC_P_k1 = (1 - R2Off_S2OnOff_PSC_P(n-1,:))/p.R2Off_S2OnOff_PSC_tauP;
  R2Off_S2OnOff_PSC_q_k1 = 0;
  R2Off_R1Off_PSC_s_k1 = ( R2Off_R1Off_PSC_scale * R2Off_R1Off_PSC_x(n-1,:) - R2Off_R1Off_PSC_s(n-1,:) )/p.R2Off_R1Off_PSC_tauR;
  R2Off_R1Off_PSC_x_k1 = -R2Off_R1Off_PSC_x(n-1,:)/p.R2Off_R1Off_PSC_tauD;
  R2Off_R1Off_PSC_F_k1 = (1 - R2Off_R1Off_PSC_F(n-1,:))/p.R2Off_R1Off_PSC_tauF;
  R2Off_R1Off_PSC_P_k1 = (1 - R2Off_R1Off_PSC_P(n-1,:))/p.R2Off_R1Off_PSC_tauP;
  R2Off_R1Off_PSC_q_k1 = 0;
  S2OnOff_R1Off_PSC_s_k1 = ( S2OnOff_R1Off_PSC_scale * S2OnOff_R1Off_PSC_x(n-1,:) - S2OnOff_R1Off_PSC_s(n-1,:) )/p.S2OnOff_R1Off_PSC_tauR;
  S2OnOff_R1Off_PSC_x_k1 = -S2OnOff_R1Off_PSC_x(n-1,:)/p.S2OnOff_R1Off_PSC_tauD;
  S2OnOff_R1Off_PSC_F_k1 = (1 - S2OnOff_R1Off_PSC_F(n-1,:))/p.S2OnOff_R1Off_PSC_tauF;
  S2OnOff_R1Off_PSC_P_k1 = (1 - S2OnOff_R1Off_PSC_P(n-1,:))/p.S2OnOff_R1Off_PSC_tauP;
  S2OnOff_R1Off_PSC_q_k1 = 0;
  R2On_R2On_iNoise_V3_sn_k1 = ( R2On_R2On_iNoise_V3_scale * R2On_R2On_iNoise_V3_xn(n-1,:) - R2On_R2On_iNoise_V3_sn(n-1,:) )/p.R2On_R2On_iNoise_V3_tauR_N;
  R2On_R2On_iNoise_V3_xn_k1 = -R2On_R2On_iNoise_V3_xn(n-1,:)/p.R2On_R2On_iNoise_V3_tauD_N + R2On_R2On_iNoise_V3_token(k,:)/p.R2On_R2On_iNoise_V3_dt;

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  On_V(n,:) = On_V(n-1,:)+dt*On_V_k1;
  On_g_ad(n,:) = On_g_ad(n-1,:)+dt*On_g_ad_k1;
  Off_V(n,:) = Off_V(n-1,:)+dt*Off_V_k1;
  Off_g_ad(n,:) = Off_g_ad(n-1,:)+dt*Off_g_ad_k1;
  R1On_V(n,:) = R1On_V(n-1,:)+dt*R1On_V_k1;
  R1On_g_ad(n,:) = R1On_g_ad(n-1,:)+dt*R1On_g_ad_k1;
  R1Off_V(n,:) = R1Off_V(n-1,:)+dt*R1Off_V_k1;
  R1Off_g_ad(n,:) = R1Off_g_ad(n-1,:)+dt*R1Off_g_ad_k1;
  S1OnOff_V(n,:) = S1OnOff_V(n-1,:)+dt*S1OnOff_V_k1;
  S1OnOff_g_ad(n,:) = S1OnOff_g_ad(n-1,:)+dt*S1OnOff_g_ad_k1;
  R2On_V(n,:) = R2On_V(n-1,:)+dt*R2On_V_k1;
  R2On_g_ad(n,:) = R2On_g_ad(n-1,:)+dt*R2On_g_ad_k1;
  R2Off_V(n,:) = R2Off_V(n-1,:)+dt*R2Off_V_k1;
  R2Off_g_ad(n,:) = R2Off_g_ad(n-1,:)+dt*R2Off_g_ad_k1;
  S2OnOff_V(n,:) = S2OnOff_V(n-1,:)+dt*S2OnOff_V_k1;
  S2OnOff_g_ad(n,:) = S2OnOff_g_ad(n-1,:)+dt*S2OnOff_g_ad_k1;
  R1On_On_PSC_s(n,:) = R1On_On_PSC_s(n-1,:)+dt*R1On_On_PSC_s_k1;
  R1On_On_PSC_x(n,:) = R1On_On_PSC_x(n-1,:)+dt*R1On_On_PSC_x_k1;
  R1On_On_PSC_F(n,:) = R1On_On_PSC_F(n-1,:)+dt*R1On_On_PSC_F_k1;
  R1On_On_PSC_P(n,:) = R1On_On_PSC_P(n-1,:)+dt*R1On_On_PSC_P_k1;
  R1On_On_PSC_q(n,:) = R1On_On_PSC_q(n-1,:)+dt*R1On_On_PSC_q_k1;
  S1OnOff_On_PSC_s(n,:) = S1OnOff_On_PSC_s(n-1,:)+dt*S1OnOff_On_PSC_s_k1;
  S1OnOff_On_PSC_x(n,:) = S1OnOff_On_PSC_x(n-1,:)+dt*S1OnOff_On_PSC_x_k1;
  S1OnOff_On_PSC_F(n,:) = S1OnOff_On_PSC_F(n-1,:)+dt*S1OnOff_On_PSC_F_k1;
  S1OnOff_On_PSC_P(n,:) = S1OnOff_On_PSC_P(n-1,:)+dt*S1OnOff_On_PSC_P_k1;
  S1OnOff_On_PSC_q(n,:) = S1OnOff_On_PSC_q(n-1,:)+dt*S1OnOff_On_PSC_q_k1;
  R1On_S1OnOff_PSC_s(n,:) = R1On_S1OnOff_PSC_s(n-1,:)+dt*R1On_S1OnOff_PSC_s_k1;
  R1On_S1OnOff_PSC_x(n,:) = R1On_S1OnOff_PSC_x(n-1,:)+dt*R1On_S1OnOff_PSC_x_k1;
  R1On_S1OnOff_PSC_F(n,:) = R1On_S1OnOff_PSC_F(n-1,:)+dt*R1On_S1OnOff_PSC_F_k1;
  R1On_S1OnOff_PSC_P(n,:) = R1On_S1OnOff_PSC_P(n-1,:)+dt*R1On_S1OnOff_PSC_P_k1;
  R1On_S1OnOff_PSC_q(n,:) = R1On_S1OnOff_PSC_q(n-1,:)+dt*R1On_S1OnOff_PSC_q_k1;
  R1Off_S1OnOff_PSC_s(n,:) = R1Off_S1OnOff_PSC_s(n-1,:)+dt*R1Off_S1OnOff_PSC_s_k1;
  R1Off_S1OnOff_PSC_x(n,:) = R1Off_S1OnOff_PSC_x(n-1,:)+dt*R1Off_S1OnOff_PSC_x_k1;
  R1Off_S1OnOff_PSC_F(n,:) = R1Off_S1OnOff_PSC_F(n-1,:)+dt*R1Off_S1OnOff_PSC_F_k1;
  R1Off_S1OnOff_PSC_P(n,:) = R1Off_S1OnOff_PSC_P(n-1,:)+dt*R1Off_S1OnOff_PSC_P_k1;
  R1Off_S1OnOff_PSC_q(n,:) = R1Off_S1OnOff_PSC_q(n-1,:)+dt*R1Off_S1OnOff_PSC_q_k1;
  R1Off_Off_PSC_s(n,:) = R1Off_Off_PSC_s(n-1,:)+dt*R1Off_Off_PSC_s_k1;
  R1Off_Off_PSC_x(n,:) = R1Off_Off_PSC_x(n-1,:)+dt*R1Off_Off_PSC_x_k1;
  R1Off_Off_PSC_F(n,:) = R1Off_Off_PSC_F(n-1,:)+dt*R1Off_Off_PSC_F_k1;
  R1Off_Off_PSC_P(n,:) = R1Off_Off_PSC_P(n-1,:)+dt*R1Off_Off_PSC_P_k1;
  R1Off_Off_PSC_q(n,:) = R1Off_Off_PSC_q(n-1,:)+dt*R1Off_Off_PSC_q_k1;
  S1OnOff_Off_PSC_s(n,:) = S1OnOff_Off_PSC_s(n-1,:)+dt*S1OnOff_Off_PSC_s_k1;
  S1OnOff_Off_PSC_x(n,:) = S1OnOff_Off_PSC_x(n-1,:)+dt*S1OnOff_Off_PSC_x_k1;
  S1OnOff_Off_PSC_F(n,:) = S1OnOff_Off_PSC_F(n-1,:)+dt*S1OnOff_Off_PSC_F_k1;
  S1OnOff_Off_PSC_P(n,:) = S1OnOff_Off_PSC_P(n-1,:)+dt*S1OnOff_Off_PSC_P_k1;
  S1OnOff_Off_PSC_q(n,:) = S1OnOff_Off_PSC_q(n-1,:)+dt*S1OnOff_Off_PSC_q_k1;
  R2On_R1On_PSC_s(n,:) = R2On_R1On_PSC_s(n-1,:)+dt*R2On_R1On_PSC_s_k1;
  R2On_R1On_PSC_x(n,:) = R2On_R1On_PSC_x(n-1,:)+dt*R2On_R1On_PSC_x_k1;
  R2On_R1On_PSC_F(n,:) = R2On_R1On_PSC_F(n-1,:)+dt*R2On_R1On_PSC_F_k1;
  R2On_R1On_PSC_P(n,:) = R2On_R1On_PSC_P(n-1,:)+dt*R2On_R1On_PSC_P_k1;
  R2On_R1On_PSC_q(n,:) = R2On_R1On_PSC_q(n-1,:)+dt*R2On_R1On_PSC_q_k1;
  S2OnOff_R1On_PSC_s(n,:) = S2OnOff_R1On_PSC_s(n-1,:)+dt*S2OnOff_R1On_PSC_s_k1;
  S2OnOff_R1On_PSC_x(n,:) = S2OnOff_R1On_PSC_x(n-1,:)+dt*S2OnOff_R1On_PSC_x_k1;
  S2OnOff_R1On_PSC_F(n,:) = S2OnOff_R1On_PSC_F(n-1,:)+dt*S2OnOff_R1On_PSC_F_k1;
  S2OnOff_R1On_PSC_P(n,:) = S2OnOff_R1On_PSC_P(n-1,:)+dt*S2OnOff_R1On_PSC_P_k1;
  S2OnOff_R1On_PSC_q(n,:) = S2OnOff_R1On_PSC_q(n-1,:)+dt*S2OnOff_R1On_PSC_q_k1;
  R2On_S2OnOff_PSC_s(n,:) = R2On_S2OnOff_PSC_s(n-1,:)+dt*R2On_S2OnOff_PSC_s_k1;
  R2On_S2OnOff_PSC_x(n,:) = R2On_S2OnOff_PSC_x(n-1,:)+dt*R2On_S2OnOff_PSC_x_k1;
  R2On_S2OnOff_PSC_F(n,:) = R2On_S2OnOff_PSC_F(n-1,:)+dt*R2On_S2OnOff_PSC_F_k1;
  R2On_S2OnOff_PSC_P(n,:) = R2On_S2OnOff_PSC_P(n-1,:)+dt*R2On_S2OnOff_PSC_P_k1;
  R2On_S2OnOff_PSC_q(n,:) = R2On_S2OnOff_PSC_q(n-1,:)+dt*R2On_S2OnOff_PSC_q_k1;
  R2Off_S2OnOff_PSC_s(n,:) = R2Off_S2OnOff_PSC_s(n-1,:)+dt*R2Off_S2OnOff_PSC_s_k1;
  R2Off_S2OnOff_PSC_x(n,:) = R2Off_S2OnOff_PSC_x(n-1,:)+dt*R2Off_S2OnOff_PSC_x_k1;
  R2Off_S2OnOff_PSC_F(n,:) = R2Off_S2OnOff_PSC_F(n-1,:)+dt*R2Off_S2OnOff_PSC_F_k1;
  R2Off_S2OnOff_PSC_P(n,:) = R2Off_S2OnOff_PSC_P(n-1,:)+dt*R2Off_S2OnOff_PSC_P_k1;
  R2Off_S2OnOff_PSC_q(n,:) = R2Off_S2OnOff_PSC_q(n-1,:)+dt*R2Off_S2OnOff_PSC_q_k1;
  R2Off_R1Off_PSC_s(n,:) = R2Off_R1Off_PSC_s(n-1,:)+dt*R2Off_R1Off_PSC_s_k1;
  R2Off_R1Off_PSC_x(n,:) = R2Off_R1Off_PSC_x(n-1,:)+dt*R2Off_R1Off_PSC_x_k1;
  R2Off_R1Off_PSC_F(n,:) = R2Off_R1Off_PSC_F(n-1,:)+dt*R2Off_R1Off_PSC_F_k1;
  R2Off_R1Off_PSC_P(n,:) = R2Off_R1Off_PSC_P(n-1,:)+dt*R2Off_R1Off_PSC_P_k1;
  R2Off_R1Off_PSC_q(n,:) = R2Off_R1Off_PSC_q(n-1,:)+dt*R2Off_R1Off_PSC_q_k1;
  S2OnOff_R1Off_PSC_s(n,:) = S2OnOff_R1Off_PSC_s(n-1,:)+dt*S2OnOff_R1Off_PSC_s_k1;
  S2OnOff_R1Off_PSC_x(n,:) = S2OnOff_R1Off_PSC_x(n-1,:)+dt*S2OnOff_R1Off_PSC_x_k1;
  S2OnOff_R1Off_PSC_F(n,:) = S2OnOff_R1Off_PSC_F(n-1,:)+dt*S2OnOff_R1Off_PSC_F_k1;
  S2OnOff_R1Off_PSC_P(n,:) = S2OnOff_R1Off_PSC_P(n-1,:)+dt*S2OnOff_R1Off_PSC_P_k1;
  S2OnOff_R1Off_PSC_q(n,:) = S2OnOff_R1Off_PSC_q(n-1,:)+dt*S2OnOff_R1Off_PSC_q_k1;
  R2On_R2On_iNoise_V3_sn(n,:) = R2On_R2On_iNoise_V3_sn(n-1,:)+dt*R2On_R2On_iNoise_V3_sn_k1;
  R2On_R2On_iNoise_V3_xn(n,:) = R2On_R2On_iNoise_V3_xn(n-1,:)+dt*R2On_R2On_iNoise_V3_xn_k1;

  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=any(S2OnOff_V(n,:)>=p.S2OnOff_V_thresh&S2OnOff_V(n-1,:)<p.S2OnOff_V_thresh);
  conditional_indx=(S2OnOff_V(n,:)>=p.S2OnOff_V_thresh&S2OnOff_V(n-1,:)<p.S2OnOff_V_thresh);
  if conditional_test, S2OnOff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); S2OnOff_tspike(S2OnOff_buffer_index(i),i)=t; S2OnOff_buffer_index(i)=mod(-1+(S2OnOff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(R2Off_V(n,:)>=p.R2Off_V_thresh&R2Off_V(n-1,:)<p.R2Off_V_thresh);
  conditional_indx=(R2Off_V(n,:)>=p.R2Off_V_thresh&R2Off_V(n-1,:)<p.R2Off_V_thresh);
  if conditional_test, R2Off_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); R2Off_tspike(R2Off_buffer_index(i),i)=t; R2Off_buffer_index(i)=mod(-1+(R2Off_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(R2On_V(n,:)>=p.R2On_V_thresh&R2On_V(n-1,:)<p.R2On_V_thresh);
  conditional_indx=(R2On_V(n,:)>=p.R2On_V_thresh&R2On_V(n-1,:)<p.R2On_V_thresh);
  if conditional_test, R2On_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); R2On_tspike(R2On_buffer_index(i),i)=t; R2On_buffer_index(i)=mod(-1+(R2On_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(S1OnOff_V(n,:)>=p.S1OnOff_V_thresh&S1OnOff_V(n-1,:)<p.S1OnOff_V_thresh);
  conditional_indx=(S1OnOff_V(n,:)>=p.S1OnOff_V_thresh&S1OnOff_V(n-1,:)<p.S1OnOff_V_thresh);
  if conditional_test, S1OnOff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); S1OnOff_tspike(S1OnOff_buffer_index(i),i)=t; S1OnOff_buffer_index(i)=mod(-1+(S1OnOff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(R1Off_V(n,:)>=p.R1Off_V_thresh&R1Off_V(n-1,:)<p.R1Off_V_thresh);
  conditional_indx=(R1Off_V(n,:)>=p.R1Off_V_thresh&R1Off_V(n-1,:)<p.R1Off_V_thresh);
  if conditional_test, R1Off_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); R1Off_tspike(R1Off_buffer_index(i),i)=t; R1Off_buffer_index(i)=mod(-1+(R1Off_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(R1On_V(n,:)>=p.R1On_V_thresh&R1On_V(n-1,:)<p.R1On_V_thresh);
  conditional_indx=(R1On_V(n,:)>=p.R1On_V_thresh&R1On_V(n-1,:)<p.R1On_V_thresh);
  if conditional_test, R1On_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); R1On_tspike(R1On_buffer_index(i),i)=t; R1On_buffer_index(i)=mod(-1+(R1On_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(Off_V(n,:)>=p.Off_V_thresh&Off_V(n-1,:)<p.Off_V_thresh);
  conditional_indx=(Off_V(n,:)>=p.Off_V_thresh&Off_V(n-1,:)<p.Off_V_thresh);
  if conditional_test, Off_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); Off_tspike(Off_buffer_index(i),i)=t; Off_buffer_index(i)=mod(-1+(Off_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(On_V(n,:)>=p.On_V_thresh&On_V(n-1,:)<p.On_V_thresh);
  conditional_indx=(On_V(n,:)>=p.On_V_thresh&On_V(n-1,:)<p.On_V_thresh);
  if conditional_test, On_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); On_tspike(On_buffer_index(i),i)=t; On_buffer_index(i)=mod(-1+(On_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(On_V(n,:) > p.On_V_thresh);
  conditional_indx=(On_V(n,:) > p.On_V_thresh);
  if conditional_test, On_V(n,conditional_indx) = p.On_V_reset; On_g_ad(n,conditional_indx) = On_g_ad(n,conditional_indx) + p.On_g_inc; end
  conditional_test=any(any(t<=On_tspike+p.On_t_ref,1));
  conditional_indx=(any(t<=On_tspike+p.On_t_ref,1));
  if conditional_test, On_V(n,conditional_indx) = p.On_V_reset; end
  conditional_test=any(Off_V(n,:) > p.Off_V_thresh);
  conditional_indx=(Off_V(n,:) > p.Off_V_thresh);
  if conditional_test, Off_V(n,conditional_indx) = p.Off_V_reset; Off_g_ad(n,conditional_indx) = Off_g_ad(n,conditional_indx) + p.Off_g_inc; end
  conditional_test=any(any(t<=Off_tspike+p.Off_t_ref,1));
  conditional_indx=(any(t<=Off_tspike+p.Off_t_ref,1));
  if conditional_test, Off_V(n,conditional_indx) = p.Off_V_reset; end
  conditional_test=any(R1On_V(n,:) > p.R1On_V_thresh);
  conditional_indx=(R1On_V(n,:) > p.R1On_V_thresh);
  if conditional_test, R1On_V(n,conditional_indx) = p.R1On_V_reset; R1On_g_ad(n,conditional_indx) = R1On_g_ad(n,conditional_indx) + p.R1On_g_inc; end
  conditional_test=any(any(t<=R1On_tspike+p.R1On_t_ref,1));
  conditional_indx=(any(t<=R1On_tspike+p.R1On_t_ref,1));
  if conditional_test, R1On_V(n,conditional_indx) = p.R1On_V_reset; end
  conditional_test=any(R1Off_V(n,:) > p.R1Off_V_thresh);
  conditional_indx=(R1Off_V(n,:) > p.R1Off_V_thresh);
  if conditional_test, R1Off_V(n,conditional_indx) = p.R1Off_V_reset; R1Off_g_ad(n,conditional_indx) = R1Off_g_ad(n,conditional_indx) + p.R1Off_g_inc; end
  conditional_test=any(any(t<=R1Off_tspike+p.R1Off_t_ref,1));
  conditional_indx=(any(t<=R1Off_tspike+p.R1Off_t_ref,1));
  if conditional_test, R1Off_V(n,conditional_indx) = p.R1Off_V_reset; end
  conditional_test=any(S1OnOff_V(n,:) > p.S1OnOff_V_thresh);
  conditional_indx=(S1OnOff_V(n,:) > p.S1OnOff_V_thresh);
  if conditional_test, S1OnOff_V(n,conditional_indx) = p.S1OnOff_V_reset; S1OnOff_g_ad(n,conditional_indx) = S1OnOff_g_ad(n,conditional_indx) + p.S1OnOff_g_inc; end
  conditional_test=any(any(t<=S1OnOff_tspike+p.S1OnOff_t_ref,1));
  conditional_indx=(any(t<=S1OnOff_tspike+p.S1OnOff_t_ref,1));
  if conditional_test, S1OnOff_V(n,conditional_indx) = p.S1OnOff_V_reset; end
  conditional_test=any(R2On_V(n,:) > p.R2On_V_thresh);
  conditional_indx=(R2On_V(n,:) > p.R2On_V_thresh);
  if conditional_test, R2On_V(n,conditional_indx) = p.R2On_V_reset; R2On_g_ad(n,conditional_indx) = R2On_g_ad(n,conditional_indx) + p.R2On_g_inc; end
  conditional_test=any(any(t<=R2On_tspike+p.R2On_t_ref,1));
  conditional_indx=(any(t<=R2On_tspike+p.R2On_t_ref,1));
  if conditional_test, R2On_V(n,conditional_indx) = p.R2On_V_reset; end
  conditional_test=any(R2Off_V(n,:) > p.R2Off_V_thresh);
  conditional_indx=(R2Off_V(n,:) > p.R2Off_V_thresh);
  if conditional_test, R2Off_V(n,conditional_indx) = p.R2Off_V_reset; R2Off_g_ad(n,conditional_indx) = R2Off_g_ad(n,conditional_indx) + p.R2Off_g_inc; end
  conditional_test=any(any(t<=R2Off_tspike+p.R2Off_t_ref,1));
  conditional_indx=(any(t<=R2Off_tspike+p.R2Off_t_ref,1));
  if conditional_test, R2Off_V(n,conditional_indx) = p.R2Off_V_reset; end
  conditional_test=any(S2OnOff_V(n,:) > p.S2OnOff_V_thresh);
  conditional_indx=(S2OnOff_V(n,:) > p.S2OnOff_V_thresh);
  if conditional_test, S2OnOff_V(n,conditional_indx) = p.S2OnOff_V_reset; S2OnOff_g_ad(n,conditional_indx) = S2OnOff_g_ad(n,conditional_indx) + p.S2OnOff_g_inc; end
  conditional_test=any(any(t<=S2OnOff_tspike+p.S2OnOff_t_ref,1));
  conditional_indx=(any(t<=S2OnOff_tspike+p.S2OnOff_t_ref,1));
  if conditional_test, S2OnOff_V(n,conditional_indx) = p.S2OnOff_V_reset; end
  conditional_test=any(any(t == On_tspike+p.R1On_On_PSC_delay,1));
  conditional_indx=(any(t == On_tspike+p.R1On_On_PSC_delay,1));
  if conditional_test, R1On_On_PSC_x(n,conditional_indx) = R1On_On_PSC_x(n,conditional_indx) + R1On_On_PSC_q(n,conditional_indx);R1On_On_PSC_q(n,conditional_indx) = R1On_On_PSC_F(n,conditional_indx).*R1On_On_PSC_P(n,conditional_indx);R1On_On_PSC_F(n,conditional_indx) = R1On_On_PSC_F(n,conditional_indx) + p.R1On_On_PSC_fF*(p.R1On_On_PSC_maxF-R1On_On_PSC_F(n,conditional_indx)); R1On_On_PSC_P(n,conditional_indx) = R1On_On_PSC_P(n,conditional_indx)*(1 - p.R1On_On_PSC_fP); end
  conditional_test=any(any(t == On_tspike+p.S1OnOff_On_PSC_delay,1));
  conditional_indx=(any(t == On_tspike+p.S1OnOff_On_PSC_delay,1));
  if conditional_test, S1OnOff_On_PSC_x(n,conditional_indx) = S1OnOff_On_PSC_x(n,conditional_indx) + S1OnOff_On_PSC_q(n,conditional_indx);S1OnOff_On_PSC_q(n,conditional_indx) = S1OnOff_On_PSC_F(n,conditional_indx).*S1OnOff_On_PSC_P(n,conditional_indx);S1OnOff_On_PSC_F(n,conditional_indx) = S1OnOff_On_PSC_F(n,conditional_indx) + p.S1OnOff_On_PSC_fF*(p.S1OnOff_On_PSC_maxF-S1OnOff_On_PSC_F(n,conditional_indx)); S1OnOff_On_PSC_P(n,conditional_indx) = S1OnOff_On_PSC_P(n,conditional_indx)*(1 - p.S1OnOff_On_PSC_fP); end
  conditional_test=any(any(t == S1OnOff_tspike+p.R1On_S1OnOff_PSC_delay,1));
  conditional_indx=(any(t == S1OnOff_tspike+p.R1On_S1OnOff_PSC_delay,1));
  if conditional_test, R1On_S1OnOff_PSC_x(n,conditional_indx) = R1On_S1OnOff_PSC_x(n,conditional_indx) + R1On_S1OnOff_PSC_q(n,conditional_indx);R1On_S1OnOff_PSC_q(n,conditional_indx) = R1On_S1OnOff_PSC_F(n,conditional_indx).*R1On_S1OnOff_PSC_P(n,conditional_indx);R1On_S1OnOff_PSC_F(n,conditional_indx) = R1On_S1OnOff_PSC_F(n,conditional_indx) + p.R1On_S1OnOff_PSC_fF*(p.R1On_S1OnOff_PSC_maxF-R1On_S1OnOff_PSC_F(n,conditional_indx)); R1On_S1OnOff_PSC_P(n,conditional_indx) = R1On_S1OnOff_PSC_P(n,conditional_indx)*(1 - p.R1On_S1OnOff_PSC_fP); end
  conditional_test=any(any(t == S1OnOff_tspike+p.R1Off_S1OnOff_PSC_delay,1));
  conditional_indx=(any(t == S1OnOff_tspike+p.R1Off_S1OnOff_PSC_delay,1));
  if conditional_test, R1Off_S1OnOff_PSC_x(n,conditional_indx) = R1Off_S1OnOff_PSC_x(n,conditional_indx) + R1Off_S1OnOff_PSC_q(n,conditional_indx);R1Off_S1OnOff_PSC_q(n,conditional_indx) = R1Off_S1OnOff_PSC_F(n,conditional_indx).*R1Off_S1OnOff_PSC_P(n,conditional_indx);R1Off_S1OnOff_PSC_F(n,conditional_indx) = R1Off_S1OnOff_PSC_F(n,conditional_indx) + p.R1Off_S1OnOff_PSC_fF*(p.R1Off_S1OnOff_PSC_maxF-R1Off_S1OnOff_PSC_F(n,conditional_indx)); R1Off_S1OnOff_PSC_P(n,conditional_indx) = R1Off_S1OnOff_PSC_P(n,conditional_indx)*(1 - p.R1Off_S1OnOff_PSC_fP); end
  conditional_test=any(any(t == Off_tspike+p.R1Off_Off_PSC_delay,1));
  conditional_indx=(any(t == Off_tspike+p.R1Off_Off_PSC_delay,1));
  if conditional_test, R1Off_Off_PSC_x(n,conditional_indx) = R1Off_Off_PSC_x(n,conditional_indx) + R1Off_Off_PSC_q(n,conditional_indx);R1Off_Off_PSC_q(n,conditional_indx) = R1Off_Off_PSC_F(n,conditional_indx).*R1Off_Off_PSC_P(n,conditional_indx);R1Off_Off_PSC_F(n,conditional_indx) = R1Off_Off_PSC_F(n,conditional_indx) + p.R1Off_Off_PSC_fF*(p.R1Off_Off_PSC_maxF-R1Off_Off_PSC_F(n,conditional_indx)); R1Off_Off_PSC_P(n,conditional_indx) = R1Off_Off_PSC_P(n,conditional_indx)*(1 - p.R1Off_Off_PSC_fP); end
  conditional_test=any(any(t == Off_tspike+p.S1OnOff_Off_PSC_delay,1));
  conditional_indx=(any(t == Off_tspike+p.S1OnOff_Off_PSC_delay,1));
  if conditional_test, S1OnOff_Off_PSC_x(n,conditional_indx) = S1OnOff_Off_PSC_x(n,conditional_indx) + S1OnOff_Off_PSC_q(n,conditional_indx);S1OnOff_Off_PSC_q(n,conditional_indx) = S1OnOff_Off_PSC_F(n,conditional_indx).*S1OnOff_Off_PSC_P(n,conditional_indx);S1OnOff_Off_PSC_F(n,conditional_indx) = S1OnOff_Off_PSC_F(n,conditional_indx) + p.S1OnOff_Off_PSC_fF*(p.S1OnOff_Off_PSC_maxF-S1OnOff_Off_PSC_F(n,conditional_indx)); S1OnOff_Off_PSC_P(n,conditional_indx) = S1OnOff_Off_PSC_P(n,conditional_indx)*(1 - p.S1OnOff_Off_PSC_fP); end
  conditional_test=any(any(t == R1On_tspike+p.R2On_R1On_PSC_delay,1));
  conditional_indx=(any(t == R1On_tspike+p.R2On_R1On_PSC_delay,1));
  if conditional_test, R2On_R1On_PSC_x(n,conditional_indx) = R2On_R1On_PSC_x(n,conditional_indx) + R2On_R1On_PSC_q(n,conditional_indx);R2On_R1On_PSC_q(n,conditional_indx) = R2On_R1On_PSC_F(n,conditional_indx).*R2On_R1On_PSC_P(n,conditional_indx);R2On_R1On_PSC_F(n,conditional_indx) = R2On_R1On_PSC_F(n,conditional_indx) + p.R2On_R1On_PSC_fF*(p.R2On_R1On_PSC_maxF-R2On_R1On_PSC_F(n,conditional_indx)); R2On_R1On_PSC_P(n,conditional_indx) = R2On_R1On_PSC_P(n,conditional_indx)*(1 - p.R2On_R1On_PSC_fP); end
  conditional_test=any(any(t == R1On_tspike+p.S2OnOff_R1On_PSC_delay,1));
  conditional_indx=(any(t == R1On_tspike+p.S2OnOff_R1On_PSC_delay,1));
  if conditional_test, S2OnOff_R1On_PSC_x(n,conditional_indx) = S2OnOff_R1On_PSC_x(n,conditional_indx) + S2OnOff_R1On_PSC_q(n,conditional_indx);S2OnOff_R1On_PSC_q(n,conditional_indx) = S2OnOff_R1On_PSC_F(n,conditional_indx).*S2OnOff_R1On_PSC_P(n,conditional_indx);S2OnOff_R1On_PSC_F(n,conditional_indx) = S2OnOff_R1On_PSC_F(n,conditional_indx) + p.S2OnOff_R1On_PSC_fF*(p.S2OnOff_R1On_PSC_maxF-S2OnOff_R1On_PSC_F(n,conditional_indx)); S2OnOff_R1On_PSC_P(n,conditional_indx) = S2OnOff_R1On_PSC_P(n,conditional_indx)*(1 - p.S2OnOff_R1On_PSC_fP); end
  conditional_test=any(any(t == S2OnOff_tspike+p.R2On_S2OnOff_PSC_delay,1));
  conditional_indx=(any(t == S2OnOff_tspike+p.R2On_S2OnOff_PSC_delay,1));
  if conditional_test, R2On_S2OnOff_PSC_x(n,conditional_indx) = R2On_S2OnOff_PSC_x(n,conditional_indx) + R2On_S2OnOff_PSC_q(n,conditional_indx);R2On_S2OnOff_PSC_q(n,conditional_indx) = R2On_S2OnOff_PSC_F(n,conditional_indx).*R2On_S2OnOff_PSC_P(n,conditional_indx);R2On_S2OnOff_PSC_F(n,conditional_indx) = R2On_S2OnOff_PSC_F(n,conditional_indx) + p.R2On_S2OnOff_PSC_fF*(p.R2On_S2OnOff_PSC_maxF-R2On_S2OnOff_PSC_F(n,conditional_indx)); R2On_S2OnOff_PSC_P(n,conditional_indx) = R2On_S2OnOff_PSC_P(n,conditional_indx)*(1 - p.R2On_S2OnOff_PSC_fP); end
  conditional_test=any(any(t == S2OnOff_tspike+p.R2Off_S2OnOff_PSC_delay,1));
  conditional_indx=(any(t == S2OnOff_tspike+p.R2Off_S2OnOff_PSC_delay,1));
  if conditional_test, R2Off_S2OnOff_PSC_x(n,conditional_indx) = R2Off_S2OnOff_PSC_x(n,conditional_indx) + R2Off_S2OnOff_PSC_q(n,conditional_indx);R2Off_S2OnOff_PSC_q(n,conditional_indx) = R2Off_S2OnOff_PSC_F(n,conditional_indx).*R2Off_S2OnOff_PSC_P(n,conditional_indx);R2Off_S2OnOff_PSC_F(n,conditional_indx) = R2Off_S2OnOff_PSC_F(n,conditional_indx) + p.R2Off_S2OnOff_PSC_fF*(p.R2Off_S2OnOff_PSC_maxF-R2Off_S2OnOff_PSC_F(n,conditional_indx)); R2Off_S2OnOff_PSC_P(n,conditional_indx) = R2Off_S2OnOff_PSC_P(n,conditional_indx)*(1 - p.R2Off_S2OnOff_PSC_fP); end
  conditional_test=any(any(t == R1Off_tspike+p.R2Off_R1Off_PSC_delay,1));
  conditional_indx=(any(t == R1Off_tspike+p.R2Off_R1Off_PSC_delay,1));
  if conditional_test, R2Off_R1Off_PSC_x(n,conditional_indx) = R2Off_R1Off_PSC_x(n,conditional_indx) + R2Off_R1Off_PSC_q(n,conditional_indx);R2Off_R1Off_PSC_q(n,conditional_indx) = R2Off_R1Off_PSC_F(n,conditional_indx).*R2Off_R1Off_PSC_P(n,conditional_indx);R2Off_R1Off_PSC_F(n,conditional_indx) = R2Off_R1Off_PSC_F(n,conditional_indx) + p.R2Off_R1Off_PSC_fF*(p.R2Off_R1Off_PSC_maxF-R2Off_R1Off_PSC_F(n,conditional_indx)); R2Off_R1Off_PSC_P(n,conditional_indx) = R2Off_R1Off_PSC_P(n,conditional_indx)*(1 - p.R2Off_R1Off_PSC_fP); end
  conditional_test=any(any(t == R1Off_tspike+p.S2OnOff_R1Off_PSC_delay,1));
  conditional_indx=(any(t == R1Off_tspike+p.S2OnOff_R1Off_PSC_delay,1));
  if conditional_test, S2OnOff_R1Off_PSC_x(n,conditional_indx) = S2OnOff_R1Off_PSC_x(n,conditional_indx) + S2OnOff_R1Off_PSC_q(n,conditional_indx);S2OnOff_R1Off_PSC_q(n,conditional_indx) = S2OnOff_R1Off_PSC_F(n,conditional_indx).*S2OnOff_R1Off_PSC_P(n,conditional_indx);S2OnOff_R1Off_PSC_F(n,conditional_indx) = S2OnOff_R1Off_PSC_F(n,conditional_indx) + p.S2OnOff_R1Off_PSC_fF*(p.S2OnOff_R1Off_PSC_maxF-S2OnOff_R1Off_PSC_F(n,conditional_indx)); S2OnOff_R1Off_PSC_P(n,conditional_indx) = S2OnOff_R1Off_PSC_P(n,conditional_indx)*(1 - p.S2OnOff_R1Off_PSC_fP); end

  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  On_On_IC_iIC(n,:)=p.On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(n,:)-p.On_On_IC_E_exc);
  Off_Off_IC_iIC(n,:)=p.Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(n,:)-p.Off_Off_IC_E_exc);
  R1On_On_PSC_syn(n,:)=p.R1On_On_PSC_gSYN.*(R1On_On_PSC_s(n,:)*R1On_On_PSC_netcon).*(R1On_V(n,:)-p.R1On_On_PSC_ESYN);
  S1OnOff_On_PSC_syn(n,:)=p.S1OnOff_On_PSC_gSYN.*(S1OnOff_On_PSC_s(n,:)*S1OnOff_On_PSC_netcon).*(S1OnOff_V(n,:)-p.S1OnOff_On_PSC_ESYN);
  R1On_S1OnOff_PSC_syn(n,:)=p.R1On_S1OnOff_PSC_gSYN.*(R1On_S1OnOff_PSC_s(n,:)*R1On_S1OnOff_PSC_netcon).*(R1On_V(n,:)-p.R1On_S1OnOff_PSC_ESYN);
  R1Off_S1OnOff_PSC_syn(n,:)=p.R1Off_S1OnOff_PSC_gSYN.*(R1Off_S1OnOff_PSC_s(n,:)*R1Off_S1OnOff_PSC_netcon).*(R1Off_V(n,:)-p.R1Off_S1OnOff_PSC_ESYN);
  R1Off_Off_PSC_syn(n,:)=p.R1Off_Off_PSC_gSYN.*(R1Off_Off_PSC_s(n,:)*R1Off_Off_PSC_netcon).*(R1Off_V(n,:)-p.R1Off_Off_PSC_ESYN);
  S1OnOff_Off_PSC_syn(n,:)=p.S1OnOff_Off_PSC_gSYN.*(S1OnOff_Off_PSC_s(n,:)*S1OnOff_Off_PSC_netcon).*(S1OnOff_V(n,:)-p.S1OnOff_Off_PSC_ESYN);
  R2On_R1On_PSC_syn(n,:)=p.R2On_R1On_PSC_gSYN.*(R2On_R1On_PSC_s(n,:)*R2On_R1On_PSC_netcon).*(R2On_V(n,:)-p.R2On_R1On_PSC_ESYN);
  S2OnOff_R1On_PSC_syn(n,:)=p.S2OnOff_R1On_PSC_gSYN.*(S2OnOff_R1On_PSC_s(n,:)*S2OnOff_R1On_PSC_netcon).*(S2OnOff_V(n,:)-p.S2OnOff_R1On_PSC_ESYN);
  R2On_S2OnOff_PSC_syn(n,:)=p.R2On_S2OnOff_PSC_gSYN.*(R2On_S2OnOff_PSC_s(n,:)*R2On_S2OnOff_PSC_netcon).*(R2On_V(n,:)-p.R2On_S2OnOff_PSC_ESYN);
  R2Off_S2OnOff_PSC_syn(n,:)=p.R2Off_S2OnOff_PSC_gSYN.*(R2Off_S2OnOff_PSC_s(n,:)*R2Off_S2OnOff_PSC_netcon).*(R2Off_V(n,:)-p.R2Off_S2OnOff_PSC_ESYN);
  R2Off_R1Off_PSC_syn(n,:)=p.R2Off_R1Off_PSC_gSYN.*(R2Off_R1Off_PSC_s(n,:)*R2Off_R1Off_PSC_netcon).*(R2Off_V(n,:)-p.R2Off_R1Off_PSC_ESYN);
  S2OnOff_R1Off_PSC_syn(n,:)=p.S2OnOff_R1Off_PSC_gSYN.*(S2OnOff_R1Off_PSC_s(n,:)*S2OnOff_R1Off_PSC_netcon).*(S2OnOff_V(n,:)-p.S2OnOff_R1Off_PSC_ESYN);
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end
