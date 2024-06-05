function [T,On_V,On_g_ad,Off_V,Off_g_ad,ROn_V,ROn_g_ad,ROff_V,ROff_g_ad,SOnOff_V,SOnOff_g_ad,TD_V,TD_g_ad,X_V,X_g_ad,C_V,C_g_ad,ROn_On_PSC3_s,ROn_On_PSC3_x,ROn_On_PSC3_F,ROn_On_PSC3_P,ROn_On_PSC3_q,SOnOff_On_PSC3_s,SOnOff_On_PSC3_x,SOnOff_On_PSC3_F,SOnOff_On_PSC3_P,SOnOff_On_PSC3_q,ROn_SOnOff_PSC3_s,ROn_SOnOff_PSC3_x,ROn_SOnOff_PSC3_F,ROn_SOnOff_PSC3_P,ROn_SOnOff_PSC3_q,ROff_SOnOff_PSC3_s,ROff_SOnOff_PSC3_x,ROff_SOnOff_PSC3_F,ROff_SOnOff_PSC3_P,ROff_SOnOff_PSC3_q,ROff_Off_PSC3_s,ROff_Off_PSC3_x,ROff_Off_PSC3_F,ROff_Off_PSC3_P,ROff_Off_PSC3_q,SOnOff_Off_PSC3_s,SOnOff_Off_PSC3_x,SOnOff_Off_PSC3_F,SOnOff_Off_PSC3_P,SOnOff_Off_PSC3_q,ROn_ROn_iNoise_V3_sn,ROn_ROn_iNoise_V3_xn,X_ROn_PSC3_s,X_ROn_PSC3_x,X_ROn_PSC3_F,X_ROn_PSC3_P,X_ROn_PSC3_q,ROn_X_PSC3_s,ROn_X_PSC3_x,ROn_X_PSC3_F,ROn_X_PSC3_P,ROn_X_PSC3_q,ROn_TD_PSC3_s,ROn_TD_PSC3_x,ROn_TD_PSC3_F,ROn_TD_PSC3_P,ROn_TD_PSC3_q,ROff_TD_PSC3_s,ROff_TD_PSC3_x,ROff_TD_PSC3_F,ROff_TD_PSC3_P,ROff_TD_PSC3_q,X_TD_PSC3_s,X_TD_PSC3_x,X_TD_PSC3_F,X_TD_PSC3_P,X_TD_PSC3_q,C_ROn_PSC3_s,C_ROn_PSC3_x,C_ROn_PSC3_F,C_ROn_PSC3_P,C_ROn_PSC3_q,On_V_spikes,Off_V_spikes,ROn_V_spikes,ROff_V_spikes,SOnOff_V_spikes,TD_V_spikes,X_V_spikes,C_V_spikes,On_On_IC_iIC,Off_Off_IC_iIC,ROn_On_PSC3_syn,SOnOff_On_PSC3_syn,ROn_SOnOff_PSC3_syn,ROff_SOnOff_PSC3_syn,ROff_Off_PSC3_syn,SOnOff_Off_PSC3_syn,X_ROn_PSC3_syn,ROn_X_PSC3_syn,ROn_TD_PSC3_syn,ROff_TD_PSC3_syn,X_TD_PSC3_syn,C_ROn_PSC3_syn,On_R,On_tau,On_Imask,Off_R,Off_tau,Off_Imask,ROn_R,ROn_tau,ROn_Imask,ROff_R,ROff_tau,ROff_Imask,SOnOff_R,SOnOff_tau,SOnOff_Imask,TD_R,TD_tau,TD_Imask,TD_Icur,X_R,X_tau,X_Imask,C_R,C_tau,C_Imask,On_On_IC_netcon,On_On_IC_input,Off_Off_IC_netcon,Off_Off_IC_input,ROn_On_PSC3_netcon,ROn_On_PSC3_scale,SOnOff_On_PSC3_netcon,SOnOff_On_PSC3_scale,ROn_SOnOff_PSC3_netcon,ROn_SOnOff_PSC3_scale,ROff_SOnOff_PSC3_netcon,ROff_SOnOff_PSC3_scale,ROff_Off_PSC3_netcon,ROff_Off_PSC3_scale,SOnOff_Off_PSC3_netcon,SOnOff_Off_PSC3_scale,ROn_ROn_iNoise_V3_netcon,ROn_ROn_iNoise_V3_token,ROn_ROn_iNoise_V3_scale,X_ROn_PSC3_netcon,X_ROn_PSC3_scale,ROn_X_PSC3_netcon,ROn_X_PSC3_scale,ROn_TD_PSC3_netcon,ROn_TD_PSC3_scale,ROff_TD_PSC3_netcon,ROff_TD_PSC3_scale,X_TD_PSC3_netcon,X_TD_PSC3_scale,C_ROn_PSC3_netcon,C_ROn_PSC3_scale]=solve_ode

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
ROn_R = 1/p.ROn_g_L;
ROn_tau = p.ROn_C*ROn_R;
ROn_Imask =  ones(1,p.ROn_Npop);
ROff_R = 1/p.ROff_g_L;
ROff_tau = p.ROff_C*ROff_R;
ROff_Imask =  ones(1,p.ROff_Npop);
SOnOff_R = 1/p.SOnOff_g_L;
SOnOff_tau = p.SOnOff_C*SOnOff_R;
SOnOff_Imask =  ones(1,p.SOnOff_Npop);
TD_R = 1/p.TD_g_L;
TD_tau = p.TD_C*TD_R;
TD_Imask =  ones(1,p.TD_Npop);
TD_Icur =  p.TD_Itonic*buildTonicCurrent(T,p.TD_Npop,dt,p.TD_numLocs);
X_R = 1/p.X_g_L;
X_tau = p.X_C*X_R;
X_Imask =  ones(1,p.X_Npop);
C_R = 1/p.C_g_L;
C_tau = p.C_C*C_R;
C_Imask =  ones(1,p.C_Npop);
On_On_IC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
On_On_IC_input =  genPoissonInputs(p.On_On_IC_trial,p.On_On_IC_locNum,p.On_On_IC_label,p.On_On_IC_t_ref,p.On_On_IC_t_ref_rel,p.On_On_IC_rec);
Off_Off_IC_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
Off_Off_IC_input =  genPoissonInputs(p.Off_Off_IC_trial,p.Off_Off_IC_locNum,p.Off_Off_IC_label,p.Off_Off_IC_t_ref,p.Off_Off_IC_t_ref_rel,p.Off_Off_IC_rec);
ROn_On_PSC3_netcon = eye(p.On_Npop,p.ROn_Npop);
ROn_On_PSC3_scale = (p.ROn_On_PSC3_tauD/p.ROn_On_PSC3_tauR)^(p.ROn_On_PSC3_tauR/(p.ROn_On_PSC3_tauD-p.ROn_On_PSC3_tauR));
SOnOff_On_PSC3_netcon = [+1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00];
SOnOff_On_PSC3_scale = (p.SOnOff_On_PSC3_tauD/p.SOnOff_On_PSC3_tauR)^(p.SOnOff_On_PSC3_tauR/(p.SOnOff_On_PSC3_tauD-p.SOnOff_On_PSC3_tauR));
ROn_SOnOff_PSC3_netcon = [+1.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00];
ROn_SOnOff_PSC3_scale = (p.ROn_SOnOff_PSC3_tauD/p.ROn_SOnOff_PSC3_tauR)^(p.ROn_SOnOff_PSC3_tauR/(p.ROn_SOnOff_PSC3_tauD-p.ROn_SOnOff_PSC3_tauR));
ROff_SOnOff_PSC3_netcon = [+1.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00; +0.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00];
ROff_SOnOff_PSC3_scale = (p.ROff_SOnOff_PSC3_tauD/p.ROff_SOnOff_PSC3_tauR)^(p.ROff_SOnOff_PSC3_tauR/(p.ROff_SOnOff_PSC3_tauD-p.ROff_SOnOff_PSC3_tauR));
ROff_Off_PSC3_netcon = eye(p.Off_Npop,p.ROff_Npop);
ROff_Off_PSC3_scale = (p.ROff_Off_PSC3_tauD/p.ROff_Off_PSC3_tauR)^(p.ROff_Off_PSC3_tauR/(p.ROff_Off_PSC3_tauD-p.ROff_Off_PSC3_tauR));
SOnOff_Off_PSC3_netcon = eye(p.Off_Npop,p.SOnOff_Npop);
SOnOff_Off_PSC3_scale = (p.SOnOff_Off_PSC3_tauD/p.SOnOff_Off_PSC3_tauR)^(p.SOnOff_Off_PSC3_tauR/(p.SOnOff_Off_PSC3_tauD-p.SOnOff_Off_PSC3_tauR));
ROn_ROn_iNoise_V3_netcon =  eye(p.ROn_Npop,p.ROn_Npop);
ROn_ROn_iNoise_V3_token = genPoissonTimes(p.ROn_Npop,p.ROn_ROn_iNoise_V3_dt,p.ROn_ROn_iNoise_V3_FR,p.ROn_ROn_iNoise_V3_sigma,p.ROn_ROn_iNoise_V3_simlen);
ROn_ROn_iNoise_V3_scale =  (p.ROn_ROn_iNoise_V3_tauD_N/p.ROn_ROn_iNoise_V3_tauR_N)^(p.ROn_ROn_iNoise_V3_tauR_N/(p.ROn_ROn_iNoise_V3_tauD_N-p.ROn_ROn_iNoise_V3_tauR_N));
X_ROn_PSC3_netcon = eye(p.ROn_Npop,p.X_Npop);
X_ROn_PSC3_scale = (p.X_ROn_PSC3_tauD/p.X_ROn_PSC3_tauR)^(p.X_ROn_PSC3_tauR/(p.X_ROn_PSC3_tauD-p.X_ROn_PSC3_tauR));
ROn_X_PSC3_netcon = [+0.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00; +1.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00; +1.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00   +1.000000000000000e+00; +1.000000000000000e+00   +1.000000000000000e+00   +1.000000000000000e+00   +0.000000000000000e+00];
ROn_X_PSC3_scale = (p.ROn_X_PSC3_tauD/p.ROn_X_PSC3_tauR)^(p.ROn_X_PSC3_tauR/(p.ROn_X_PSC3_tauD-p.ROn_X_PSC3_tauR));
ROn_TD_PSC3_netcon = eye(p.TD_Npop,p.ROn_Npop);
ROn_TD_PSC3_scale = (p.ROn_TD_PSC3_tauD/p.ROn_TD_PSC3_tauR)^(p.ROn_TD_PSC3_tauR/(p.ROn_TD_PSC3_tauD-p.ROn_TD_PSC3_tauR));
ROff_TD_PSC3_netcon = eye(p.TD_Npop,p.ROff_Npop);
ROff_TD_PSC3_scale = (p.ROff_TD_PSC3_tauD/p.ROff_TD_PSC3_tauR)^(p.ROff_TD_PSC3_tauR/(p.ROff_TD_PSC3_tauD-p.ROff_TD_PSC3_tauR));
X_TD_PSC3_netcon = eye(p.TD_Npop,p.X_Npop);
X_TD_PSC3_scale = (p.X_TD_PSC3_tauD/p.X_TD_PSC3_tauR)^(p.X_TD_PSC3_tauR/(p.X_TD_PSC3_tauD-p.X_TD_PSC3_tauR));
C_ROn_PSC3_netcon = [+1.000000000000000e+00; +1.000000000000000e+00; +1.000000000000000e+00; +1.000000000000000e+00];
C_ROn_PSC3_scale = (p.C_ROn_PSC3_tauD/p.C_ROn_PSC3_tauR)^(p.C_ROn_PSC3_tauR/(p.C_ROn_PSC3_tauD-p.C_ROn_PSC3_tauR));

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
ROn_V = zeros(nsamp,p.ROn_Npop);
ROn_V(1,:) =  p.ROn_E_L*ones(1,p.ROn_Npop);
ROn_g_ad = zeros(nsamp,p.ROn_Npop);
ROn_g_ad(1,:) =  zeros(1,p.ROn_Npop);
ROff_V = zeros(nsamp,p.ROff_Npop);
ROff_V(1,:) =  p.ROff_E_L*ones(1,p.ROff_Npop);
ROff_g_ad = zeros(nsamp,p.ROff_Npop);
ROff_g_ad(1,:) =  zeros(1,p.ROff_Npop);
SOnOff_V = zeros(nsamp,p.SOnOff_Npop);
SOnOff_V(1,:) =  p.SOnOff_E_L*ones(1,p.SOnOff_Npop);
SOnOff_g_ad = zeros(nsamp,p.SOnOff_Npop);
SOnOff_g_ad(1,:) =  zeros(1,p.SOnOff_Npop);
TD_V = zeros(nsamp,p.TD_Npop);
TD_V(1,:) =  p.TD_E_L*ones(1,p.TD_Npop);
TD_g_ad = zeros(nsamp,p.TD_Npop);
TD_g_ad(1,:) =  zeros(1,p.TD_Npop);
X_V = zeros(nsamp,p.X_Npop);
X_V(1,:) =  p.X_E_L*ones(1,p.X_Npop);
X_g_ad = zeros(nsamp,p.X_Npop);
X_g_ad(1,:) =  zeros(1,p.X_Npop);
C_V = zeros(nsamp,p.C_Npop);
C_V(1,:) =  p.C_E_L*ones(1,p.C_Npop);
C_g_ad = zeros(nsamp,p.C_Npop);
C_g_ad(1,:) =  zeros(1,p.C_Npop);
ROn_On_PSC3_s = zeros(nsamp,p.On_Npop);
ROn_On_PSC3_s(1,:) =  zeros(1,p.On_Npop);
ROn_On_PSC3_x = zeros(nsamp,p.On_Npop);
ROn_On_PSC3_x(1,:) =  zeros(1,p.On_Npop);
ROn_On_PSC3_F = zeros(nsamp,p.On_Npop);
ROn_On_PSC3_F(1,:) =  ones(1,p.On_Npop);
ROn_On_PSC3_P = zeros(nsamp,p.On_Npop);
ROn_On_PSC3_P(1,:) =  ones(1,p.On_Npop);
ROn_On_PSC3_q = zeros(nsamp,p.On_Npop);
ROn_On_PSC3_q(1,:) =  ones(1,p.On_Npop);
SOnOff_On_PSC3_s = zeros(nsamp,p.On_Npop);
SOnOff_On_PSC3_s(1,:) =  zeros(1,p.On_Npop);
SOnOff_On_PSC3_x = zeros(nsamp,p.On_Npop);
SOnOff_On_PSC3_x(1,:) =  zeros(1,p.On_Npop);
SOnOff_On_PSC3_F = zeros(nsamp,p.On_Npop);
SOnOff_On_PSC3_F(1,:) =  ones(1,p.On_Npop);
SOnOff_On_PSC3_P = zeros(nsamp,p.On_Npop);
SOnOff_On_PSC3_P(1,:) =  ones(1,p.On_Npop);
SOnOff_On_PSC3_q = zeros(nsamp,p.On_Npop);
SOnOff_On_PSC3_q(1,:) =  ones(1,p.On_Npop);
ROn_SOnOff_PSC3_s = zeros(nsamp,p.SOnOff_Npop);
ROn_SOnOff_PSC3_s(1,:) =  zeros(1,p.SOnOff_Npop);
ROn_SOnOff_PSC3_x = zeros(nsamp,p.SOnOff_Npop);
ROn_SOnOff_PSC3_x(1,:) =  zeros(1,p.SOnOff_Npop);
ROn_SOnOff_PSC3_F = zeros(nsamp,p.SOnOff_Npop);
ROn_SOnOff_PSC3_F(1,:) =  ones(1,p.SOnOff_Npop);
ROn_SOnOff_PSC3_P = zeros(nsamp,p.SOnOff_Npop);
ROn_SOnOff_PSC3_P(1,:) =  ones(1,p.SOnOff_Npop);
ROn_SOnOff_PSC3_q = zeros(nsamp,p.SOnOff_Npop);
ROn_SOnOff_PSC3_q(1,:) =  ones(1,p.SOnOff_Npop);
ROff_SOnOff_PSC3_s = zeros(nsamp,p.SOnOff_Npop);
ROff_SOnOff_PSC3_s(1,:) =  zeros(1,p.SOnOff_Npop);
ROff_SOnOff_PSC3_x = zeros(nsamp,p.SOnOff_Npop);
ROff_SOnOff_PSC3_x(1,:) =  zeros(1,p.SOnOff_Npop);
ROff_SOnOff_PSC3_F = zeros(nsamp,p.SOnOff_Npop);
ROff_SOnOff_PSC3_F(1,:) =  ones(1,p.SOnOff_Npop);
ROff_SOnOff_PSC3_P = zeros(nsamp,p.SOnOff_Npop);
ROff_SOnOff_PSC3_P(1,:) =  ones(1,p.SOnOff_Npop);
ROff_SOnOff_PSC3_q = zeros(nsamp,p.SOnOff_Npop);
ROff_SOnOff_PSC3_q(1,:) =  ones(1,p.SOnOff_Npop);
ROff_Off_PSC3_s = zeros(nsamp,p.Off_Npop);
ROff_Off_PSC3_s(1,:) =  zeros(1,p.Off_Npop);
ROff_Off_PSC3_x = zeros(nsamp,p.Off_Npop);
ROff_Off_PSC3_x(1,:) =  zeros(1,p.Off_Npop);
ROff_Off_PSC3_F = zeros(nsamp,p.Off_Npop);
ROff_Off_PSC3_F(1,:) =  ones(1,p.Off_Npop);
ROff_Off_PSC3_P = zeros(nsamp,p.Off_Npop);
ROff_Off_PSC3_P(1,:) =  ones(1,p.Off_Npop);
ROff_Off_PSC3_q = zeros(nsamp,p.Off_Npop);
ROff_Off_PSC3_q(1,:) =  ones(1,p.Off_Npop);
SOnOff_Off_PSC3_s = zeros(nsamp,p.Off_Npop);
SOnOff_Off_PSC3_s(1,:) =  zeros(1,p.Off_Npop);
SOnOff_Off_PSC3_x = zeros(nsamp,p.Off_Npop);
SOnOff_Off_PSC3_x(1,:) =  zeros(1,p.Off_Npop);
SOnOff_Off_PSC3_F = zeros(nsamp,p.Off_Npop);
SOnOff_Off_PSC3_F(1,:) =  ones(1,p.Off_Npop);
SOnOff_Off_PSC3_P = zeros(nsamp,p.Off_Npop);
SOnOff_Off_PSC3_P(1,:) =  ones(1,p.Off_Npop);
SOnOff_Off_PSC3_q = zeros(nsamp,p.Off_Npop);
SOnOff_Off_PSC3_q(1,:) =  ones(1,p.Off_Npop);
ROn_ROn_iNoise_V3_sn = zeros(nsamp,p.ROn_Npop);
ROn_ROn_iNoise_V3_sn(1,:) =  0 * ones(1,p.ROn_Npop);
ROn_ROn_iNoise_V3_xn = zeros(nsamp,p.ROn_Npop);
ROn_ROn_iNoise_V3_xn(1,:) =  0 * ones(1,p.ROn_Npop);
X_ROn_PSC3_s = zeros(nsamp,p.ROn_Npop);
X_ROn_PSC3_s(1,:) =  zeros(1,p.ROn_Npop);
X_ROn_PSC3_x = zeros(nsamp,p.ROn_Npop);
X_ROn_PSC3_x(1,:) =  zeros(1,p.ROn_Npop);
X_ROn_PSC3_F = zeros(nsamp,p.ROn_Npop);
X_ROn_PSC3_F(1,:) =  ones(1,p.ROn_Npop);
X_ROn_PSC3_P = zeros(nsamp,p.ROn_Npop);
X_ROn_PSC3_P(1,:) =  ones(1,p.ROn_Npop);
X_ROn_PSC3_q = zeros(nsamp,p.ROn_Npop);
X_ROn_PSC3_q(1,:) =  ones(1,p.ROn_Npop);
ROn_X_PSC3_s = zeros(nsamp,p.X_Npop);
ROn_X_PSC3_s(1,:) =  zeros(1,p.X_Npop);
ROn_X_PSC3_x = zeros(nsamp,p.X_Npop);
ROn_X_PSC3_x(1,:) =  zeros(1,p.X_Npop);
ROn_X_PSC3_F = zeros(nsamp,p.X_Npop);
ROn_X_PSC3_F(1,:) =  ones(1,p.X_Npop);
ROn_X_PSC3_P = zeros(nsamp,p.X_Npop);
ROn_X_PSC3_P(1,:) =  ones(1,p.X_Npop);
ROn_X_PSC3_q = zeros(nsamp,p.X_Npop);
ROn_X_PSC3_q(1,:) =  ones(1,p.X_Npop);
ROn_TD_PSC3_s = zeros(nsamp,p.TD_Npop);
ROn_TD_PSC3_s(1,:) =  zeros(1,p.TD_Npop);
ROn_TD_PSC3_x = zeros(nsamp,p.TD_Npop);
ROn_TD_PSC3_x(1,:) =  zeros(1,p.TD_Npop);
ROn_TD_PSC3_F = zeros(nsamp,p.TD_Npop);
ROn_TD_PSC3_F(1,:) =  ones(1,p.TD_Npop);
ROn_TD_PSC3_P = zeros(nsamp,p.TD_Npop);
ROn_TD_PSC3_P(1,:) =  ones(1,p.TD_Npop);
ROn_TD_PSC3_q = zeros(nsamp,p.TD_Npop);
ROn_TD_PSC3_q(1,:) =  ones(1,p.TD_Npop);
ROff_TD_PSC3_s = zeros(nsamp,p.TD_Npop);
ROff_TD_PSC3_s(1,:) =  zeros(1,p.TD_Npop);
ROff_TD_PSC3_x = zeros(nsamp,p.TD_Npop);
ROff_TD_PSC3_x(1,:) =  zeros(1,p.TD_Npop);
ROff_TD_PSC3_F = zeros(nsamp,p.TD_Npop);
ROff_TD_PSC3_F(1,:) =  ones(1,p.TD_Npop);
ROff_TD_PSC3_P = zeros(nsamp,p.TD_Npop);
ROff_TD_PSC3_P(1,:) =  ones(1,p.TD_Npop);
ROff_TD_PSC3_q = zeros(nsamp,p.TD_Npop);
ROff_TD_PSC3_q(1,:) =  ones(1,p.TD_Npop);
X_TD_PSC3_s = zeros(nsamp,p.TD_Npop);
X_TD_PSC3_s(1,:) =  zeros(1,p.TD_Npop);
X_TD_PSC3_x = zeros(nsamp,p.TD_Npop);
X_TD_PSC3_x(1,:) =  zeros(1,p.TD_Npop);
X_TD_PSC3_F = zeros(nsamp,p.TD_Npop);
X_TD_PSC3_F(1,:) =  ones(1,p.TD_Npop);
X_TD_PSC3_P = zeros(nsamp,p.TD_Npop);
X_TD_PSC3_P(1,:) =  ones(1,p.TD_Npop);
X_TD_PSC3_q = zeros(nsamp,p.TD_Npop);
X_TD_PSC3_q(1,:) =  ones(1,p.TD_Npop);
C_ROn_PSC3_s = zeros(nsamp,p.ROn_Npop);
C_ROn_PSC3_s(1,:) =  zeros(1,p.ROn_Npop);
C_ROn_PSC3_x = zeros(nsamp,p.ROn_Npop);
C_ROn_PSC3_x(1,:) =  zeros(1,p.ROn_Npop);
C_ROn_PSC3_F = zeros(nsamp,p.ROn_Npop);
C_ROn_PSC3_F(1,:) =  ones(1,p.ROn_Npop);
C_ROn_PSC3_P = zeros(nsamp,p.ROn_Npop);
C_ROn_PSC3_P(1,:) =  ones(1,p.ROn_Npop);
C_ROn_PSC3_q = zeros(nsamp,p.ROn_Npop);
C_ROn_PSC3_q(1,:) =  ones(1,p.ROn_Npop);

% MONITORS:
On_tspike = -1e32*ones(5,p.On_Npop);
On_buffer_index = ones(1,p.On_Npop);
On_V_spikes = zeros(nsamp,p.On_Npop);
Off_tspike = -1e32*ones(5,p.Off_Npop);
Off_buffer_index = ones(1,p.Off_Npop);
Off_V_spikes = zeros(nsamp,p.Off_Npop);
ROn_tspike = -1e32*ones(5,p.ROn_Npop);
ROn_buffer_index = ones(1,p.ROn_Npop);
ROn_V_spikes = zeros(nsamp,p.ROn_Npop);
ROff_tspike = -1e32*ones(5,p.ROff_Npop);
ROff_buffer_index = ones(1,p.ROff_Npop);
ROff_V_spikes = zeros(nsamp,p.ROff_Npop);
SOnOff_tspike = -1e32*ones(5,p.SOnOff_Npop);
SOnOff_buffer_index = ones(1,p.SOnOff_Npop);
SOnOff_V_spikes = zeros(nsamp,p.SOnOff_Npop);
TD_tspike = -1e32*ones(5,p.TD_Npop);
TD_buffer_index = ones(1,p.TD_Npop);
TD_V_spikes = zeros(nsamp,p.TD_Npop);
X_tspike = -1e32*ones(5,p.X_Npop);
X_buffer_index = ones(1,p.X_Npop);
X_V_spikes = zeros(nsamp,p.X_Npop);
C_tspike = -1e32*ones(5,p.C_Npop);
C_buffer_index = ones(1,p.C_Npop);
C_V_spikes = zeros(nsamp,p.C_Npop);
On_On_IC_iIC = zeros(nsamp,p.On_Npop);
On_On_IC_iIC(1,:)=p.On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(1,:)-p.On_On_IC_E_exc);
Off_Off_IC_iIC = zeros(nsamp,p.Off_Npop);
Off_Off_IC_iIC(1,:)=p.Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(1,:)-p.Off_Off_IC_E_exc);
ROn_On_PSC3_syn = zeros(nsamp,p.ROn_Npop);
ROn_On_PSC3_syn(1,:)=((ROn_On_PSC3_s(1,:)*ROn_On_PSC3_netcon).*(ROn_V(1,:)-p.ROn_On_PSC3_ESYN))*p.ROn_On_PSC3_gSYN;
SOnOff_On_PSC3_syn = zeros(nsamp,p.SOnOff_Npop);
SOnOff_On_PSC3_syn(1,:)=((SOnOff_On_PSC3_s(1,:)*SOnOff_On_PSC3_netcon).*(SOnOff_V(1,:)-p.SOnOff_On_PSC3_ESYN))*p.SOnOff_On_PSC3_gSYN;
ROn_SOnOff_PSC3_syn = zeros(nsamp,p.ROn_Npop);
ROn_SOnOff_PSC3_syn(1,:)=((ROn_SOnOff_PSC3_s(1,:)*ROn_SOnOff_PSC3_netcon).*(ROn_V(1,:)-p.ROn_SOnOff_PSC3_ESYN))*p.ROn_SOnOff_PSC3_gSYN;
ROff_SOnOff_PSC3_syn = zeros(nsamp,p.ROff_Npop);
ROff_SOnOff_PSC3_syn(1,:)=((ROff_SOnOff_PSC3_s(1,:)*ROff_SOnOff_PSC3_netcon).*(ROff_V(1,:)-p.ROff_SOnOff_PSC3_ESYN))*p.ROff_SOnOff_PSC3_gSYN;
ROff_Off_PSC3_syn = zeros(nsamp,p.ROff_Npop);
ROff_Off_PSC3_syn(1,:)=((ROff_Off_PSC3_s(1,:)*ROff_Off_PSC3_netcon).*(ROff_V(1,:)-p.ROff_Off_PSC3_ESYN))*p.ROff_Off_PSC3_gSYN;
SOnOff_Off_PSC3_syn = zeros(nsamp,p.SOnOff_Npop);
SOnOff_Off_PSC3_syn(1,:)=((SOnOff_Off_PSC3_s(1,:)*SOnOff_Off_PSC3_netcon).*(SOnOff_V(1,:)-p.SOnOff_Off_PSC3_ESYN))*p.SOnOff_Off_PSC3_gSYN;
X_ROn_PSC3_syn = zeros(nsamp,p.X_Npop);
X_ROn_PSC3_syn(1,:)=((X_ROn_PSC3_s(1,:)*X_ROn_PSC3_netcon).*(X_V(1,:)-p.X_ROn_PSC3_ESYN))*p.X_ROn_PSC3_gSYN;
ROn_X_PSC3_syn = zeros(nsamp,p.ROn_Npop);
ROn_X_PSC3_syn(1,:)=((ROn_X_PSC3_s(1,:)*ROn_X_PSC3_netcon).*(ROn_V(1,:)-p.ROn_X_PSC3_ESYN))*p.ROn_X_PSC3_gSYN;
ROn_TD_PSC3_syn = zeros(nsamp,p.ROn_Npop);
ROn_TD_PSC3_syn(1,:)=((ROn_TD_PSC3_s(1,:)*ROn_TD_PSC3_netcon).*(ROn_V(1,:)-p.ROn_TD_PSC3_ESYN))*p.ROn_TD_PSC3_gSYN;
ROff_TD_PSC3_syn = zeros(nsamp,p.ROff_Npop);
ROff_TD_PSC3_syn(1,:)=((ROff_TD_PSC3_s(1,:)*ROff_TD_PSC3_netcon).*(ROff_V(1,:)-p.ROff_TD_PSC3_ESYN))*p.ROff_TD_PSC3_gSYN;
X_TD_PSC3_syn = zeros(nsamp,p.X_Npop);
X_TD_PSC3_syn(1,:)=((X_TD_PSC3_s(1,:)*X_TD_PSC3_netcon).*(X_V(1,:)-p.X_TD_PSC3_ESYN))*p.X_TD_PSC3_gSYN;
C_ROn_PSC3_syn = zeros(nsamp,p.C_Npop);
C_ROn_PSC3_syn(1,:)=((C_ROn_PSC3_s(1,:)*C_ROn_PSC3_netcon).*(C_V(1,:)-p.C_ROn_PSC3_ESYN))*p.C_ROn_PSC3_gSYN;

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
  ROn_V_k1 = ( (p.ROn_E_L-ROn_V(n-1,:)) - ROn_R*ROn_g_ad(n-1,:).*(ROn_V(n-1,:)-p.ROn_E_k) - ROn_R*((((((ROn_On_PSC3_s(n-1,:)*ROn_On_PSC3_netcon).*(ROn_V(n-1,:)-p.ROn_On_PSC3_ESYN))*p.ROn_On_PSC3_gSYN)))+((((((ROn_SOnOff_PSC3_s(n-1,:)*ROn_SOnOff_PSC3_netcon).*(ROn_V(n-1,:)-p.ROn_SOnOff_PSC3_ESYN))*p.ROn_SOnOff_PSC3_gSYN)))+((((p.ROn_ROn_iNoise_V3_nSYN.*(ROn_ROn_iNoise_V3_sn(n-1,:)*ROn_ROn_iNoise_V3_netcon).*(ROn_V(n-1,:)-p.ROn_ROn_iNoise_V3_E_exc))))+((((((ROn_X_PSC3_s(n-1,:)*ROn_X_PSC3_netcon).*(ROn_V(n-1,:)-p.ROn_X_PSC3_ESYN))*p.ROn_X_PSC3_gSYN)))+((((((ROn_TD_PSC3_s(n-1,:)*ROn_TD_PSC3_netcon).*(ROn_V(n-1,:)-p.ROn_TD_PSC3_ESYN))*p.ROn_TD_PSC3_gSYN)))))))) + ROn_R*p.ROn_Itonic.*ROn_Imask + ROn_R*p.ROn_noise.*randn(1,p.ROn_Npop) ) / ROn_tau;
  ROn_g_ad_k1 = -ROn_g_ad(n-1,:) / p.ROn_tau_ad;
  ROff_V_k1 = ( (p.ROff_E_L-ROff_V(n-1,:)) - ROff_R*ROff_g_ad(n-1,:).*(ROff_V(n-1,:)-p.ROff_E_k) - ROff_R*((((((ROff_SOnOff_PSC3_s(n-1,:)*ROff_SOnOff_PSC3_netcon).*(ROff_V(n-1,:)-p.ROff_SOnOff_PSC3_ESYN))*p.ROff_SOnOff_PSC3_gSYN)))+((((((ROff_Off_PSC3_s(n-1,:)*ROff_Off_PSC3_netcon).*(ROff_V(n-1,:)-p.ROff_Off_PSC3_ESYN))*p.ROff_Off_PSC3_gSYN)))+((((((ROff_TD_PSC3_s(n-1,:)*ROff_TD_PSC3_netcon).*(ROff_V(n-1,:)-p.ROff_TD_PSC3_ESYN))*p.ROff_TD_PSC3_gSYN)))))) + ROff_R*p.ROff_Itonic.*ROff_Imask + ROff_R*p.ROff_noise.*randn(1,p.ROff_Npop) ) / ROff_tau;
  ROff_g_ad_k1 = -ROff_g_ad(n-1,:) / p.ROff_tau_ad;
  SOnOff_V_k1 = ( (p.SOnOff_E_L-SOnOff_V(n-1,:)) - SOnOff_R*SOnOff_g_ad(n-1,:).*(SOnOff_V(n-1,:)-p.SOnOff_E_k) - SOnOff_R*((((((SOnOff_On_PSC3_s(n-1,:)*SOnOff_On_PSC3_netcon).*(SOnOff_V(n-1,:)-p.SOnOff_On_PSC3_ESYN))*p.SOnOff_On_PSC3_gSYN)))+((((((SOnOff_Off_PSC3_s(n-1,:)*SOnOff_Off_PSC3_netcon).*(SOnOff_V(n-1,:)-p.SOnOff_Off_PSC3_ESYN))*p.SOnOff_Off_PSC3_gSYN))))) + SOnOff_R*p.SOnOff_Itonic.*SOnOff_Imask + SOnOff_R*p.SOnOff_noise.*randn(1,p.SOnOff_Npop) ) / SOnOff_tau;
  SOnOff_g_ad_k1 = -SOnOff_g_ad(n-1,:) / p.SOnOff_tau_ad;
  TD_V_k1 = ( (p.TD_E_L-TD_V(n-1,:)) - TD_R*TD_g_ad(n-1,:).*(TD_V(n-1,:)-p.TD_E_k) + TD_R*TD_Icur(k,:).*TD_Imask + TD_R*p.TD_noise.*randn(1,p.TD_Npop) ) / TD_tau;
  TD_g_ad_k1 = -TD_g_ad(n-1,:) / p.TD_tau_ad;
  X_V_k1 = ( (p.X_E_L-X_V(n-1,:)) - X_R*X_g_ad(n-1,:).*(X_V(n-1,:)-p.X_E_k) - X_R*((((((X_ROn_PSC3_s(n-1,:)*X_ROn_PSC3_netcon).*(X_V(n-1,:)-p.X_ROn_PSC3_ESYN))*p.X_ROn_PSC3_gSYN)))+((((((X_TD_PSC3_s(n-1,:)*X_TD_PSC3_netcon).*(X_V(n-1,:)-p.X_TD_PSC3_ESYN))*p.X_TD_PSC3_gSYN))))) + X_R*p.X_Itonic.*X_Imask + X_R*p.X_noise.*randn(1,p.X_Npop) ) / X_tau;
  X_g_ad_k1 = -X_g_ad(n-1,:) / p.X_tau_ad;
  C_V_k1 = ( (p.C_E_L-C_V(n-1,:)) - C_R*C_g_ad(n-1,:).*(C_V(n-1,:)-p.C_E_k) - C_R*((((((C_ROn_PSC3_s(n-1,:)*C_ROn_PSC3_netcon).*(C_V(n-1,:)-p.C_ROn_PSC3_ESYN))*p.C_ROn_PSC3_gSYN)))) + C_R*p.C_Itonic.*C_Imask + C_R*p.C_noise.*randn(1,p.C_Npop) ) / C_tau;
  C_g_ad_k1 = -C_g_ad(n-1,:) / p.C_tau_ad;
  ROn_On_PSC3_s_k1 = ( ROn_On_PSC3_scale * ROn_On_PSC3_x(n-1,:) - ROn_On_PSC3_s(n-1,:) )/p.ROn_On_PSC3_tauR;
  ROn_On_PSC3_x_k1 = -ROn_On_PSC3_x(n-1,:)/p.ROn_On_PSC3_tauD;
  ROn_On_PSC3_F_k1 = (1 - ROn_On_PSC3_F(n-1,:))/p.ROn_On_PSC3_tauF;
  ROn_On_PSC3_P_k1 = (1 - ROn_On_PSC3_P(n-1,:))/p.ROn_On_PSC3_tauP;
  ROn_On_PSC3_q_k1 = 0;
  SOnOff_On_PSC3_s_k1 = ( SOnOff_On_PSC3_scale * SOnOff_On_PSC3_x(n-1,:) - SOnOff_On_PSC3_s(n-1,:) )/p.SOnOff_On_PSC3_tauR;
  SOnOff_On_PSC3_x_k1 = -SOnOff_On_PSC3_x(n-1,:)/p.SOnOff_On_PSC3_tauD;
  SOnOff_On_PSC3_F_k1 = (1 - SOnOff_On_PSC3_F(n-1,:))/p.SOnOff_On_PSC3_tauF;
  SOnOff_On_PSC3_P_k1 = (1 - SOnOff_On_PSC3_P(n-1,:))/p.SOnOff_On_PSC3_tauP;
  SOnOff_On_PSC3_q_k1 = 0;
  ROn_SOnOff_PSC3_s_k1 = ( ROn_SOnOff_PSC3_scale * ROn_SOnOff_PSC3_x(n-1,:) - ROn_SOnOff_PSC3_s(n-1,:) )/p.ROn_SOnOff_PSC3_tauR;
  ROn_SOnOff_PSC3_x_k1 = -ROn_SOnOff_PSC3_x(n-1,:)/p.ROn_SOnOff_PSC3_tauD;
  ROn_SOnOff_PSC3_F_k1 = (1 - ROn_SOnOff_PSC3_F(n-1,:))/p.ROn_SOnOff_PSC3_tauF;
  ROn_SOnOff_PSC3_P_k1 = (1 - ROn_SOnOff_PSC3_P(n-1,:))/p.ROn_SOnOff_PSC3_tauP;
  ROn_SOnOff_PSC3_q_k1 = 0;
  ROff_SOnOff_PSC3_s_k1 = ( ROff_SOnOff_PSC3_scale * ROff_SOnOff_PSC3_x(n-1,:) - ROff_SOnOff_PSC3_s(n-1,:) )/p.ROff_SOnOff_PSC3_tauR;
  ROff_SOnOff_PSC3_x_k1 = -ROff_SOnOff_PSC3_x(n-1,:)/p.ROff_SOnOff_PSC3_tauD;
  ROff_SOnOff_PSC3_F_k1 = (1 - ROff_SOnOff_PSC3_F(n-1,:))/p.ROff_SOnOff_PSC3_tauF;
  ROff_SOnOff_PSC3_P_k1 = (1 - ROff_SOnOff_PSC3_P(n-1,:))/p.ROff_SOnOff_PSC3_tauP;
  ROff_SOnOff_PSC3_q_k1 = 0;
  ROff_Off_PSC3_s_k1 = ( ROff_Off_PSC3_scale * ROff_Off_PSC3_x(n-1,:) - ROff_Off_PSC3_s(n-1,:) )/p.ROff_Off_PSC3_tauR;
  ROff_Off_PSC3_x_k1 = -ROff_Off_PSC3_x(n-1,:)/p.ROff_Off_PSC3_tauD;
  ROff_Off_PSC3_F_k1 = (1 - ROff_Off_PSC3_F(n-1,:))/p.ROff_Off_PSC3_tauF;
  ROff_Off_PSC3_P_k1 = (1 - ROff_Off_PSC3_P(n-1,:))/p.ROff_Off_PSC3_tauP;
  ROff_Off_PSC3_q_k1 = 0;
  SOnOff_Off_PSC3_s_k1 = ( SOnOff_Off_PSC3_scale * SOnOff_Off_PSC3_x(n-1,:) - SOnOff_Off_PSC3_s(n-1,:) )/p.SOnOff_Off_PSC3_tauR;
  SOnOff_Off_PSC3_x_k1 = -SOnOff_Off_PSC3_x(n-1,:)/p.SOnOff_Off_PSC3_tauD;
  SOnOff_Off_PSC3_F_k1 = (1 - SOnOff_Off_PSC3_F(n-1,:))/p.SOnOff_Off_PSC3_tauF;
  SOnOff_Off_PSC3_P_k1 = (1 - SOnOff_Off_PSC3_P(n-1,:))/p.SOnOff_Off_PSC3_tauP;
  SOnOff_Off_PSC3_q_k1 = 0;
  ROn_ROn_iNoise_V3_sn_k1 = ( ROn_ROn_iNoise_V3_scale * ROn_ROn_iNoise_V3_xn(n-1,:) - ROn_ROn_iNoise_V3_sn(n-1,:) )/p.ROn_ROn_iNoise_V3_tauR_N;
  ROn_ROn_iNoise_V3_xn_k1 = -ROn_ROn_iNoise_V3_xn(n-1,:)/p.ROn_ROn_iNoise_V3_tauD_N + ROn_ROn_iNoise_V3_token(k,:)/p.ROn_ROn_iNoise_V3_dt;
  X_ROn_PSC3_s_k1 = ( X_ROn_PSC3_scale * X_ROn_PSC3_x(n-1,:) - X_ROn_PSC3_s(n-1,:) )/p.X_ROn_PSC3_tauR;
  X_ROn_PSC3_x_k1 = -X_ROn_PSC3_x(n-1,:)/p.X_ROn_PSC3_tauD;
  X_ROn_PSC3_F_k1 = (1 - X_ROn_PSC3_F(n-1,:))/p.X_ROn_PSC3_tauF;
  X_ROn_PSC3_P_k1 = (1 - X_ROn_PSC3_P(n-1,:))/p.X_ROn_PSC3_tauP;
  X_ROn_PSC3_q_k1 = 0;
  ROn_X_PSC3_s_k1 = ( ROn_X_PSC3_scale * ROn_X_PSC3_x(n-1,:) - ROn_X_PSC3_s(n-1,:) )/p.ROn_X_PSC3_tauR;
  ROn_X_PSC3_x_k1 = -ROn_X_PSC3_x(n-1,:)/p.ROn_X_PSC3_tauD;
  ROn_X_PSC3_F_k1 = (1 - ROn_X_PSC3_F(n-1,:))/p.ROn_X_PSC3_tauF;
  ROn_X_PSC3_P_k1 = (1 - ROn_X_PSC3_P(n-1,:))/p.ROn_X_PSC3_tauP;
  ROn_X_PSC3_q_k1 = 0;
  ROn_TD_PSC3_s_k1 = ( ROn_TD_PSC3_scale * ROn_TD_PSC3_x(n-1,:) - ROn_TD_PSC3_s(n-1,:) )/p.ROn_TD_PSC3_tauR;
  ROn_TD_PSC3_x_k1 = -ROn_TD_PSC3_x(n-1,:)/p.ROn_TD_PSC3_tauD;
  ROn_TD_PSC3_F_k1 = (1 - ROn_TD_PSC3_F(n-1,:))/p.ROn_TD_PSC3_tauF;
  ROn_TD_PSC3_P_k1 = (1 - ROn_TD_PSC3_P(n-1,:))/p.ROn_TD_PSC3_tauP;
  ROn_TD_PSC3_q_k1 = 0;
  ROff_TD_PSC3_s_k1 = ( ROff_TD_PSC3_scale * ROff_TD_PSC3_x(n-1,:) - ROff_TD_PSC3_s(n-1,:) )/p.ROff_TD_PSC3_tauR;
  ROff_TD_PSC3_x_k1 = -ROff_TD_PSC3_x(n-1,:)/p.ROff_TD_PSC3_tauD;
  ROff_TD_PSC3_F_k1 = (1 - ROff_TD_PSC3_F(n-1,:))/p.ROff_TD_PSC3_tauF;
  ROff_TD_PSC3_P_k1 = (1 - ROff_TD_PSC3_P(n-1,:))/p.ROff_TD_PSC3_tauP;
  ROff_TD_PSC3_q_k1 = 0;
  X_TD_PSC3_s_k1 = ( X_TD_PSC3_scale * X_TD_PSC3_x(n-1,:) - X_TD_PSC3_s(n-1,:) )/p.X_TD_PSC3_tauR;
  X_TD_PSC3_x_k1 = -X_TD_PSC3_x(n-1,:)/p.X_TD_PSC3_tauD;
  X_TD_PSC3_F_k1 = (1 - X_TD_PSC3_F(n-1,:))/p.X_TD_PSC3_tauF;
  X_TD_PSC3_P_k1 = (1 - X_TD_PSC3_P(n-1,:))/p.X_TD_PSC3_tauP;
  X_TD_PSC3_q_k1 = 0;
  C_ROn_PSC3_s_k1 = ( C_ROn_PSC3_scale * C_ROn_PSC3_x(n-1,:) - C_ROn_PSC3_s(n-1,:) )/p.C_ROn_PSC3_tauR;
  C_ROn_PSC3_x_k1 = -C_ROn_PSC3_x(n-1,:)/p.C_ROn_PSC3_tauD;
  C_ROn_PSC3_F_k1 = (1 - C_ROn_PSC3_F(n-1,:))/p.C_ROn_PSC3_tauF;
  C_ROn_PSC3_P_k1 = (1 - C_ROn_PSC3_P(n-1,:))/p.C_ROn_PSC3_tauP;
  C_ROn_PSC3_q_k1 = 0;

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  On_V(n,:) = On_V(n-1,:)+dt*On_V_k1;
  On_g_ad(n,:) = On_g_ad(n-1,:)+dt*On_g_ad_k1;
  Off_V(n,:) = Off_V(n-1,:)+dt*Off_V_k1;
  Off_g_ad(n,:) = Off_g_ad(n-1,:)+dt*Off_g_ad_k1;
  ROn_V(n,:) = ROn_V(n-1,:)+dt*ROn_V_k1;
  ROn_g_ad(n,:) = ROn_g_ad(n-1,:)+dt*ROn_g_ad_k1;
  ROff_V(n,:) = ROff_V(n-1,:)+dt*ROff_V_k1;
  ROff_g_ad(n,:) = ROff_g_ad(n-1,:)+dt*ROff_g_ad_k1;
  SOnOff_V(n,:) = SOnOff_V(n-1,:)+dt*SOnOff_V_k1;
  SOnOff_g_ad(n,:) = SOnOff_g_ad(n-1,:)+dt*SOnOff_g_ad_k1;
  TD_V(n,:) = TD_V(n-1,:)+dt*TD_V_k1;
  TD_g_ad(n,:) = TD_g_ad(n-1,:)+dt*TD_g_ad_k1;
  X_V(n,:) = X_V(n-1,:)+dt*X_V_k1;
  X_g_ad(n,:) = X_g_ad(n-1,:)+dt*X_g_ad_k1;
  C_V(n,:) = C_V(n-1,:)+dt*C_V_k1;
  C_g_ad(n,:) = C_g_ad(n-1,:)+dt*C_g_ad_k1;
  ROn_On_PSC3_s(n,:) = ROn_On_PSC3_s(n-1,:)+dt*ROn_On_PSC3_s_k1;
  ROn_On_PSC3_x(n,:) = ROn_On_PSC3_x(n-1,:)+dt*ROn_On_PSC3_x_k1;
  ROn_On_PSC3_F(n,:) = ROn_On_PSC3_F(n-1,:)+dt*ROn_On_PSC3_F_k1;
  ROn_On_PSC3_P(n,:) = ROn_On_PSC3_P(n-1,:)+dt*ROn_On_PSC3_P_k1;
  ROn_On_PSC3_q(n,:) = ROn_On_PSC3_q(n-1,:)+dt*ROn_On_PSC3_q_k1;
  SOnOff_On_PSC3_s(n,:) = SOnOff_On_PSC3_s(n-1,:)+dt*SOnOff_On_PSC3_s_k1;
  SOnOff_On_PSC3_x(n,:) = SOnOff_On_PSC3_x(n-1,:)+dt*SOnOff_On_PSC3_x_k1;
  SOnOff_On_PSC3_F(n,:) = SOnOff_On_PSC3_F(n-1,:)+dt*SOnOff_On_PSC3_F_k1;
  SOnOff_On_PSC3_P(n,:) = SOnOff_On_PSC3_P(n-1,:)+dt*SOnOff_On_PSC3_P_k1;
  SOnOff_On_PSC3_q(n,:) = SOnOff_On_PSC3_q(n-1,:)+dt*SOnOff_On_PSC3_q_k1;
  ROn_SOnOff_PSC3_s(n,:) = ROn_SOnOff_PSC3_s(n-1,:)+dt*ROn_SOnOff_PSC3_s_k1;
  ROn_SOnOff_PSC3_x(n,:) = ROn_SOnOff_PSC3_x(n-1,:)+dt*ROn_SOnOff_PSC3_x_k1;
  ROn_SOnOff_PSC3_F(n,:) = ROn_SOnOff_PSC3_F(n-1,:)+dt*ROn_SOnOff_PSC3_F_k1;
  ROn_SOnOff_PSC3_P(n,:) = ROn_SOnOff_PSC3_P(n-1,:)+dt*ROn_SOnOff_PSC3_P_k1;
  ROn_SOnOff_PSC3_q(n,:) = ROn_SOnOff_PSC3_q(n-1,:)+dt*ROn_SOnOff_PSC3_q_k1;
  ROff_SOnOff_PSC3_s(n,:) = ROff_SOnOff_PSC3_s(n-1,:)+dt*ROff_SOnOff_PSC3_s_k1;
  ROff_SOnOff_PSC3_x(n,:) = ROff_SOnOff_PSC3_x(n-1,:)+dt*ROff_SOnOff_PSC3_x_k1;
  ROff_SOnOff_PSC3_F(n,:) = ROff_SOnOff_PSC3_F(n-1,:)+dt*ROff_SOnOff_PSC3_F_k1;
  ROff_SOnOff_PSC3_P(n,:) = ROff_SOnOff_PSC3_P(n-1,:)+dt*ROff_SOnOff_PSC3_P_k1;
  ROff_SOnOff_PSC3_q(n,:) = ROff_SOnOff_PSC3_q(n-1,:)+dt*ROff_SOnOff_PSC3_q_k1;
  ROff_Off_PSC3_s(n,:) = ROff_Off_PSC3_s(n-1,:)+dt*ROff_Off_PSC3_s_k1;
  ROff_Off_PSC3_x(n,:) = ROff_Off_PSC3_x(n-1,:)+dt*ROff_Off_PSC3_x_k1;
  ROff_Off_PSC3_F(n,:) = ROff_Off_PSC3_F(n-1,:)+dt*ROff_Off_PSC3_F_k1;
  ROff_Off_PSC3_P(n,:) = ROff_Off_PSC3_P(n-1,:)+dt*ROff_Off_PSC3_P_k1;
  ROff_Off_PSC3_q(n,:) = ROff_Off_PSC3_q(n-1,:)+dt*ROff_Off_PSC3_q_k1;
  SOnOff_Off_PSC3_s(n,:) = SOnOff_Off_PSC3_s(n-1,:)+dt*SOnOff_Off_PSC3_s_k1;
  SOnOff_Off_PSC3_x(n,:) = SOnOff_Off_PSC3_x(n-1,:)+dt*SOnOff_Off_PSC3_x_k1;
  SOnOff_Off_PSC3_F(n,:) = SOnOff_Off_PSC3_F(n-1,:)+dt*SOnOff_Off_PSC3_F_k1;
  SOnOff_Off_PSC3_P(n,:) = SOnOff_Off_PSC3_P(n-1,:)+dt*SOnOff_Off_PSC3_P_k1;
  SOnOff_Off_PSC3_q(n,:) = SOnOff_Off_PSC3_q(n-1,:)+dt*SOnOff_Off_PSC3_q_k1;
  ROn_ROn_iNoise_V3_sn(n,:) = ROn_ROn_iNoise_V3_sn(n-1,:)+dt*ROn_ROn_iNoise_V3_sn_k1;
  ROn_ROn_iNoise_V3_xn(n,:) = ROn_ROn_iNoise_V3_xn(n-1,:)+dt*ROn_ROn_iNoise_V3_xn_k1;
  X_ROn_PSC3_s(n,:) = X_ROn_PSC3_s(n-1,:)+dt*X_ROn_PSC3_s_k1;
  X_ROn_PSC3_x(n,:) = X_ROn_PSC3_x(n-1,:)+dt*X_ROn_PSC3_x_k1;
  X_ROn_PSC3_F(n,:) = X_ROn_PSC3_F(n-1,:)+dt*X_ROn_PSC3_F_k1;
  X_ROn_PSC3_P(n,:) = X_ROn_PSC3_P(n-1,:)+dt*X_ROn_PSC3_P_k1;
  X_ROn_PSC3_q(n,:) = X_ROn_PSC3_q(n-1,:)+dt*X_ROn_PSC3_q_k1;
  ROn_X_PSC3_s(n,:) = ROn_X_PSC3_s(n-1,:)+dt*ROn_X_PSC3_s_k1;
  ROn_X_PSC3_x(n,:) = ROn_X_PSC3_x(n-1,:)+dt*ROn_X_PSC3_x_k1;
  ROn_X_PSC3_F(n,:) = ROn_X_PSC3_F(n-1,:)+dt*ROn_X_PSC3_F_k1;
  ROn_X_PSC3_P(n,:) = ROn_X_PSC3_P(n-1,:)+dt*ROn_X_PSC3_P_k1;
  ROn_X_PSC3_q(n,:) = ROn_X_PSC3_q(n-1,:)+dt*ROn_X_PSC3_q_k1;
  ROn_TD_PSC3_s(n,:) = ROn_TD_PSC3_s(n-1,:)+dt*ROn_TD_PSC3_s_k1;
  ROn_TD_PSC3_x(n,:) = ROn_TD_PSC3_x(n-1,:)+dt*ROn_TD_PSC3_x_k1;
  ROn_TD_PSC3_F(n,:) = ROn_TD_PSC3_F(n-1,:)+dt*ROn_TD_PSC3_F_k1;
  ROn_TD_PSC3_P(n,:) = ROn_TD_PSC3_P(n-1,:)+dt*ROn_TD_PSC3_P_k1;
  ROn_TD_PSC3_q(n,:) = ROn_TD_PSC3_q(n-1,:)+dt*ROn_TD_PSC3_q_k1;
  ROff_TD_PSC3_s(n,:) = ROff_TD_PSC3_s(n-1,:)+dt*ROff_TD_PSC3_s_k1;
  ROff_TD_PSC3_x(n,:) = ROff_TD_PSC3_x(n-1,:)+dt*ROff_TD_PSC3_x_k1;
  ROff_TD_PSC3_F(n,:) = ROff_TD_PSC3_F(n-1,:)+dt*ROff_TD_PSC3_F_k1;
  ROff_TD_PSC3_P(n,:) = ROff_TD_PSC3_P(n-1,:)+dt*ROff_TD_PSC3_P_k1;
  ROff_TD_PSC3_q(n,:) = ROff_TD_PSC3_q(n-1,:)+dt*ROff_TD_PSC3_q_k1;
  X_TD_PSC3_s(n,:) = X_TD_PSC3_s(n-1,:)+dt*X_TD_PSC3_s_k1;
  X_TD_PSC3_x(n,:) = X_TD_PSC3_x(n-1,:)+dt*X_TD_PSC3_x_k1;
  X_TD_PSC3_F(n,:) = X_TD_PSC3_F(n-1,:)+dt*X_TD_PSC3_F_k1;
  X_TD_PSC3_P(n,:) = X_TD_PSC3_P(n-1,:)+dt*X_TD_PSC3_P_k1;
  X_TD_PSC3_q(n,:) = X_TD_PSC3_q(n-1,:)+dt*X_TD_PSC3_q_k1;
  C_ROn_PSC3_s(n,:) = C_ROn_PSC3_s(n-1,:)+dt*C_ROn_PSC3_s_k1;
  C_ROn_PSC3_x(n,:) = C_ROn_PSC3_x(n-1,:)+dt*C_ROn_PSC3_x_k1;
  C_ROn_PSC3_F(n,:) = C_ROn_PSC3_F(n-1,:)+dt*C_ROn_PSC3_F_k1;
  C_ROn_PSC3_P(n,:) = C_ROn_PSC3_P(n-1,:)+dt*C_ROn_PSC3_P_k1;
  C_ROn_PSC3_q(n,:) = C_ROn_PSC3_q(n-1,:)+dt*C_ROn_PSC3_q_k1;

  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=any(C_V(n,:)>=p.C_V_thresh&C_V(n-1,:)<p.C_V_thresh);
  conditional_indx=(C_V(n,:)>=p.C_V_thresh&C_V(n-1,:)<p.C_V_thresh);
  if conditional_test, C_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); C_tspike(C_buffer_index(i),i)=t; C_buffer_index(i)=mod(-1+(C_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(X_V(n,:)>=p.X_V_thresh&X_V(n-1,:)<p.X_V_thresh);
  conditional_indx=(X_V(n,:)>=p.X_V_thresh&X_V(n-1,:)<p.X_V_thresh);
  if conditional_test, X_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); X_tspike(X_buffer_index(i),i)=t; X_buffer_index(i)=mod(-1+(X_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(TD_V(n,:)>=p.TD_V_thresh&TD_V(n-1,:)<p.TD_V_thresh);
  conditional_indx=(TD_V(n,:)>=p.TD_V_thresh&TD_V(n-1,:)<p.TD_V_thresh);
  if conditional_test, TD_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); TD_tspike(TD_buffer_index(i),i)=t; TD_buffer_index(i)=mod(-1+(TD_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(SOnOff_V(n,:)>=p.SOnOff_V_thresh&SOnOff_V(n-1,:)<p.SOnOff_V_thresh);
  conditional_indx=(SOnOff_V(n,:)>=p.SOnOff_V_thresh&SOnOff_V(n-1,:)<p.SOnOff_V_thresh);
  if conditional_test, SOnOff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); SOnOff_tspike(SOnOff_buffer_index(i),i)=t; SOnOff_buffer_index(i)=mod(-1+(SOnOff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(ROff_V(n,:)>=p.ROff_V_thresh&ROff_V(n-1,:)<p.ROff_V_thresh);
  conditional_indx=(ROff_V(n,:)>=p.ROff_V_thresh&ROff_V(n-1,:)<p.ROff_V_thresh);
  if conditional_test, ROff_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); ROff_tspike(ROff_buffer_index(i),i)=t; ROff_buffer_index(i)=mod(-1+(ROff_buffer_index(i)+1),5)+1; end; end
  conditional_test=any(ROn_V(n,:)>=p.ROn_V_thresh&ROn_V(n-1,:)<p.ROn_V_thresh);
  conditional_indx=(ROn_V(n,:)>=p.ROn_V_thresh&ROn_V(n-1,:)<p.ROn_V_thresh);
  if conditional_test, ROn_V_spikes(n,conditional_indx)=1;inds=find(conditional_indx); for j=1:length(inds), i=inds(j); ROn_tspike(ROn_buffer_index(i),i)=t; ROn_buffer_index(i)=mod(-1+(ROn_buffer_index(i)+1),5)+1; end; end
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
  conditional_test=any(ROn_V(n,:) > p.ROn_V_thresh);
  conditional_indx=(ROn_V(n,:) > p.ROn_V_thresh);
  if conditional_test, ROn_V(n,conditional_indx) = p.ROn_V_reset; ROn_g_ad(n,conditional_indx) = ROn_g_ad(n,conditional_indx) + p.ROn_g_inc; end
  conditional_test=any(any(t<=ROn_tspike+p.ROn_t_ref,1));
  conditional_indx=(any(t<=ROn_tspike+p.ROn_t_ref,1));
  if conditional_test, ROn_V(n,conditional_indx) = p.ROn_V_reset; end
  conditional_test=any(ROff_V(n,:) > p.ROff_V_thresh);
  conditional_indx=(ROff_V(n,:) > p.ROff_V_thresh);
  if conditional_test, ROff_V(n,conditional_indx) = p.ROff_V_reset; ROff_g_ad(n,conditional_indx) = ROff_g_ad(n,conditional_indx) + p.ROff_g_inc; end
  conditional_test=any(any(t<=ROff_tspike+p.ROff_t_ref,1));
  conditional_indx=(any(t<=ROff_tspike+p.ROff_t_ref,1));
  if conditional_test, ROff_V(n,conditional_indx) = p.ROff_V_reset; end
  conditional_test=any(SOnOff_V(n,:) > p.SOnOff_V_thresh);
  conditional_indx=(SOnOff_V(n,:) > p.SOnOff_V_thresh);
  if conditional_test, SOnOff_V(n,conditional_indx) = p.SOnOff_V_reset; SOnOff_g_ad(n,conditional_indx) = SOnOff_g_ad(n,conditional_indx) + p.SOnOff_g_inc; end
  conditional_test=any(any(t<=SOnOff_tspike+p.SOnOff_t_ref,1));
  conditional_indx=(any(t<=SOnOff_tspike+p.SOnOff_t_ref,1));
  if conditional_test, SOnOff_V(n,conditional_indx) = p.SOnOff_V_reset; end
  conditional_test=any(TD_V(n,:) > p.TD_V_thresh);
  conditional_indx=(TD_V(n,:) > p.TD_V_thresh);
  if conditional_test, TD_V(n,conditional_indx) = p.TD_V_reset; TD_g_ad(n,conditional_indx) = TD_g_ad(n,conditional_indx) + p.TD_g_inc; end
  conditional_test=any(any(t<=TD_tspike+p.TD_t_ref,1));
  conditional_indx=(any(t<=TD_tspike+p.TD_t_ref,1));
  if conditional_test, TD_V(n,conditional_indx) = p.TD_V_reset; end
  conditional_test=any(X_V(n,:) > p.X_V_thresh);
  conditional_indx=(X_V(n,:) > p.X_V_thresh);
  if conditional_test, X_V(n,conditional_indx) = p.X_V_reset; X_g_ad(n,conditional_indx) = X_g_ad(n,conditional_indx) + p.X_g_inc; end
  conditional_test=any(any(t<=X_tspike+p.X_t_ref,1));
  conditional_indx=(any(t<=X_tspike+p.X_t_ref,1));
  if conditional_test, X_V(n,conditional_indx) = p.X_V_reset; end
  conditional_test=any(C_V(n,:) > p.C_V_thresh);
  conditional_indx=(C_V(n,:) > p.C_V_thresh);
  if conditional_test, C_V(n,conditional_indx) = p.C_V_reset; C_g_ad(n,conditional_indx) = C_g_ad(n,conditional_indx) + p.C_g_inc; end
  conditional_test=any(any(t<=C_tspike+p.C_t_ref,1));
  conditional_indx=(any(t<=C_tspike+p.C_t_ref,1));
  if conditional_test, C_V(n,conditional_indx) = p.C_V_reset; end
  conditional_test=any(any(t == On_tspike+p.ROn_On_PSC3_delay,1));
  conditional_indx=(any(t == On_tspike+p.ROn_On_PSC3_delay,1));
  if conditional_test, ROn_On_PSC3_x(n,conditional_indx) = ROn_On_PSC3_x(n,conditional_indx) + ROn_On_PSC3_q(n,conditional_indx);ROn_On_PSC3_q(n,conditional_indx) = ROn_On_PSC3_F(n,conditional_indx).*ROn_On_PSC3_P(n,conditional_indx);ROn_On_PSC3_F(n,conditional_indx) = ROn_On_PSC3_F(n,conditional_indx) + p.ROn_On_PSC3_fF*(p.ROn_On_PSC3_maxF-ROn_On_PSC3_F(n,conditional_indx)); ROn_On_PSC3_P(n,conditional_indx) = ROn_On_PSC3_P(n,conditional_indx)*(1 - p.ROn_On_PSC3_fP); end
  conditional_test=any(any(t == On_tspike+p.SOnOff_On_PSC3_delay,1));
  conditional_indx=(any(t == On_tspike+p.SOnOff_On_PSC3_delay,1));
  if conditional_test, SOnOff_On_PSC3_x(n,conditional_indx) = SOnOff_On_PSC3_x(n,conditional_indx) + SOnOff_On_PSC3_q(n,conditional_indx);SOnOff_On_PSC3_q(n,conditional_indx) = SOnOff_On_PSC3_F(n,conditional_indx).*SOnOff_On_PSC3_P(n,conditional_indx);SOnOff_On_PSC3_F(n,conditional_indx) = SOnOff_On_PSC3_F(n,conditional_indx) + p.SOnOff_On_PSC3_fF*(p.SOnOff_On_PSC3_maxF-SOnOff_On_PSC3_F(n,conditional_indx)); SOnOff_On_PSC3_P(n,conditional_indx) = SOnOff_On_PSC3_P(n,conditional_indx)*(1 - p.SOnOff_On_PSC3_fP); end
  conditional_test=any(any(t == SOnOff_tspike+p.ROn_SOnOff_PSC3_delay,1));
  conditional_indx=(any(t == SOnOff_tspike+p.ROn_SOnOff_PSC3_delay,1));
  if conditional_test, ROn_SOnOff_PSC3_x(n,conditional_indx) = ROn_SOnOff_PSC3_x(n,conditional_indx) + ROn_SOnOff_PSC3_q(n,conditional_indx);ROn_SOnOff_PSC3_q(n,conditional_indx) = ROn_SOnOff_PSC3_F(n,conditional_indx).*ROn_SOnOff_PSC3_P(n,conditional_indx);ROn_SOnOff_PSC3_F(n,conditional_indx) = ROn_SOnOff_PSC3_F(n,conditional_indx) + p.ROn_SOnOff_PSC3_fF*(p.ROn_SOnOff_PSC3_maxF-ROn_SOnOff_PSC3_F(n,conditional_indx)); ROn_SOnOff_PSC3_P(n,conditional_indx) = ROn_SOnOff_PSC3_P(n,conditional_indx)*(1 - p.ROn_SOnOff_PSC3_fP); end
  conditional_test=any(any(t == SOnOff_tspike+p.ROff_SOnOff_PSC3_delay,1));
  conditional_indx=(any(t == SOnOff_tspike+p.ROff_SOnOff_PSC3_delay,1));
  if conditional_test, ROff_SOnOff_PSC3_x(n,conditional_indx) = ROff_SOnOff_PSC3_x(n,conditional_indx) + ROff_SOnOff_PSC3_q(n,conditional_indx);ROff_SOnOff_PSC3_q(n,conditional_indx) = ROff_SOnOff_PSC3_F(n,conditional_indx).*ROff_SOnOff_PSC3_P(n,conditional_indx);ROff_SOnOff_PSC3_F(n,conditional_indx) = ROff_SOnOff_PSC3_F(n,conditional_indx) + p.ROff_SOnOff_PSC3_fF*(p.ROff_SOnOff_PSC3_maxF-ROff_SOnOff_PSC3_F(n,conditional_indx)); ROff_SOnOff_PSC3_P(n,conditional_indx) = ROff_SOnOff_PSC3_P(n,conditional_indx)*(1 - p.ROff_SOnOff_PSC3_fP); end
  conditional_test=any(any(t == Off_tspike+p.ROff_Off_PSC3_delay,1));
  conditional_indx=(any(t == Off_tspike+p.ROff_Off_PSC3_delay,1));
  if conditional_test, ROff_Off_PSC3_x(n,conditional_indx) = ROff_Off_PSC3_x(n,conditional_indx) + ROff_Off_PSC3_q(n,conditional_indx);ROff_Off_PSC3_q(n,conditional_indx) = ROff_Off_PSC3_F(n,conditional_indx).*ROff_Off_PSC3_P(n,conditional_indx);ROff_Off_PSC3_F(n,conditional_indx) = ROff_Off_PSC3_F(n,conditional_indx) + p.ROff_Off_PSC3_fF*(p.ROff_Off_PSC3_maxF-ROff_Off_PSC3_F(n,conditional_indx)); ROff_Off_PSC3_P(n,conditional_indx) = ROff_Off_PSC3_P(n,conditional_indx)*(1 - p.ROff_Off_PSC3_fP); end
  conditional_test=any(any(t == Off_tspike+p.SOnOff_Off_PSC3_delay,1));
  conditional_indx=(any(t == Off_tspike+p.SOnOff_Off_PSC3_delay,1));
  if conditional_test, SOnOff_Off_PSC3_x(n,conditional_indx) = SOnOff_Off_PSC3_x(n,conditional_indx) + SOnOff_Off_PSC3_q(n,conditional_indx);SOnOff_Off_PSC3_q(n,conditional_indx) = SOnOff_Off_PSC3_F(n,conditional_indx).*SOnOff_Off_PSC3_P(n,conditional_indx);SOnOff_Off_PSC3_F(n,conditional_indx) = SOnOff_Off_PSC3_F(n,conditional_indx) + p.SOnOff_Off_PSC3_fF*(p.SOnOff_Off_PSC3_maxF-SOnOff_Off_PSC3_F(n,conditional_indx)); SOnOff_Off_PSC3_P(n,conditional_indx) = SOnOff_Off_PSC3_P(n,conditional_indx)*(1 - p.SOnOff_Off_PSC3_fP); end
  conditional_test=any(any(t == ROn_tspike+p.X_ROn_PSC3_delay,1));
  conditional_indx=(any(t == ROn_tspike+p.X_ROn_PSC3_delay,1));
  if conditional_test, X_ROn_PSC3_x(n,conditional_indx) = X_ROn_PSC3_x(n,conditional_indx) + X_ROn_PSC3_q(n,conditional_indx);X_ROn_PSC3_q(n,conditional_indx) = X_ROn_PSC3_F(n,conditional_indx).*X_ROn_PSC3_P(n,conditional_indx);X_ROn_PSC3_F(n,conditional_indx) = X_ROn_PSC3_F(n,conditional_indx) + p.X_ROn_PSC3_fF*(p.X_ROn_PSC3_maxF-X_ROn_PSC3_F(n,conditional_indx)); X_ROn_PSC3_P(n,conditional_indx) = X_ROn_PSC3_P(n,conditional_indx)*(1 - p.X_ROn_PSC3_fP); end
  conditional_test=any(any(t == X_tspike+p.ROn_X_PSC3_delay,1));
  conditional_indx=(any(t == X_tspike+p.ROn_X_PSC3_delay,1));
  if conditional_test, ROn_X_PSC3_x(n,conditional_indx) = ROn_X_PSC3_x(n,conditional_indx) + ROn_X_PSC3_q(n,conditional_indx);ROn_X_PSC3_q(n,conditional_indx) = ROn_X_PSC3_F(n,conditional_indx).*ROn_X_PSC3_P(n,conditional_indx);ROn_X_PSC3_F(n,conditional_indx) = ROn_X_PSC3_F(n,conditional_indx) + p.ROn_X_PSC3_fF*(p.ROn_X_PSC3_maxF-ROn_X_PSC3_F(n,conditional_indx)); ROn_X_PSC3_P(n,conditional_indx) = ROn_X_PSC3_P(n,conditional_indx)*(1 - p.ROn_X_PSC3_fP); end
  conditional_test=any(any(t == TD_tspike+p.ROn_TD_PSC3_delay,1));
  conditional_indx=(any(t == TD_tspike+p.ROn_TD_PSC3_delay,1));
  if conditional_test, ROn_TD_PSC3_x(n,conditional_indx) = ROn_TD_PSC3_x(n,conditional_indx) + ROn_TD_PSC3_q(n,conditional_indx);ROn_TD_PSC3_q(n,conditional_indx) = ROn_TD_PSC3_F(n,conditional_indx).*ROn_TD_PSC3_P(n,conditional_indx);ROn_TD_PSC3_F(n,conditional_indx) = ROn_TD_PSC3_F(n,conditional_indx) + p.ROn_TD_PSC3_fF*(p.ROn_TD_PSC3_maxF-ROn_TD_PSC3_F(n,conditional_indx)); ROn_TD_PSC3_P(n,conditional_indx) = ROn_TD_PSC3_P(n,conditional_indx)*(1 - p.ROn_TD_PSC3_fP); end
  conditional_test=any(any(t == TD_tspike+p.ROff_TD_PSC3_delay,1));
  conditional_indx=(any(t == TD_tspike+p.ROff_TD_PSC3_delay,1));
  if conditional_test, ROff_TD_PSC3_x(n,conditional_indx) = ROff_TD_PSC3_x(n,conditional_indx) + ROff_TD_PSC3_q(n,conditional_indx);ROff_TD_PSC3_q(n,conditional_indx) = ROff_TD_PSC3_F(n,conditional_indx).*ROff_TD_PSC3_P(n,conditional_indx);ROff_TD_PSC3_F(n,conditional_indx) = ROff_TD_PSC3_F(n,conditional_indx) + p.ROff_TD_PSC3_fF*(p.ROff_TD_PSC3_maxF-ROff_TD_PSC3_F(n,conditional_indx)); ROff_TD_PSC3_P(n,conditional_indx) = ROff_TD_PSC3_P(n,conditional_indx)*(1 - p.ROff_TD_PSC3_fP); end
  conditional_test=any(any(t == TD_tspike+p.X_TD_PSC3_delay,1));
  conditional_indx=(any(t == TD_tspike+p.X_TD_PSC3_delay,1));
  if conditional_test, X_TD_PSC3_x(n,conditional_indx) = X_TD_PSC3_x(n,conditional_indx) + X_TD_PSC3_q(n,conditional_indx);X_TD_PSC3_q(n,conditional_indx) = X_TD_PSC3_F(n,conditional_indx).*X_TD_PSC3_P(n,conditional_indx);X_TD_PSC3_F(n,conditional_indx) = X_TD_PSC3_F(n,conditional_indx) + p.X_TD_PSC3_fF*(p.X_TD_PSC3_maxF-X_TD_PSC3_F(n,conditional_indx)); X_TD_PSC3_P(n,conditional_indx) = X_TD_PSC3_P(n,conditional_indx)*(1 - p.X_TD_PSC3_fP); end
  conditional_test=any(any(t == ROn_tspike+p.C_ROn_PSC3_delay,1));
  conditional_indx=(any(t == ROn_tspike+p.C_ROn_PSC3_delay,1));
  if conditional_test, C_ROn_PSC3_x(n,conditional_indx) = C_ROn_PSC3_x(n,conditional_indx) + C_ROn_PSC3_q(n,conditional_indx);C_ROn_PSC3_q(n,conditional_indx) = C_ROn_PSC3_F(n,conditional_indx).*C_ROn_PSC3_P(n,conditional_indx);C_ROn_PSC3_F(n,conditional_indx) = C_ROn_PSC3_F(n,conditional_indx) + p.C_ROn_PSC3_fF*(p.C_ROn_PSC3_maxF-C_ROn_PSC3_F(n,conditional_indx)); C_ROn_PSC3_P(n,conditional_indx) = C_ROn_PSC3_P(n,conditional_indx)*(1 - p.C_ROn_PSC3_fP); end

  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  On_On_IC_iIC(n,:)=p.On_On_IC_g_postIC*(On_On_IC_input(k,:)*On_On_IC_netcon).*(On_V(n,:)-p.On_On_IC_E_exc);
  Off_Off_IC_iIC(n,:)=p.Off_Off_IC_g_postIC*(Off_Off_IC_input(k,:)*Off_Off_IC_netcon).*(Off_V(n,:)-p.Off_Off_IC_E_exc);
  ROn_On_PSC3_syn(n,:)=((ROn_On_PSC3_s(n,:)*ROn_On_PSC3_netcon).*(ROn_V(n,:)-p.ROn_On_PSC3_ESYN))*p.ROn_On_PSC3_gSYN;
  SOnOff_On_PSC3_syn(n,:)=((SOnOff_On_PSC3_s(n,:)*SOnOff_On_PSC3_netcon).*(SOnOff_V(n,:)-p.SOnOff_On_PSC3_ESYN))*p.SOnOff_On_PSC3_gSYN;
  ROn_SOnOff_PSC3_syn(n,:)=((ROn_SOnOff_PSC3_s(n,:)*ROn_SOnOff_PSC3_netcon).*(ROn_V(n,:)-p.ROn_SOnOff_PSC3_ESYN))*p.ROn_SOnOff_PSC3_gSYN;
  ROff_SOnOff_PSC3_syn(n,:)=((ROff_SOnOff_PSC3_s(n,:)*ROff_SOnOff_PSC3_netcon).*(ROff_V(n,:)-p.ROff_SOnOff_PSC3_ESYN))*p.ROff_SOnOff_PSC3_gSYN;
  ROff_Off_PSC3_syn(n,:)=((ROff_Off_PSC3_s(n,:)*ROff_Off_PSC3_netcon).*(ROff_V(n,:)-p.ROff_Off_PSC3_ESYN))*p.ROff_Off_PSC3_gSYN;
  SOnOff_Off_PSC3_syn(n,:)=((SOnOff_Off_PSC3_s(n,:)*SOnOff_Off_PSC3_netcon).*(SOnOff_V(n,:)-p.SOnOff_Off_PSC3_ESYN))*p.SOnOff_Off_PSC3_gSYN;
  X_ROn_PSC3_syn(n,:)=((X_ROn_PSC3_s(n,:)*X_ROn_PSC3_netcon).*(X_V(n,:)-p.X_ROn_PSC3_ESYN))*p.X_ROn_PSC3_gSYN;
  ROn_X_PSC3_syn(n,:)=((ROn_X_PSC3_s(n,:)*ROn_X_PSC3_netcon).*(ROn_V(n,:)-p.ROn_X_PSC3_ESYN))*p.ROn_X_PSC3_gSYN;
  ROn_TD_PSC3_syn(n,:)=((ROn_TD_PSC3_s(n,:)*ROn_TD_PSC3_netcon).*(ROn_V(n,:)-p.ROn_TD_PSC3_ESYN))*p.ROn_TD_PSC3_gSYN;
  ROff_TD_PSC3_syn(n,:)=((ROff_TD_PSC3_s(n,:)*ROff_TD_PSC3_netcon).*(ROff_V(n,:)-p.ROff_TD_PSC3_ESYN))*p.ROff_TD_PSC3_gSYN;
  X_TD_PSC3_syn(n,:)=((X_TD_PSC3_s(n,:)*X_TD_PSC3_netcon).*(X_V(n,:)-p.X_TD_PSC3_ESYN))*p.X_TD_PSC3_gSYN;
  C_ROn_PSC3_syn(n,:)=((C_ROn_PSC3_s(n,:)*C_ROn_PSC3_netcon).*(C_V(n,:)-p.C_ROn_PSC3_ESYN))*p.C_ROn_PSC3_gSYN;
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end
