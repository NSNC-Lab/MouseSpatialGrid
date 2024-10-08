% z = struct;
% 
% z.On_E_L = 10;
% z.On_E_k = 10;
% z.On_On_IC_g_postIC = [1,1;1,1];
% z.On_On_IC_E_exc = [1,1;1,1];
% z.On_Itonic = 10;
% z.On_noise = 10;
% z.On_Npop = 2;
% 
% On_V = [1,1;1,1];
% On_R = [1,1;1,1];
% On_g_ad = [1,1;1,1];
% On_tau = 5;
% On_On_IC_input = [1,1;1,1];
% On_On_IC_netcon = [1,1;1,1];
% On_Imask = 10;
% 
% 
% tic
% 
% On_E_L = gpuArray(z.On_E_L);
% On_E_k = gpuArray(z.On_E_k);
% On_On_IC_g_postIC = gpuArray(z.On_On_IC_g_postIC);
% On_On_IC_E_exc = gpuArray(z.On_On_IC_E_exc);
% On_Itonic = gpuArray(z.On_Itonic);
% On_noise = gpuArray(z.On_noise);
% On_Npop = gpuArray(z.On_Npop);
% 
% On_V2 = gpuArray([1,1;1,1]);
% On_R2 = gpuArray([1,1;1,1]);
% On_g_ad2 = gpuArray([1,1;1,1]);
% On_tau2 = gpuArray(5);
% On_On_IC_input2 = gpuArray([1,1;1,1]);
% On_On_IC_netcon2 = gpuArray([1,1;1,1]);
% On_Imask2 = gpuArray(10);
% On_V_k12 = gpuArray([0,0;0,0]);
% 
% toc;
% 
% tic;
% 
% for trial = 1:(168000)
% 
% 
% On_V_k12 = ( (On_E_L) - On_R2*On_g_ad2(:,:).*(On_V2(:,:)-On_E_k) -...
%     On_R2*((((On_On_IC_g_postIC*(On_On_IC_input2(:,:)*On_On_IC_netcon2).*(On_V2(:,:)-On_On_IC_E_exc)))))...
%     + On_R2*On_Itonic.*On_Imask2 + On_R2*On_noise.*randn(On_Npop,On_Npop) ) / On_tau2;
% 
% end
% 
% toc
% 
% tic
% 
% for trial = 1:(168000)
% 
% On_V_k1 = ( (z.On_E_L-On_V(:,:)) - On_R*On_g_ad(:,:).*(On_V(:,:)-z.On_E_k) -...
%     On_R*((((z.On_On_IC_g_postIC*(On_On_IC_input(:,:)*On_On_IC_netcon).*(On_V(:,:)-z.On_On_IC_E_exc)))))...
%     + On_R*z.On_Itonic.*On_Imask + On_R*z.On_noise.*randn(z.On_Npop,z.On_Npop) ) / On_tau;
% 
% end
% 
% toc
z= struct;

z.a = 1;
z.b = 2;

thingy = fieldnames(z)

thingy{1}

