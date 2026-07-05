%% Numerical 
clc; close all; clear all;


%　rmpath('../CS Linear Inverse Model/')
Type = "Colored";


% Noise correlation time
Gam = 0.5;

% Mean state of dynamics
Amean = 1.2*[ -1.0, 0.5, 1.0;
               0.5, -2.0, 0.0;
               0.0, 1.0, -1.0];

% Mean state of white-noise diffusion matrix
Dmean = 0.3 * [ 1.0, 0.0, 0.5; 
                0.0, 2.0, 0.0;
                0.5, 0.0, 1.0];

% Colored-noise diffusion matrix
Q = 0.5 * [ 1.0, 0.5, 0.0; 
            0.5, 1.0, 0.0;
            0.0, 0.0, 1.0];

n = 3; % Dimension
N = 100; % Resolution
dt = 1/N; % Timestep
tp = (0:(N-1))/N; % Time coordinate

A = zeros(n,n,N);
D = zeros(n,n,N);


% Enfore periodicity
tCp = dt*(0:(N-1));
for k = 1:N

    w = ones(n);
    w(1,:) = 1+0.2*sin(2*pi*tCp(k));
    w(2,:) = 1+0.4*sin(2*pi*tCp(k)+0.4);
    w(3,3) = 1+0.5*sin(2*pi*tCp(k)+0.8);
    A(:,:,k) = w.*Amean;

    w = ones(n);
    w(1,1) = 1+0.3*cos(2*pi*tCp(k));
    D(:,:,k) = w.*Dmean;
    
end

% Make system parameters functions
At = mat2fun(A,tp);
Dt = mat2fun(D,tp);


tt_Ratio = 1;
dt_integration = dt/tt_Ratio; % Integration timestep


%%% Simulation
rng('twister')
myopt = struct();
xvec = CS_Generate_timeseries(Type,At,Q,Dt,Gam,10000,dt,myopt);
xvec = xvec(1:n,1:tt_Ratio:end);

%%
clc;

% Analytic correlation function
maxlag = N;

[K,h_K] = CS_CorrFcn(Type,dt,maxlag,At,Q,Dt,Gam,tp);

% Observed correlation function
K_CS_obs = CorrFcn_CS_Obs(xvec,dt,maxlag);

disp('Relative error in lagged covariance function')
rel_err(K_CS_obs,K)


% Smooth the correlation function
myopt = struct();
myopt.fs = N;
myopt.fpass = 0.05*N;
K_CS_obs = Phase_Averaging(K_CS_obs,myopt);

K_CS_obs = K; % What happenes if we have infinite observation? The estimation error mainly arises from numerical errors

% Observed covariance
C_CS_obs = squeeze(K_CS_obs(:,:,:,1));

% Derivatives
N0 = K_CS_obs(:,:,:,1);
N1 = ( -K_CS_obs(:,:,:,1) + K_CS_obs(:,:,:,2) ) / dt;
N2 = ( K_CS_obs(:,:,:,1) -2*K_CS_obs(:,:,:,2) + K_CS_obs(:,:,:,3) ) / dt^2;

%% Linear inverse modeling 
clc;

myopt = struct();
myopt.Solve_Diffusion = "Analytic";
myopt.FM = 1;

myopt.A = A;
myopt.Q = Q;

[ALIM,QLIM,DLIM,h] = CS_LIM(Type,N0,N1,N2,Gam,myopt);


disp('Relative error in estimate')
rel_err(ALIM,A)
rel_err(QLIM,Q)
rel_err(DLIM,D)

%%

r1 = pagemtimes(sqrtm(2*Q),h_K.C(n+1:end,1:n,:));
r2 = pagemtimes(Q,pagetranspose(h.B));

rel_err(r1,r2)

r1 = pagemtimes(sqrtm(2*Q),h_K.C(n+1:end,1:n,:));
r2 = pagemtimes(Q,h.B);

rel_err(r1,r2)

% %%
% 
% close all;
% 
% k = 0;
% for i = 1:n
%     for j = 1:n
%         k = k + 1;
%         subplot(n,n,k)
%         v = A(i,j,:);
%         plot(tp,v(:)); hold on;
%         v = ALIM(i,j,:);
%         plot(tp,v(:)); hold on;
%     end
% end
% 
% figure
% k = 0;
% for i = 1:n
%     for j = 1:n
%         k = k + 1;
%         subplot(n,n,k)
%         v = N0(i,j,:);
%         plot(v(:)); hold on;
%         v = K(i,j,:,1);
%         plot(v(:)); hold on;
%     end
% end
% 
% figure
% k = 0;
% for i = 1:n
%     for j = 1:n
%         k = k + 1;
%         subplot(n,n,k)
%         v = DLIM(i,j,:);
%         plot(v(:)); hold on;
%         v = D(i,j,:);
%         plot(v(:)); hold on;
%     end
% end
% 
% 
