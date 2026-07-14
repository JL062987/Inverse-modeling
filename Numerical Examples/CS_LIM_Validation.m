%% Numerical 
clc; close all; clear all;

restoredefaultpath 
addpath("../Function/")

Type = "Colored"

% Noise correlation time
Gam = 0.2;

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


% Enforce periodicity
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


tt_Ratio = 10;
dt_integration = dt/tt_Ratio; % Integration timestep


%%% Simulation
rng('twister')
myopt = struct();
xvec = CS_Generate_timeseries(Type,At,Q,Dt,Gam,10000,dt_integration,myopt);
xvec = xvec(1:n,1:tt_Ratio:end); % Sparse observation

%%
clc;

% Analytic correlation function
maxlag = N;

[K,h_K] = CS_CorrFcn(Type,dt,maxlag,At,Q,Dt,Gam,tp);

% Observed correlation function
K_CS_obs = CorrFcn_CS_Obs(xvec,dt,maxlag);

disp('Relative error in lagged covariance function')
rel_err(K_CS_obs,K)

% K_CS_obs = K; % What happenes if we have infinite observation? The estimation error mainly arises from numerical errors

% Observed lagged covariance and its derivatives
N0 = squeeze( K_CS_obs(:,:,:,1) );
N1 = squeeze( -K_CS_obs(:,:,:,1) + K_CS_obs(:,:,:,2) ) / dt;
N2 = squeeze( K_CS_obs(:,:,:,1) -2*K_CS_obs(:,:,:,2) + K_CS_obs(:,:,:,3) ) / dt^2;

% Smooth the statistics
myopt = struct();
myopt.fs = N;
myopt.fpass = 0.05*N;

N0 = Phase_Averaging(N0,myopt);
N1 = Phase_Averaging(N1,myopt);
N2 = Phase_Averaging(N2,myopt);

%% Linear inverse modeling 
clc;

myopt = struct();
myopt.Solve_Diffusion = "Analytic";
myopt.FM = 1;

[ALIM,QLIM,DLIM,h] = CS_LIM(Type,N0,N1,N2,Gam,myopt);

disp("The relative errors are: ")
str = strcat("A: ",  num2str(100*rel_err(ALIM,A),'%3.2f'),"%");
disp(str)
str = strcat("Q: ",  num2str(100*rel_err(QLIM,Q),'%3.2f'),"%");
disp(str)
str = strcat("D: ",  num2str(100*rel_err(DLIM,D),'%3.2f'),"%");
disp(str)



