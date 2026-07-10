%% This file provides a demo code for Section III.
% Created by Justin Lien on 2025/11/25.
% Implemented on Matlab 2024a.

clc; close all; clear all;

% To secure the folder path
restoredefaultpath 
addpath("../Helper_Functions/")

set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultTextFontName','Times New Roman')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultLegendFontName','Times New Roman')

%% Set up

%%% Define system parameters
n = 2; % Dimension
B = zeros(n,n,n); 
B(1,1,2) = -1;
B(2,1,1) = 1;
B = @(t) (1+0.2*cos(sin(2*pi*t))) * B; % Quadratic
A = @(t) (1+0.3*sin(2*pi*t))*[-1 +2; -1, -2]; % Linear 
C = @(t) [0.5;0]; % Constant

Q = @(t) [1, 0.5; 0.5, 1]; % Stochastic matrix, diffusion matrix, noise covariance matrix
sq2Q = sqrtm(2*Q(0)); 

Gam = 0.50; % Noise correlation time


%%% Numerics
TT = 10000; % Total duration
dt = 0.001; % Integration timestep
ttRatio = 10; 
Dt = dt*ttRatio; % Sampling timestep
N = round(1/Dt); % Sampling frequency

% Time grid
tp = Dt*(0:(N-1));

% Point-wise evalution (true answer)
B0 = zeros(n,n,n,N);
A0 = zeros(n,n,N);
C0 = zeros(n,N);
Q0 = zeros(n,n,N);

for j = 1:length(tp)
    t = tp(j);
    B0(:,:,:,j) = B(t);
    A0(:,:,j)   = A(t);
    C0(:,j)     = C(t);
    Q0(:,:,j)   = Q(t);
end


%% Forward integration: Generate synthetic data
clc; close all;

rng('twister') % For reproducability


% Generate the colored noise 
tic
eta = zeros(n,(TT+10)/dt+1);
for j = 1:n
    eta(j,:) = Generate_Colored_Noise(Gam,TT+10,dt);
end
toc

% Generate the state vector
tic
f = BAC2F(B,A,C); % Define dynamics by packing B, A, and C
O = zeros(n,1); 
Walls = @(t,X) O; % The artificial wall
xvec = Generate_CS_Colored_nLIM_timeseries(B,A,C,Q(0),eta,dt,ttRatio,TT,Walls);
toc


%% Inverse method: Reconstruct system parameters (B,A,C,Q,gamma)
clc; close all;

%%% Step 1: compute the statistics
[E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta] = CS_nLIM_Observation(xvec,Dt);

% Setup for Step 2
myopts.Method = "Phase_average"; % "Fixed_phase" or "Phase_average"
myopts.Method = "Fixed_phase";
myopts.FM = 2; % This is for phase_average
myopts.ShowLassoPath = 0; % Make a figure or not (0/1). When this is 1, make sure that "Gam_range = 0.5" or any other desired scalar. Also make sure using for-loop instead of parfor-loop.
myopts.ShowRuntime = 0; % Make a figure or not (0/1). 

% Prelocate variables for Step 5
Gam_range = 0.41:0.01:0.60;
% Gam_range = 0.5; % Uncomment this for reproducing Figures 1, 2, and S1.
err_range = nan(size(Gam_range));
Bnum_range = cell(size(Gam_range));
Anum_range = cell(size(Gam_range));
Cnum_range = cell(size(Gam_range));
Qnum_range = cell(size(Gam_range));
Qnum2_range = cell(size(Gam_range));

% parfor i = 1:length(Gam_range)
for i = 1:length(Gam_range)

    Gam = Gam_range(i);

    %%% Step 2: Sparse identification  
    % idx_nonzeroB = zeros(n^3,1);
    % idx_nonzeroB([2,5]) = 1; % Prescribe active nonlinear components (When uncomment this and the previous line, comment the following 2 lines.)
    idx_nonzeroB = Sparse_Identification(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,myopts) % Data-driven
    idx_nonzeroB = idx_nonzeroB(:); % Use active nonlinear components identified by Lasso regression
    idx_nonzeroA = ones(n^2,1);
    idx_nonzeroC = ones(n^1,1); 

    %%% Step 3: Solve the dynamical operators
    [~,Bnum,Anum,Cnum,~] = CS_Colored_nLIM_Dynamics(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzeroB,idx_nonzeroA,idx_nonzeroC,myopts);
    

    %%% Step 4: Solve the stochastic matrix.
    Qnum = CS_Colored_nLIM_Diffusion(Bnum,Anum,Cnum,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta);
    Qnum2 = CS_Colored_nLIM_Diffusion_Residual_Method(xvec,Bnum,Anum,Cnum,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,"Colored");
    
    % Computing the model-reproduced correlation function for Step 5.
    tic
    t2idx = @(t) mymod( round(N*t+0.5), N );
    O = zeros(n,1);
    Walls = @(t,X) O;
    Bnumfun = @(t) Bnum(:,:,:,t2idx(t));
    Anumfun = @(t) Anum(:,:,t2idx(t));
    Cnumfun = @(t) Cnum(:,t2idx(t));
    xvec_model = Generate_CS_Colored_nLIM_timeseries(Bnumfun,Anumfun,Cnumfun,Qnum,eta,dt,ttRatio,TT,Walls);
    toc
    
    % Loss function: Eq. (10)
    maxlag = 100; % Corresponding to 1 unit of time
    K_model = CorrelationFunction_CS_obs(xvec_model,Dt,maxlag);
    K_input = CorrelationFunction_CS_obs(xvec,Dt,maxlag);
    err_range(i) = rel_err(K_model,K_input);

    % Store variables
    Bnum_range{i} = Bnum;
    Anum_range{i} = Anum;
    Cnum_range{i} = Cnum;
    Qnum_range{i} = Qnum;
    Qnum2_range{i} = Qnum2;

end

%% Step 5: Determine the noise correlation time
[~,opt_gam] = min(err_range);

ShowLossFunction = 1;
if ShowLossFunction == 1
    figure
    plot(Gam_range,err_range); hold on;
    scatter(Gam_range(opt_gam),err_range(opt_gam))
    
    ax = gca;
    ax.YLim = [0.005,0.055];
    lgd = [];
    
    MakeGoodFigure4SyntheticData(ax,lgd,"Loss_Func",strcat(myopts.Method,"_Loss_Function"))
end

%% Display results

% Pick up the estimated parameters for CS-Colored-nLIM
Bnum = Bnum_range{opt_gam};
Anum = Anum_range{opt_gam};
Cnum = Cnum_range{opt_gam};
Qnum = Qnum_range{opt_gam};
Qnum2 = Qnum2_range{opt_gam};

disp("The active nonlinear components are identified by Lasso regression. The relative errors are: ")
str = strcat("B: ",  num2str(100*rel_err( Bnum, B0 ),'%3.2f'),"%");
disp(str)
str = strcat("A: ",  num2str(100*rel_err( Anum, A0 ),'%3.2f'),"%");
disp(str)
str = strcat("C: ",  num2str(100*rel_err( Cnum, C0 ),'%3.2f'),"%");
disp(str)
str = strcat("Q (By Eq. A8): ",  num2str(100*rel_err( Qnum, Q(0) ),'%3.2f'),"%");
disp(str)
str = strcat("Q2 (By residual approach): ",  num2str(100*rel_err( Qnum2, Q(0) ),'%3.2f'),"%");
disp(str)

