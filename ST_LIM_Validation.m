clc; close all; clear all;

rng("twister") % For reproducability


% System dimension
n = 3;

% Dynamical matrix
A = [-1 1/2 1; 1/2 -2 0; 0 1 -1]; 
eig(A);

% Diffusion matrix for colored noise
Q = 0.5*[1 0 1/2; 0 2 1 ; 1/2 1 2];
eig(Q);

% Diffusion matrix for white noise
D = 0.5*eye(3);
eig(D);

% Select the LIM configuration for validation
Type = "CW";

% Noise menory
Gam = 0.5;



%% Simulation

Tf = 10000; % Timespan
Delta_t = 0.001; % Integration timestep
tt_ratio = 100;
dt = tt_ratio*Delta_t; % Sampling timestep
maxlag = round(1/(dt)); % The lag specifying the interval over which the minimization of the Gam-selection algorithm is taken.

xvec = ST_Generate_timeseries(Type,A,Q,D,Gam,Tf,Delta_t);
xvec = xvec(:,1:tt_ratio:end); % Make a sparse observation
    
%% LIM implementation
clc; close all;

% The observed correlation function up to maxlag.
Corr_Obs = CorrFcn_ST_Obs(xvec,maxlag);

% Covariance
N0 = Corr_Obs(:,:,1);

% Try a lower-order scheme
N1 = (Corr_Obs(:,:,2)-Corr_Obs(:,:,1))/dt;
N2 = (Corr_Obs(:,:,3)-2*Corr_Obs(:,:,2)+Corr_Obs(:,:,1))/dt^2;

% Try a higher-order scheme
N1 = (-0.5*Corr_Obs(:,:,3)+2*Corr_Obs(:,:,2)-1.5*Corr_Obs(:,:,1))/dt;
N2 = (-1*Corr_Obs(:,:,4)+4*Corr_Obs(:,:,3)-5*Corr_Obs(:,:,2)+2*Corr_Obs(:,:,1))/dt^2;


[A_Num,Q_Num,D_Num] = ST_LIM(Type,N0,N1,N2,Gam);
KLIM = ST_CorrFcn(Type,dt,maxlag,A_Num,Q_Num,D_Num,Gam);

disp('Relative error of the LIM estimates')
rel_err(A_Num,A)
rel_err(Q_Num,Q)
rel_err(D_Num,D)


%%% True lagged covariance and its derivatives


K = ST_CorrFcn(Type,dt,maxlag,A,Q,D,Gam);

if Type == "White"; B = 0; else; B = inv(eye(n)-Gam*A); end

C = K(:,:,1);
N0_true = C;
N1_true = A*C+Q*B';
N2_true = (-N1_true+A*(Gam*N1_true+C))/Gam;

disp('Relative error of the derivative estimates')
rel_err(N0,N0_true)
rel_err(N1,N1_true)
rel_err(N2,N2_true)


disp('Relative error of the observed and true lagged covariance functions')
rel_err(Corr_Obs,K)




%% Display lagged covariance functions


k = 0;
for i = 1:n
    for j = 1:n
        k = k + 1;
        subplot(n,n,k);
        v = K(i,j,:);
        plot(v(:)); hold on;
        v = Corr_Obs(i,j,:);
        plot(v(:)); hold on;
        v = KLIM(i,j,:);
        plot(v(:)); hold on;
    end
end