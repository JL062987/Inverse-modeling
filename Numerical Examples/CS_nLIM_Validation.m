% Created by Justin Lien on 2026/07/14.
% Implemented on Matlab 2024a.

clc; close all; clear all;

% To secure the folder path
restoredefaultpath 
addpath("../Function/")

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


% Active nonlinear components
idx_nonzeroB = zeros(n^3,1);
idx_nonzeroB([2,5]) = 1;

Type = "Colored";
% Type = "White";

%%% Numerics
TT = 5000; % Total duration
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


% Generate the state vector
tic
rng('twister') % For reproducability
opt_sim = struct();
opt_sim.NoiseGeneration = "precompute";
f = BAC2F(B,A,C); % Define dynamics by packing B, A, and C
xvec = CS_Nonlinear_Generate_timeseries(Type,f,Q(0),Gam,TT,dt,opt_sim);
xvec = xvec(:,1:ttRatio:end); % Make sparse observation
toc


%% Inverse method: Reconstruct system parameters (B,A,C,Q,gamma)
clc; close all;

%%% Compute correlation functions
[Eobs,Kobs,Mobs,Sobs] = CorrFcn_CS_Obs(xvec,Dt,3,true);

%%% Lag-0 correlations and their derivatives (you may use advanced numerical differentiation)
E0 = Eobs;

K0 = Kobs(:,:,:,1);
K1 = ( -1*Kobs(:,:,:,1) + 1*Kobs(:,:,:,2) ) / Dt ;
K2 = ( +1*Kobs(:,:,:,1) + -2*Kobs(:,:,:,2) + 1*Kobs(:,:,:,3) ) / Dt^2 ;

M0 = Mobs(:,:,:,:,1);
M1 = ( -1*Mobs(:,:,:,:,1) + 1*Mobs(:,:,:,:,2) ) / Dt ;
M2 = ( +1*Mobs(:,:,:,:,1) + -2*Mobs(:,:,:,:,2) + 1*Mobs(:,:,:,:,3) ) / Dt^2 ;

S0 = Sobs(:,:,:,:,:,1);
S1 = ( -1*Sobs(:,:,:,:,:,1) + 1*Sobs(:,:,:,:,:,2) ) / Dt ;

% Phase averaging (Optional; you may define your own phase average process, or not apply it)
myopt_PA = struct();
myopt_PA.fs = N;
myopt_PA.fpass = 0.05*N;

% First-order
myopt_PA.dim = 2;
E0 = Phase_Averaging(E0,myopt_PA);

% Second-order
myopt_PA.dim = 3;
K0 = Phase_Averaging(K0,myopt_PA);
K1 = Phase_Averaging(K1,myopt_PA);
K2 = Phase_Averaging(K2,myopt_PA);

% Third-order
myopt_PA.dim = 4;
M0 = Phase_Averaging(M0,myopt_PA);
M1 = Phase_Averaging(M1,myopt_PA);
M2 = Phase_Averaging(M2,myopt_PA);

% Fourth-order
myopt_PA.dim = 5;
S0 = Phase_Averaging(S0,myopt_PA);
S1 = Phase_Averaging(S1,myopt_PA);

% Derivative with respect to phase
Etheta = (circshift(E0,-1,2) - E0) / Dt;
Ktheta = (circshift(K0,-1,3) - K0) / Dt;


% Prelocate variables for Step 5
if Type == "White"
    Gam_range = 0; % If noise type is white, uncomment this.
else
    Gam_range = 0.455:0.005:0.55;
end
err_range = nan(size(Gam_range));
Bnum_range = cell(size(Gam_range));
Anum_range = cell(size(Gam_range));
Cnum_range = cell(size(Gam_range));
Qnum_range = cell(size(Gam_range));


% Options for parameter estimation
opt_nLIM = struct();
opt_nLIM.Sparse_Identification.Type = "Empirical"; % "Empirical" or "Prescribed"
opt_nLIM.Sparse_Identification.Solve_Lambda = "Grid_Search"; % "Grid_Search" or a function handle
opt_nLIM.Sparse_Identification.user_opts = struct();
opt_nLIM.Sparse_Identification.idx_nonzeroB = idx_nonzeroB;
opt.Sparse_Identification.FISTA_opts.Display_Details = false;

opt_nLIM.Fourier_Mode.Type = "Fixed_phase"; % "Fixed_phase" or "Phase_average"
opt_nLIM.Fourier_Mode.FM = 2;

opt_nLIM.Fourier_Mode.Solve_Diffusion = "Analytic";

opt_nLIM.xvec = xvec;

% Options for simulations
opt_sim = struct();
opt_sim.NoiseGeneration = "augmented"; % augmented/precompute;

parfor i = 1:length(Gam_range)
% for i = 1:length(Gam_range)

    Gam = Gam_range(i);
   
    %%% Estimating system parameters
    [Bnum,Anum,Cnum,Qnum,h] = CS_nLIM(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,Ktheta,opt_nLIM);
    
    %%% Computing the model-reproduced correlation function for Step 5.
    t2idx = @(t) mymod( round(N*t+0.5), N );
    O = zeros(n,1);
    Walls = @(t,X) O;
    Bnumfun = @(t) Bnum(:,:,:,t2idx(t));
    Anumfun = @(t) Anum(:,:,t2idx(t));
    Cnumfun = @(t) Cnum(:,t2idx(t));
    Fnum = BAC2F(Bnumfun,Anumfun,Cnumfun,Walls); % Define dynamics
    
    rng(1000,'twister')
    tic
    xvec_model = CS_Nonlinear_Generate_timeseries(Type,Fnum,Qnum,Gam,TT,dt,opt_sim);
    xvec_model = xvec_model(:,1:ttRatio:end); % Make sparse observation
    toc
    
    % Loss function: Eq. (10)
    maxlag = 100; % Corresponding to 1 unit of time
    K_model = CorrFcn_CS_Obs(xvec_model,Dt,maxlag);
    K_input = CorrFcn_CS_Obs(xvec,Dt,maxlag);
    err_range(i) = rel_err(K_model,K_input);

    % Store variables
    Bnum_range{i} = Bnum;
    Anum_range{i} = Anum;
    Cnum_range{i} = Cnum;
    Qnum_range{i} = Qnum;

end

%%% Determine the noise correlation time
[~,opt_gam] = min(err_range);

%% Display results

% Pick up the estimated parameters for CS-Colored-nLIM
Bnum = Bnum_range{opt_gam};
Anum = Anum_range{opt_gam};
Cnum = Cnum_range{opt_gam};
Qnum = Qnum_range{opt_gam};

disp("The relative errors are: ")
str = strcat("B: ",  num2str(100*rel_err( Bnum, B0 ),'%3.2f'),"%");
disp(str)
str = strcat("A: ",  num2str(100*rel_err( Anum, A0 ),'%3.2f'),"%");
disp(str)
str = strcat("C: ",  num2str(100*rel_err( Cnum, C0 ),'%3.2f'),"%");
disp(str)
str = strcat("Q: ",  num2str(100*rel_err( Qnum, Q(0) ),'%3.2f'),"%");
disp(str)


%% Make figures
close all;

ShowLossFunction = 1;
if ShowLossFunction == 1
    figure
    plot(Gam_range,err_range); hold on;
    scatter(Gam_range(opt_gam),err_range(opt_gam))
    
    ax = gca;
    if opt_sim.NoiseGeneration == "augmented"
        ax.YLim = [0.038,0.042];
        ax.YAxis.TickValues = 0.038:0.002:0.042;
    else
        ax.YLim = [0.0147,0.017];
        ax.YAxis.TickValues = 0.015:0.001:0.017;
    end
    lgd = [];
    
    MakeGoodFigure4SyntheticData(ax,lgd,"Loss_Func",strcat(opt_nLIM.Fourier_Mode.Type,"_Loss_Function"))
end

ShowLassoPath = 0;
if ShowLassoPath == 1
    
    figure

    P = plot(h_SI.Lambda_range,h_SI.Lasso_path');
    ax = gca;
    
    lgd = legend('$B_{1,1,1}$','$B_{2,1,1}$','$B_{1,1,2}$', ...
        '$B_{2,1,2}$','$B_{1,2,2}$','$B_{2,2,2}$');
    
    if opt_nLIM.Fourier_Mode.Type == "Phase_average"
        set(gca, 'XScale', 'log'); 
        ax.YLim = [1.3e-3, 1.1e+2];
    end
    set(gca, 'YScale', 'log'); 
    
    MakeGoodFigure4SyntheticData(ax,lgd,"Lasso_Path",strcat(opt_nLIM.Fourier_Mode.Type,"_Lasso_Path"))
end

ShowRuntime = 0;
if ShowRuntime == 1
    figure

    H = bar(h_SI.Lambda_range,h_SI.iter_count);

    yyaxis right
    P = plot(h_SI.Lambda_range,h_SI.t_count);
    P.LineWidth = 1;
    P.Marker = 'o';
    P.MarkerSize = 4;

    ax = gca;

    lgd = legend("Iteration","Runtime");

    MakeGoodFigure4SyntheticData(ax,lgd,"Runtime",strcat(opt_nLIM.Fourier_Mode.Type,"_Runtime"))

end

%%

function MakeGoodFigure4SyntheticData(ax,Lgd,Type,Filename)

    x0 = 0.105; y0 = 0.140;
    dx = 0.862; dy = 0.818;
    myFigSize = [1, 1, 3.40, 2.0];
    myfs = 6.5;
    myFS = 8;

    set(gcf, 'Units', 'inch', 'InnerPosition', myFigSize);
    
    fs = myfs; FS = myFS;
    
    ax.Position = [x0,y0,dx,dy];

    ax.XAxis.FontSize = fs;
    ax.YAxis(1).FontSize = fs;

    if Type == "Lasso_Path"
        myXLabel = "Lasso parameter $\lambda$";
        myYLabel = "Vector norm $\left\Vert \mathrm{vec} B_{ijk} \right\Vert$";
    elseif Type == "Loss_Func"
        myXLabel = "Noise correlation time $\gamma$";
        myYLabel = "Loss function";
    elseif Type == "Runtime"
        ax.YAxis(2).FontSize = fs;
        myXLabel = "Lasso parameter $\lambda$";
        myYLabel = "Iteration";
        myYLabel2 = "Runtime (ms)";
        x0 = 0.110; y0 = 0.140;
        dx = 0.772; dy = 0.818;
        ax.Position = [x0,y0,dx,dy];
    end

    if Type == "Lasso_Path"
        ax.XAxis.TickLabels = arrayfun(@(a)strcat(num2str(a,'%3.1f')),ax.XAxis.TickValues,'uni',0);
    
        Lgd.ItemTokenSize = [7,7];
        Lgd.Box = 'off';
        Lgd.FontSize = fs;
        Lgd.NumColumns = 3;
        Lgd.Position(1) = x0+0.75*dx-0.5*Lgd.Position(3);
        Lgd.Position(2) = 0.8;

    elseif Type == "Loss_Func"

        


        ax.XAxis.TickLabels = arrayfun(@(a)strcat(num2str(a,'%3.2f')),ax.XAxis.TickValues,'uni',0);
        ax.YAxis.TickLabels = arrayfun(@(a)strcat(num2str(100*a,'%3.1f'),"\%"),ax.YAxis.TickValues,'uni',0);

    elseif Type == "Runtime"

        ax.XAxis.TickLabels = arrayfun(@(a)strcat(num2str(a,'%3.1f')),ax.XAxis.TickValues,'uni',0);
        ax.YAxis(2).Color = [0 0 0];
        ax.YAxis(2).TickLabels = arrayfun(@(a)strcat(num2str(1000*a,'%3.0f')),ax.YAxis(2).TickValues,'uni',0);

        Lgd.ItemTokenSize = [8,8];
        Lgd.Box = 'off';
        Lgd.FontSize = fs;
        Lgd.NumColumns = 2;
        Lgd.Position(1) = x0+0.75*dx-0.5*Lgd.Position(3);
        Lgd.Position(2) = 0.8;
  
    end
    
    ax = axes('Position', [0 0 1 1], 'Visible', 'off');    
    t1 = text(x0+0.5*dx,0, myXLabel, 'HorizontalAlignment', 'center', 'Parent', ax);
    t1.VerticalAlignment = 'bottom';
    t1.FontSize = FS;

    t1 = text(0, y0+0.5*dy, myYLabel, 'HorizontalAlignment', 'center', 'Parent', ax);
    t1.VerticalAlignment = 'top';
    t1.Rotation = 90;
    t1.FontSize = FS;

    if Type == "Runtime" || Type == "Runtime2"
        t1 = text(1, y0+0.5*dy, myYLabel2, 'HorizontalAlignment', 'center', 'Parent', ax);
        t1.VerticalAlignment = 'bottom';
        t1.Rotation = 90;
        t1.FontSize = FS;
    end

    % mysavegcf(Filename)


end
