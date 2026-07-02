%% This is for the Figure 4.2 in my PhD thesis.
% Please refer to section 4.2.5.
clc; close all; clear all;

set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultTextFontName','Times New Roman')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultLegendFontName','Times New Roman')


%% Demostration of the basic use of Colored-LIM

% Dynamical matrix
A = [-1 1/2 1; 1/2 -2 0; 0 1 -1]; 
eig(A);

% Diffusion matrix
Q = [1 0 1/2; 0 2 1 ; 1/2 1 2];
eig(Q);

% Noise menory
Gam = 0.5;

% System dimension
n = 3;

B = inv(eye(n)-Gam*A);

Tf = 1000; % Timespan
Delta_t = 0.001; % Integration timestep
tt_ratio = 10;
dt = tt_ratio*Delta_t; % Sampling timestep
maxlag = round(1/(dt)); % The lag specifying the interval over which the minimization of the Gam-selection algorithm is taken.

% True correlation function
Corr_Colored_True = ST_CorrFcn('Colored',dt,maxlag,A,Q,[],Gam);
C = Corr_Colored_True(:,:,1);
N1_true = A*C+Q*B';
N2_true = (-N1_true+A*(Gam*N1_true+C))/Gam;

rng("twister") % For reproducability

N = 1; % Generate an N-member ensemble
errA   = zeros(N,1);
errQ   = zeros(N,1);
errGam = zeros(N,1);
E0 = zeros(N,1);
E1 = zeros(N,1);

err_fun = @(x,y) norm( vec(x-y) ); % Evaluate the absolute error between two objects
rel_err_fun = @(x,y) norm( x-y )/norm(y); % Evaluate the relative error between two objects

for i = 1:N
    str = strcat(num2str(i),"-th ensemble member");
    disp(str)
    
    %%% Generatng synthetic data
    tic
    myopt.Method = "Milstein2";
    [xvec,~] = ST_Generate_timeseries("Colored",A,Q,[],Gam,Tf,Delta_t,myopt); % Generating data
    xvec = xvec(:,1:tt_ratio:end); % Make a sparse observation
    t0 = toc;
    str = strcat("Generating synthetic data took ",num2str(t0,'%3.2f'),"s.");
    disp(str)
    
    %%% Implementation of Colored-LIM
    tic
    % Computing the observed correlation function up to maxlag.
    Corr_Obs = CorrFcn_ST_Obs(xvec,maxlag);

    N0 = Corr_Obs(:,:,1);
    N1 = (Corr_Obs(:,:,2)-Corr_Obs(:,:,1))/dt;
    N2 = (Corr_Obs(:,:,3)-2*Corr_Obs(:,:,2)+Corr_Obs(:,:,1))/dt^2;
    
    E0(i) = err_fun( Corr_Colored_True, Corr_Obs ); % The absolute error due to observation
    str = strcat("Numerical error due to data Generation: E0 = ", num2str(E0(i),'%3.3f'));
    disp(str);
    
    % The Gam-selection algorithm
    Gam1 = (0.6*Gam):0.003:(1.4*Gam);
    % Gam1 = 0.5;
    err = zeros(size(Gam1)); % Record the absolute error bewteen a Colored-LIM correlation and the observed correlation.
    for j = 1:length(Gam1)
        [A_Num,Q_Num] = ST_LIM('Colored',N0,N1,N2,Gam1(j));
        Corr_Colored = ST_CorrFcn('Colored',dt,maxlag,A_Num,Q_Num,[],Gam1(j));
        err(j) =  err_fun( Corr_Colored, Corr_Obs );
    end
    [E1(i),b] = min(err); % Pick up the best Gam that gives the smallest error.
    str = strcat("Numerical error in minimization: E1 = ", num2str(E1(i),'%3.3f'));
    disp(str);
    
    errGam(i) = rel_err_fun(Gam1(b),Gam); % The relative error of Gam

    [A_Num,Q_Num] = ST_LIM('Colored',N0,N1,N2,Gam1(b));
    Corr_Colored = ST_CorrFcn('Colored',dt,maxlag,A_Num,Q_Num,[],Gam1(b));
    errA(i) = rel_err_fun( A_Num, A );
    errQ(i) = rel_err_fun( Q_Num, Q );

    t1 = toc;

    str = strcat("Implementation time of Colored-LIM: ", num2str(t1,'%3.2f'),"s");
    disp(str);
end


%% Show results %%

disp("The statistics of Colored-LIM performance")
str = strcat("e_Gam: ", strcat(num2str(100*median(errGam),'%3.1f'),"%"),"+",strcat(num2str(100*std(errGam,1),'%3.1f'),"%"));
disp(str);
str = strcat("e_A: ", strcat(num2str(100*median(errA),'%3.1f'),"%"),"+",strcat(num2str(100*std(errA,1),'%3.1f'),"%"));
disp(str);
str = strcat("e_Q: ", strcat(num2str(100*median(errQ),'%3.1f'),"%"),"+",strcat(num2str(100*std(errQ,1),'%3.1f'),"%"));
disp(str);
str = strcat("E0: ", strcat(num2str(median(E0),'%3.4f')),"+",strcat(num2str(std(E0,1),'%3.4f')));
disp(str);
str = strcat("E1: ", strcat(num2str(median(E1),'%3.4f')),"+",strcat(num2str(std(E1,1),'%3.4f')));
disp(str);
%str = strcat("Numerical error in minimization: ", strcat(num2str(100*errQ2,'%3.1f'),"%"));
%disp(str);

%%% The demostration of the basic use of Colored-LIM ends here %%% 

%% Making Figure 1 
% Run this only if N = 1!
close all;

fig = figure;
set(gcf, 'Units', 'inch', 'InnerPosition', [1, 1, 3.40, 2.40]);

rho = (0:maxlag)*dt;
i = 0;
for j = 1:3
    for k = 1:3
        i = i + 1;

        subplot(3,3,i)
        y1 = vec(Corr_Colored_True(j,k,:));
        P = plot(rho,y1); hold on;
        P.Color = [0.6350 0.0780 0.1840];
        P.LineWidth = 1.0;
        y2 = vec(Corr_Obs(j,k,:));
        P = plot(rho,y2); hold on;
        P.Color = 0.5020*[1 1 1];
        P.LineWidth = 2.0;
        P.LineStyle = ":";
        % P.Marker = 'o';
        % P.MarkerSize = 5;
        y3 = vec(Corr_Colored(j,k,:));
        P = plot(rho,y3); hold on;
        P.Color = [0 0.4470 0.7410];
        P.LineWidth = 1.0;

        ax = gca;
        ym = min([y1;y2;y3]);
        yM = max([y1;y2;y3]);
        yd = 0.5*(yM-ym)/2;
        ax.YLim = [ym-yd,yM+yd];
        if j == 3 && k == 3
            lgd = legend("$\textbf{\textrm{K}}$","$K^\mathrm{obs}$","$K^\mathrm{ST-C}$");
            lgd.Location = 'southwest';
            lgd.ItemTokenSize = [10,10];
            lgd.Box = 'off';
        end
        AX(i) = ax;
    end
end

x0 = 0.085; y0 = 0.105;
dx = 0.267; Dx = 0.315;
dy = 0.280; Dy = 0.300;
xx = x0:Dx:(x0+2*Dx);
yy = flip( y0:Dy:(y0+2*Dy) );
i = 0;
for j = 1:3
    for k = 1:3
        i = i + 1;
        ax = AX(i);
        ax.Position = [xx(k),yy(j),dx,dy];
        ax.XLim = dt*[-1,maxlag];
        
        ax.XAxis.FontSize = 6;
        if j == 3
            ax.XAxis.TickValues = 0:0.5:1;
            ax.XAxis.TickLabels = arrayfun(@(a)strcat(num2str(a,'%3.1f')),ax.XAxis.TickValues,'uni',0);
        else
            ax.XAxis.TickValues = [ ];
        end
        

        ax.YAxis.FontSize = 6;
        


        if length( ax.YAxis.TickValues ) >= 4
            ax.YAxis.TickValues = ax.YAxis.TickValues(1: ceil( length( ax.YAxis.TickValues )/3 ):end);            
        end
        ax.YAxis.TickLabels = arrayfun(@(a)strcat(num2str(a,'%3.1f')),ax.YAxis.TickValues,'uni',0);


        ax.YTickMode = 'manual';
        ax.YTickLabelMode = 'manual';

    end
end


% New axis
ax = axes('Position', [0 0 1 1], 'Visible', 'off'); 
% Add text using the 'text' function, relative to the whole figure
t1 = text(x0+(2*Dx+dx)/2, 0.000, 'Lag variable $\gamma$', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Parent', ax);
t1.VerticalAlignment = 'bottom';
t1.FontSize = 8;
t1 = text(0, y0+(2*Dy+dy)/2, 'Correlation function', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Parent', ax);
t1.VerticalAlignment = 'top';
t1.Rotation = 90;
t1.FontSize = 8;

lgd.NumColumns = 3;
lgd.Position(1) = 0.02;
lgd.Position(2) = -0.015;

for i = 1:3
    for j = 1:3
        str = strcat( "$i = ", num2str(i), "$, $j = ", num2str(j), "$" );
        t1 = text(xx(i)+0.4*dx, yy(j)+0.8*dy, str, 'Parent', ax);
        % t1.Rotation = 90;
        t1.HorizontalAlignment = 'left';
        t1.VerticalAlignment = 'bottom';
        t1.FontSize = 7;
    end
end

mysavegcf('ST_Colored_Fig1_thesis')

%%
for i = 1:3
    for j = 1:3

        set(fig,'CurrentAxes',ax)
        
        str = strcat( "$i = ", num2str(i), "$, $j = ", num2str(j), "$" );
        t = text( xx(i)+0.4*dx, yy(j)+0.8*dy, str, Units = "normalized");
        t.FontSize = 7;
        t.VerticalAlignment = 'bottom';

        % Example position in normalized units
        norm_pos = [t.Position(1), t.Position(2), 0.1, 0.1];
        
        % Get the current axis limits
        %ax = gca; % Get current axis handle
        xlim = get(ax, 'XLim');
        ylim = get(ax, 'YLim');
        
        % Get the position of the axis in normalized units
        ax_pos = get(ax, 'Position');
        
        % Convert normalized units to data units
        data_x = (norm_pos(1) - ax_pos(1)) / ax_pos(3) * (xlim(2) - xlim(1)) + xlim(1);
        data_y = (norm_pos(2) - ax_pos(2)) / ax_pos(4) * (ylim(2) - ylim(1)) + ylim(1);
        data_w = norm_pos(3) / ax_pos(3) * (xlim(2) - xlim(1));
        data_h = norm_pos(4) / ax_pos(4) * (ylim(2) - ylim(1));
        
        % The result in data units
        data_pos = [data_x, data_y, data_w, data_h];

    end
end

drawnow

mysavegcf('ST_Colored_Fig1_thesis')
