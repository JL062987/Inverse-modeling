function [B,A,C,Q,h] = CS_nLIM(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,Ktheta,opt)

    [n,~] = size(K0);

    %% Set up

    % opt
    % ├── Sparse_Identification
    % │   ├── Type
    % │   │   └── "Empirical" (default) or "Prescribed"
    % │   ├── Solve_Lambda
    % │   │   └── "Grid_Search" (default) or a function handle
    % │   ├── user_opt
    % │   │   └── Optional parameters passed to Solve_Lambda
    % │   └── FISTA_opt
    % │       └── Optional parameters passed to fista_general
    % │
    % ├── Fourier_Mode
    % │   ├── Type
    % │   │   └── "Fixed_phase" (default) or "Phase_average"
    % │   └── FM
    % │       └── Fourier truncation mode; required only when
    % │           Type = "Phase_average" (default: 2)
    % │
    % ├── Solve_Diffusion
    % │   └── "Analytic" (default), "Residual", or "Loss_Function"
    % │
    % └── xvec
    %     └── State data required only when Solve_Diffusion = "Residual"

    %%% Sparse_Identification
    if ~isfield(opt,"Sparse_Identification")
        opt.Sparse_Identification = struct;
    end

    if ~isfield(opt.Sparse_Identification,"Type")
        opt.Sparse_Identification.Type = "Empirical"; % Empirical/Prescribed
    end

    if ~isfield(opt.Sparse_Identification,"Solve_Lambda")
        opt.Sparse_Identification.Solve_Lambda = "Grid_Search"; % Should be "Grid_Search" or a function handle
        opt.Sparse_Identification.user_opts = struct();
    else
        if ~isfield(opt.Sparse_Identification,"user_opt")
            opt.Sparse_Identification.user_opts = struct();
        end
    end

    if ~isfield(opt.Sparse_Identification,"FISTA_opts")
        opt.Sparse_Identification.FISTA_opts = struct();
        opt.Sparse_Identification.FISTA_opts.Display_Details = true;
    end

    %%% Fourier mode truncation
    if ~isfield(opt,"Fourier_Mode")
        opt.Fourier_Mode = struct;
    end

    if ~isfield(opt.Fourier_Mode,"Type")
        opt.Fourier_Mode.Type = "Fixed_phase"; % Fixed_phase/Phase_average
    end

    if opt.Fourier_Mode.Type == "Phase_average"
        if ~isfield(opt.Fourier_Mode,"FM")
            opt.Fourier_Mode.FM = 2;
            disp("Fourier mode truncation without specifying FM! Default: FM = 2.")
        end
    else 
        opt.Fourier_Mode.FM = [];
    end

    %%% Diffusion
    if ~isfield(opt,"Solve_Diffusion") % Analytic/Residual/Loss_Function
        opt.Solve_Diffusion = "Analytic";
    end

    if opt.Solve_Diffusion == "Residual"
        if isfield(opt,"xvec")
            xvec = opt.xvec;
        else
            error('Missing input: opt.xvec is required when Solve_Diffusion is "Residual".')
        end
    end

    if Type == "White"; Gam = 0; end

    %% Pass opt to functions
    myopts.Method = opt.Fourier_Mode.Type;
    myopts.FM = opt.Fourier_Mode.FM;
    myopts.Solve_Lambda = opt.Sparse_Identification.Solve_Lambda;
    myopts.user_opts = opt.Sparse_Identification.user_opts;
    myopts.FISTA_opts = opt.Sparse_Identification.FISTA_opts;

    %% Sparse identification

    if opt.Sparse_Identification.Type == "Empirical"
        [idx_nonzeroB,h_SI] = CS_nLIM_Sparse_Identification(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,myopts); % Data-driven; if not successfully identifying the correct components, consider making longer simulations or consider a larger range in "opts_in.Solve_Lambda.Lambda_range = 0:0.1:0.5".
        for j = 1:size(h_SI.ActiveComponents,1) % Mask this forloop if active nonlinear components are prescribed!
            ActiveComp = h_SI.ActiveComponents(j,:);
            str = strcat("Active components: (i,j,k) = (",num2str(ActiveComp(1)),",",num2str(ActiveComp(2)),",",num2str(ActiveComp(3)),")");
            disp(str)
        end
        idx_nonzeroB = idx_nonzeroB(:); % Use active nonlinear components identified by Lasso regression
        idx_nonzeroA = ones(n^2,1);
        idx_nonzeroC = ones(n^1,1); 

    else % Prescribed

        [~,JJ,KK] = ndgrid(1:n,1:n,1:n);
        if isfield(opt.Sparse_Identification,"idx_nonzeroB")
            idx_nonzeroB = opt.Sparse_Identification.idx_nonzeroB(:); 
        else
            idx_nonzeroB = zeros(n^3,1);
            idx_nonzeroB(vec(JJ<=KK))=1;
        end
        if isfield(opt.Sparse_Identification,"idx_nonzeroA"); idx_nonzeroA = opt.Sparse_Identification.idx_nonzeroA(:); else; idx_nonzeroA = ones(n^2,1) ; end
        if isfield(opt.Sparse_Identification,"idx_nonzeroC"); idx_nonzeroC = opt.Sparse_Identification.idx_nonzeroC(:); else; idx_nonzeroC = ones(n^1,1) ; end
        
    end

    %% Solve dynamics
    [~,B,A,C,~] = CS_nLIM_Dynamics(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzeroB,idx_nonzeroA,idx_nonzeroC,myopts);
    

    %% Solve diffusion
    

    switch opt.Solve_Diffusion
        case "Analytic"
            Q = CS_nLIM_Diffusion(Type,B,A,C,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,Ktheta);
        case "Residual"
            Q = CS_nLIM_Diffusion_Residual_Method(Type,xvec,B,A,C,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta);
        case "Loss_Function"
            % Define you own loss function
       
    end


    %% 
    h = struct();
    

end