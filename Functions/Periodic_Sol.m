function [xp,h] = Periodic_Sol(A,b,tp)


    % Find the periodic solution xp of the linear equation dX = A(t)Xdt + b(t)dt.
    % To ensure uniqueness, dX = A(t)Xdt is assumed to have no nonzero periodic solution.
    % The period is assumed to be 1.

    n = length(b(1)); % Dimension
    N = length(tp); % Number of sampling

    tspan = (0:N)/N;

    [~,y] = ode45(@(s,y) myfun(A,s,y,n), tspan, eye(n));
    y = reshape(y,length(tspan),n,[]);
    Psi = permute(y,[2 3 1]);
    PsiT = Psi(:,:,end); % Monodromy matrix


    integrand = zeros(n,length(tspan));
    for j = 1:N+1
        integrand(:,j) = Psi(:,:,j)\b(tspan(j));
    end

    int_integrand = trapz(tspan,integrand,2);

    x0 = ( eye(n) - PsiT ) \ ( PsiT * int_integrand ); % Initial condition at t = 0

    h.x0 = x0;

    if norm(imag(x0)) >0 % Numerical instability
        % ( eye(n) - PsiT )
        % ( PsiT * int_integrand )
        % ( eye(n) - PsiT ) \ ( PsiT * int_integrand )
        disp('Numerical instability may occur!')
        x0 = zeros(n,1);
        pause(100)
    end

    myMaxStep = 1e-2;
    options = odeset(MaxStep=myMaxStep);
    [~,xp] = ode45(@(s,y) A(s)*y+b(s), [-1,tspan], x0, options);
    xp = xp';
    xp = xp(:,2:end);



    % Some criteria to determine whether the resulting periodic function is valid.
    k = 0; myflag = 0;
    notperiodic = @(xx) norm( xx(:,1) - xx(:,end) ) / norm( (xx(:,1)+xx(:,end))/2 ) > 0.005; 
    not_periodic_but_well = @(xx) norm( xx(:,1) - xx(:,end) ) / norm( (xx(:,1)+xx(:,end))/2 ) > 0.02; % 1% for CS-Colored-LIM paper.
    not_periodic_terrible = @(xx) norm( (xx(:,1)+xx(:,end))/2 )/n > 10;

    while notperiodic( xp ) && myflag == 0 

        k = k + 1;
        
        x0 = xp(:,end);
        [t,xp] = ode45(@(s,y) A(s)*y+b(s), tspan, x0, options);
        xp = xp';

        if not_periodic_terrible( xp )
            myflag = 1; % Probably no attractor!
        end

        if mymod(k,10) == 10 % If not converge over time, decrease maxstep.
            myMaxStep = 0.5*myMaxStep;
            options = odeset(MaxStep=myMaxStep);
        end

        if k > 32

            if not_periodic_but_well( xp )
                xp
                disp('Periodic sol. cannot be found!')

                close all;
                figure
                plot(xp(1,:))
                pause(1)
                
            else 
                myflag = 1;
            end

        end
    end

    if myflag == 1 % Does not find the periodic solution
        xp = nan;
    else % Accept the periodic solution
        xp = xp(:,1:end-1);
    end

    % disp('-------------!!!!')
    h.flag = myflag;


    
end

function y = myfun(A,t,x,n)

    x = reshape(x,n,[]);
    y = A(t) * x;
    y = y(:);
end
