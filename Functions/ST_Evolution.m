function [yt,yd,yc] = ST_Evolution(Type,C,A,Q,D,Gam,x0,tp)

[n,N] = size(x0);
yd = zeros(n,N,length(tp)); % Dynamical
yc = zeros(n,N,length(tp)); % Coupling
    

if Type == "Colored" || Type == "CW"

    [~,Gd,Gs] = ST_Evolution_Operator(Type,C,A,Q,D,Gam,tp);
    
    yd = Gd*x0;
    yc = Gs*x0;
    yt = yd + yc; % Total

end