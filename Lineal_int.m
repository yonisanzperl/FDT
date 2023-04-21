function [FC,CV,A]=Lineal_int(gC,sigma)
N=size(gC,1);
a=-1;  %% -0.4 for awake sleep  %% -1 for HCP

% Jacobian:
A = a*eye(N) + gC;
Qn = (sigma^2)*eye(N);

CV=sylvester(A,A',-Qn);
CV=(CV+CV')/2;
FC=corrcov(CV);

