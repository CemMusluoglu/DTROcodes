clear all
close all

% Example code to run the TI-DTRO algorithm.

% Number of nodes
nbnodes=30;
% Number of channels per node
nbsensnode=15*ones(1,nbnodes);
% Number of time samples per data block
nbsamples=10000;

% Total number of sensors
nbsens=sum(nbsensnode); 

% Number of desired point sources
Q_data=5;
% Projection dimension (number of filters)
Q=5;
% Total number of point sources
D=10;

% Desired and interfering signals respectively
[Ufull,Vfull]=create_data(nbsamples,nbsensnode,nbnodes,Q_data,D);

% Covariance matrices
Ruu=1/nbsamples*conj(Ufull'*Ufull);
Rvv=1/nbsamples*conj(Vfull'*Vfull);
Ruu=make_sym(Ruu);
Rvv=make_sym(Rvv);
data=struct;
data.R_first=Ruu;
data.R_second=Rvv;

params=struct;
params.nbsens=nbsens;
params.Q=Q;
params.nbnodes=nbnodes;
params.nbsensnode=nbsensnode;
params.denom_sum=0;

% Convergence criteria for the centralized trace ratio
% Stops at whichever is achieved last
% If one of them is non-positive, the positive one will be taken into account
conv=struct;
conv.tol_rho=1e-12;
conv.nbiter=-1;

% If debug==1, there is a plot showing dynamically the comparison
% between the first column of a theoretical optimal vector
% and the estimation.
debug=1;

% Centralized trace ratio, used as ground truth
[W_star,rho]=trace_ratio(Q,nbsens,Ruu,Rvv);

% Convergence criteria for the distributed trace ratio
conv=struct;
conv.tol_rho=-1;
conv.nbiter=300;

% Create an adjacency matrix
% Example: Fully connected
graph_adj=(ones(nbnodes,nbnodes)-eye(nbnodes));

% TI-DTRO
[W,rho_track,norm_track,norm_star_track]=ti_dtro(params,data,graph_adj,conv,debug,W_star);

% Plot the MSE
loglog(norm_star_track,'b','LineWidth',2)
xlabel('Iterations','Interpreter','latex')
ylabel('$\frac{1}{MQ}||W^{i}-W^*||_2^F$','Interpreter','latex')
grid on



