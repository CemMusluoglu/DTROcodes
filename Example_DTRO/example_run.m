% Example script to run the TI-DTRO algorithm.

% Author: Cem Musluoglu, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics
% Correspondence: cemates.musluoglu@esat.kuleuven.be

clear all
close all

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
% Number of Monte-Carlo runs
mc_runs=10;

data_cell=cell(mc_runs,1);

for n_runs=1:mc_runs
    
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
    debug=0;

    % Centralized trace ratio, used as ground truth
    [X_star,rho]=trace_ratio(Q,nbsens,Ruu,Rvv,params.denom_sum);

    % Convergence criteria for the distributed trace ratio
    conv=struct;
    conv.tol_rho=-1;
    conv.nbiter=300;

    % Create an adjacency matrix
    % Example: Fully connected
    graph_adj=(ones(nbnodes,nbnodes)-eye(nbnodes));

    % TI-DTRO
    [X,rho_track,norm_track,norm_star_track]=ti_dtro(params,data,graph_adj,conv,debug,X_star);
    
    data_cell{n_runs}=norm_star_track;

end

% Plot the MSE

x_int=[1:conv.nbiter];
q_5=quantile(cell2mat(data_cell),0.5);
q_25=quantile(cell2mat(data_cell),0.25);
q_75=quantile(cell2mat(data_cell),0.75);
loglog(q_5,'b','LineWidth',2);

hold on
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_75)],'b','FaceAlpha','0.2','LineStyle','none')
fill([x_int,fliplr(x_int)],[q_5,fliplr(q_25)],'b','FaceAlpha','0.2','LineStyle','none')
xlim([1,inf])
ylim([1e-10,inf])

xlabel('Iterations','Interpreter','latex')
ylabel('$\frac{1}{MQ}||X^{i}-X^*||_2^F$','Interpreter','latex')
grid on


