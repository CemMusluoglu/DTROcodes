clear all
close all


tellernodes=0;
nbrep=1;
nbnodes=30;
nbsensnode=15*ones(1,nbnodes); %number of sensors at each node (one element for each node)
nbsamples=10000; %number of observations per iteration of the DACMEE algorithm

nbsens=sum(nbsensnode); %total number of sensors

Q_data=5;

[Ufull,Vfull]=create_data(nbsensnode,nbsamples,nbnodes,Q_data);

Ruu=1/nbsamples*conj(Ufull'*Ufull);
Rvv=1/nbsamples*conj(Vfull'*Vfull);

Ruu=make_sym(Ruu);
Rvv=make_sym(Rvv);

Q=5;

params=struct;
params.nbsens=nbsens;
params.Q=Q;
params.nbnodes=nbnodes;
params.nbsensnode=nbsensnode;
params.nbsamples=nbsamples;
params.denom_sum=0;

conv=struct;
conv.tol_rho=-1;
conv.nbiter=30*30;

debug=0;

%% Estimate filter

conv.tol_rho=1e-12;
conv.nbiter=-1;

data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
[W1_star,~,~,~]=trace_ratio(params,data,conv);


%% Distributed trace ratio

conv.tol_rho=-1;
conv.nbiter=3000;

data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
[w_1,rho_track_1,norm_1_track,norm1_star_track]=distrib_trace_func(params,data,conv,debug,W1_star);

%% Distributed trace ratio tree

[tree]=branch();
tree=sort(tree);
while length(tree)<nbnodes
    [tree]=branch();
    tree=sort(tree);
end
tree=tree(1:nbnodes);
trimtreeplot(tree)
%%
leafs=setdiff([1:nbnodes],tree);
path_tree=comp_path(tree,leafs);

%%
conv.tol_rho=-1;
conv.nbiter=30000;
data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
debug=0;
%debug.follow_path=1;
%[w_1_tree,rho_track_1_tree,norm1_track_tree,norm1_star_track_tree]=distrib_trace_func_tree(params,data,tree,conv,debug,W1_star);

params.follow_path=0;
[w_1_tree_nop,rho_track_1_tree_nop,norm1_track_tree_nop,norm1_star_track_tree_nop]=distrib_trace_func_tree(params,data,tree,conv,debug,W1_star);

%% Hybrid Graphs

adj=comp_adj(tree);
hybrid=make_hybrid(tree,adj);
%%
data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
[w_1_h,rho_track_1_h,norm1_track_h,norm1_star_track_h]=distrib_trace_func_graph(params,data,hybrid,conv,debug,W1_star);

%% Erdos Renyi graphs
p=struct;
p.connected=1;
G=gsp_erdos_renyi(nbnodes,0.6,p);
er_adj=double(full(G.A));
er_graph=graph(er_adj);
bins=conncomp(er_graph);
%%
data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
[w_1_er,rho_track_1_er,norm1_track_er,norm1_star_track_er]=distrib_trace_func_graph(params,data,er_adj,conv,debug,W1_star);


%% Trace ratio

data=struct;
data.R_first=Ruu;
data.R_second=Rvv;
[W_1,rho_track_1_full,norm1_track_full,norm1_star_track_full]=trace_ratio(params,data,conv,W1_star);


%%
rho_star=trace(W1_star'*Ruu*W1_star)/trace(W1_star'*(Rvv)*W1_star);
loglog(rho_star-rho_track_1,'r','LineWidth',2)
hold on
loglog(rho_star-rho_track_1_tree,'g','LineWidth',2)
%semilogx(rho_track_1_full,'k','LineWidth',2)
legend('FC','TG')
xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$\rho^*-\rho^i$','Interpreter','latex','FontSize', 20)
grid on

figure
x_int=[2:conv.nbiter];
loglog(x_int,diff(rho_track_1),'r','LineWidth',2)
hold on
loglog(x_int,diff(rho_track_1_tree),'g','LineWidth',2)
%semilogy(x_int,abs(diff(rho_track_1_full)),'k','LineWidth',2)
legend('FC','TG')
xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$\rho^i-\rho^{i-1}$','Interpreter','latex','FontSize', 20)

figure
plot(x_int,norm_1_track,'r','LineWidth',2)
hold on
plot(x_int,norm1_track_tree,'g','LineWidth',2)
%plot(x_int,norm1_track_full,'k','LineWidth',2)
legend('FC','TG')
xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$C||X^i-X^{i-1}||_F$','Interpreter','latex','FontSize', 20)

figure
x_int=[1:conv.nbiter];
loglog(x_int,norm1_star_track,'r','LineWidth',2)
hold on
loglog(x_int,norm1_star_track_tree,'g','LineWidth',2)
%plot(x_int,norm1_star_track_full,'k','LineWidth',2)
legend('FC','TG')
xlabel('Iteration $i$','Interpreter','latex','FontSize', 20)
ylabel('$C||X^i-X^*||^2_F$','Interpreter','latex','FontSize', 20)
grid on
%{
%% Trace ratio

W1 = randn(size(Ruu,1),Q);
W2 = randn(size(Rvv,1),Q);

[ W , ~,~,~,~,~] = traincsp( X , [1,2] , 2*Q , 'traceratio' ,W1,W2);

W1=W(:,1:Q);
W2=W(:,Q+1:end);

%% GEVD

[W1_eig,D1_eig]=eig(Ruu,Ruu+Rvv);
[~,ind]=sort(diag(D1_eig),'descend');
D1_eig=diag(D1_eig);
D1_eig=diag(D1_eig(ind));
W1_eig=W1_eig(:,ind(1:Q));

[W2_eig,D2_eig]=eig(Rvv,Rvv+Ruu);
[~,ind]=sort(diag(D2_eig),'descend');
D2_eig=diag(D2_eig);
D2_eig=diag(D2_eig(ind));
W2_eig=W2_eig(:,ind(1:Q));

%% Compare maxima

rho1_max_eig=trace(W1_eig'*Ruu*W1_eig)/trace(W1_eig'*(Ruu+Rvv)*W1_eig); %=sum(lambdaGEVD_1^Q)/Q
rho1_max_tr=trace(W(:,1:Q)'*Ruu*W(:,1:Q))/trace(W(:,1:Q)'*(Ruu+Rvv)*W(:,1:Q));
rho1_max_tr>=rho1_max_eig
[test1_eigv,test1_eigvl]=eig(Ruu-rho1_max_eig*(Ruu+Rvv));
[~,ind]=sort(diag(test1_eigvl),'descend');
test1_eigvl=diag(test1_eigvl);
test1_eigvl=diag(test1_eigvl(ind(1:Q)));
test1_eigv=test1_eigv(:,ind(1:Q));
rho1_test=trace(test1_eigv'*Ruu*test1_eigv)/trace(test1_eigv'*(Ruu+Rvv)*test1_eigv);

%% Distributed trace ratio solving GEVD

[V_ruu,D_ruu,~]=svd(Ruu+Rvv);
[V_rvv,D_rvv,~]=svd(Rvv+Ruu);
inv_D_ruu=diag(1./sqrt(diag(D_ruu)));
inv_D_rvv=diag(1./sqrt(diag(D_rvv)));

data=struct;
data.R_first=inv_D_rvv*V_rvv'*Ruu*V_rvv*inv_D_rvv;
data.R_first=(data.R_first+data.R_first')/2;
%data.R_second=inv_D_rvv*V_rvv'*Rvv*V_rvv*inv_D_rvv;
data.R_second=eye(size(Rvv));
[w_1_eig,rho_1_eig,rho_track_1_eig]=distrib_trace_func(params,data,conv,debug);

data=struct;
data.R_first=inv_D_ruu*V_ruu'*Rvv*V_ruu*inv_D_ruu;
data.R_first=(data.R_first+data.R_first')/2;
%data.R_second=inv_D_ruu*V_ruu'*Ruu*V_ruu*inv_D_ruu;
data.R_second=eye(size(Ruu));
[w_2_eig,rho_2_eig,rho_track_2_eig]=distrib_trace_func(params,data,conv,debug);

w_1_eig2=V_rvv*inv_D_rvv*w_1_eig; % =W1_eig=GEVC(Ruu,Ruu+Rvv) !!
w_2_eig2=V_ruu*inv_D_ruu*w_2_eig; % =W2_eig=GEVC(Rvv,Rvv+Ruu) !!


%%
%X=cat(3,Ufull',Vfull');
Y = tmprod(X,w_tot',1);
F = log(var(Y,[],2));
F = squeeze(F);
            
ldamodel = fitcdiscr(F',labels,'Gamma',0.05);
%ldamodel = fitcdiscr(F',[1,-1]);
labest = predict(ldamodel,F')';
results = mean( labels == labest );

%%

%}