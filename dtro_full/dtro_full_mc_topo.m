clear all
close all

mc_runs=200;

nbnodes=31;
nbsensnode=15*ones(1,nbnodes);
nbsamples=10000;

nbsens=sum(nbsensnode);

Q_data=5;
Q=5;

h=waitbar(0,'Computing');

for n_runs=1:mc_runs

    [Ufull,Vfull]=create_data(nbsensnode,nbsamples,nbnodes,Q_data);

    Ruu=1/nbsamples*conj(Ufull'*Ufull);
    Rvv=1/nbsamples*conj(Vfull'*Vfull);

    Ruu=make_sym(Ruu);
    Rvv=make_sym(Rvv);

    params=struct;
    params.nbsens=nbsens;
    params.Q=Q;
    params.nbnodes=nbnodes;
    params.nbsensnode=nbsensnode;
    params.denom_sum=0;
    params.sgn_sync=1;

    conv=struct;
    conv.tol_rho=1e-12;
    conv.nbiter=-1;

    debug=0;

    % Estimate filter

    data=struct;
    data.R_first=Ruu;
    data.R_second=Rvv;
    [W_star,~,~,~]=trace_ratio(params,data,conv);
    rho_star=trace(W_star'*Ruu*W_star)./trace(W_star'*(Rvv)*W_star);

    % Distributed trace ratio

    conv.tol_rho=-1;
    conv.nbiter=3000;
    params.follow_path=2;

    [rho_track_ER_03,~,norm_star_ER_03,alg_connect_ER_03]=run_tro('ER',params,data,conv,debug,W_star,0.3);
    [rho_track_ER_06,~,norm_star_ER_06,alg_connect_ER_06]=run_tro('ER',params,data,conv,debug,W_star,0.6);
    [rho_track_ER_09,~,norm_star_ER_09,alg_connect_ER_09]=run_tro('ER',params,data,conv,debug,W_star,0.9);
    
    [rho_track_FC,~,norm_star_FC]=run_tro('FC',params,data,conv,debug,W_star);
    [rho_track_TP,~,norm_star_TP]=run_tro('TP',params,data,conv,debug,W_star);
    [rho_track_TS,~,norm_star_TS]=run_tro('TS',params,data,conv,debug,W_star);
    [rho_track_TR2,~,norm_star_TR2]=run_tro('TR',params,data,conv,debug,W_star,2);
    [rho_track_TR5,~,norm_star_TR5]=run_tro('TR',params,data,conv,debug,W_star,5);
    [rho_track_RT,~,norm_star_RT,alg_connect_RT]=run_tro('RT',params,data,conv,debug,W_star);

    norm_star_cell_ER_03{n_runs}=norm_star_ER_03;
    norm_star_cell_ER_06{n_runs}=norm_star_ER_06;
    norm_star_cell_ER_09{n_runs}=norm_star_ER_09;

    norm_star_cell_FC{n_runs}=norm_star_FC;
    norm_star_cell_TP{n_runs}=norm_star_TP;
    norm_star_cell_TS{n_runs}=norm_star_TS;
    norm_star_cell_TR2{n_runs}=norm_star_TR2;
    norm_star_cell_TR5{n_runs}=norm_star_TR5;
    norm_star_cell_RT{n_runs}=norm_star_RT;
    
    alg_c_ER_03{n_runs}=alg_connect_ER_03;
    alg_c_ER_06{n_runs}=alg_connect_ER_06;
    alg_c_ER_09{n_runs}=alg_connect_ER_09;
    
    alg_c_RT{n_runs}=alg_connect_RT;
    
    waitbar(n_runs/mc_runs,h,[sprintf('%3.2f',100*n_runs/mc_runs),'%'])
    
end

save topo.mat norm_star_cell_RT norm_star_cell_ER_03 norm_star_cell_ER_06 ...
    norm_star_cell_ER_09 norm_star_cell_FC norm_star_cell_TP norm_star_cell_TS ...
    norm_star_cell_TR2 norm_star_cell_TR5 alg_c_ER_03 alg_c_ER_06 alg_c_ER_09 alg_c_RT

close(h)

function [rho_track,norm_track,norm_star_track,alg_connect]=run_tro(type,params,data,conv,debug,W_star,reg_param)
    if isequal(type,'FC')
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func(params,data,conv,debug,W_star);

    elseif isequal(type,'RT')
        [tree]=branch();
        tree=sort(tree);
        while length(tree)<params.nbnodes
            [tree]=branch();
            tree=sort(tree);
        end
        tree=tree(1:params.nbnodes);
        adj=comp_adj(tree);
        alg_connect=algebraic_connectivity(adj);
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_tree(params,data,tree,conv,debug,W_star);
    elseif isequal(type,'TP')
        tree=create_tree(params.nbnodes,'path');
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_tree(params,data,tree,conv,debug,W_star);
    elseif isequal(type,'TS')
        tree=create_tree(params.nbnodes,'star');
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_tree(params,data,tree,conv,debug,W_star);
    elseif isequal(type,'TR')
        tree=create_tree(params.nbnodes,'regular',reg_param);
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_tree(params,data,tree,conv,debug,W_star);
    elseif isequal(type,'TH')
        [tree]=branch();
        tree=sort(tree);
        while length(tree)<params.nbnodes
            [tree]=branch();
            tree=sort(tree);
        end
        tree=tree(1:params.nbnodes);
        adj=comp_adj(tree);
        hybrid=make_hybrid(tree,adj);
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_graph(params,data,hybrid,conv,debug,W_star);

    elseif isequal(type,'ER')
        p=struct;
        p.connected=1;
        G=gsp_erdos_renyi(params.nbnodes,reg_param,p);
        er_adj=double(full(G.A));
        er_graph=graph(er_adj);
        bins=conncomp(er_graph);
        alg_connect=algebraic_connectivity(er_adj);
        [~,rho_track,norm_track,norm_star_track]=distrib_trace_func_graph(params,data,er_adj,conv,debug,W_star);

    end
end

function alg_connect=algebraic_connectivity(adj)
    nbnodes=size(adj,1);
    L=sum(adj).*eye(nbnodes)-adj;
    l=eig(L);
    l=sort(l,'ascend');
    alg_connect=l(2);
end

