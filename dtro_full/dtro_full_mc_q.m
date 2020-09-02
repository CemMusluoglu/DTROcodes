clear all
close all

mc_runs=200;

nbnodes=30;
nbsensnode=15*ones(1,nbnodes);
nbsamples=10000;

nbsens=sum(nbsensnode);

D=10;

for Q=[1,3,5,7]
    
    Q_data=Q;

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
        conv.nbiter=30*100;
        params.follow_path=2;

        [rho_track_FC,~,norm_star_FC]=run_tro('FC',params,data,conv,debug,W_star);
        [rho_track_RT,~,norm_star_RT,alg_connect_RT]=run_tro('RT',params,data,conv,debug,W_star);

        %rho_star_cell{n_runs}=rho_star;

        norm_star_cell_FC{n_runs}=norm_star_FC;
        norm_star_cell_RT{n_runs}=norm_star_RT;



        %rho_cell_RT{n_runs}=rho_track_RT;
        %rho_cell_FC{n_runs}=rho_track_FC;

        n_runs

    end
    
    fname=strcat('Q',num2str(Q));
    save(fname,'norm_star_cell_FC','norm_star_cell_RT');

    
end

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

