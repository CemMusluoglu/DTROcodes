function [w,rho_track,norm_track,norm_star_track]=ti_dtro(params,data,graph_adj,conv,debug,W_star)

% Function running the distributed trace ratio algorithm on a random
% connected graph.
% INPUTS :
% params : Structure containing the following fields:
%         Q : Number of filters to use (dimension of projected space)
%         nbnodes : Number of nodes in the network
%         nbsensnode (nbnodes x 1): Vector containing the number of sensors for each
%         node
%         nbsens : Sum of the number of sensors for each node (dimension of
%         the network-wide signals). Is equal to sum(nbsensnode)
%         denom_sum : Binary value for determining the trace ratio
%         denominator. If equals 1, the denominator is the sum of the
%         covariance matrices. If 0 not
% data   : Structure containing the following fields:
%         R_first (nbsens x nbsens): Covariance matrix in the numerator
%         R_second (nbsens x nbsens): Covariance matrix in the denominator
% conv   : Parameters concerning the stopping criteria of the algorithm
%         tol_rho : Tolerance in objective: |rho^(i+1)-rho^(i)|>tol_rho
%         nbiter : Max. nb. of iterations.
% If both values are valid, the algorithm will continue until the stricter
% condition (OR). One of the criteria can be chosen explicitly by
% initializing the other to a negative value.
%
% graph_adj (nbnodes x nbnodes): Adjacency (binary) matrix, with 
%           graph_adj(i,j)=1 if i and j are connected. Otherwise 0. 
%           graph_adj(i,i)=0.
% debug  : If debug equals 1, dynamically plot first projection vector
%          across iterations
% W_star (nbsens x Q): (Optional) True projection matrix, computed 
%          for example with the centralized algorithm. Allows to compare 
%          convergence, if it is not provided, there might be a difference in 
%          the signs of the columns of the output of this algorithm and W_star
%
% OUTPUTS : 
% (We denote by nbiter the total number of iterations, even though the number
% might be larger if the other stopping criterion is also used)
% w (nbsens x Q)              : Projection matrix
% rho_track (nbiter x 1)      : Sequence of objective values across iterations
% norm_track (nbiter x 1)     : Sequence of ||W^(i+1)-W^(i)||_F^2
% norm_star_track (nbiter x 1): Sequence of ||W^(i)-W^*||_F^2
%
% EXTERNAL FUNCTIONS
% dijkstra.m : Joseph Kirk (2020). Dijkstra's Shortest Path Algorithm 
% https://www.mathworks.com/matlabcentral/fileexchange/20025-dijkstra-s-minimum-cost-path-algorithm?s_tid=prof_contriblnk
% MATLAB Central File Exchange.
%
% Can be replaced by other algorithms pruning the graph to a tree.



nbsens=params.nbsens;
Q=params.Q;
nbnodes=params.nbnodes;
nbsensnode=params.nbsensnode;
denom_sum=params.denom_sum;

Ruu=data.R_first;
Rvv=data.R_second;

tol_rho=conv.tol_rho;
nbiter=conv.nbiter;

winit=randn(nbsens,Q);
winit=normc(winit);
rho=trace(winit'*Ruu*winit)/trace(winit'*denom_mat(Ruu,Rvv,denom_sum)*winit);

teller=1;
for k=1:nbnodes
    W{k}=winit(teller:teller+nbsensnode(k)-1,:);
    w=winit;
    teller=teller+nbsensnode(k);
end
clear teller

i=0;
rho_old=rho+1;

rho_track=[];
rho_track_test=[];
norm_track=[];
norm_star_track=[];

adj=graph_adj;

while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    w_old=w;
    
    q=rem(i,nbnodes)+1;
    
    [neighbors,path]=find_path(q,adj);
    
    Nu=constr_Nu(neighbors,path);
    
    C_q=constr_C(W,Q,q,nbsensnode,nbnodes,neighbors,Nu);

    [Ruucompressed,Rvvcompressed,Compressor]=compress(Ruu,Rvv,C_q);

    WGnew=comp_WG(Ruucompressed,rho,Rvvcompressed,Compressor,Q,denom_sum);

    rho_old=rho;
    
    if(debug==1)
        clear B_int
        clear s_int
        clear ind_int
        [B_int,s_int]=eig(Ruu-rho*denom_mat(Ruu,Rvv,denom_sum));
        [~,ind_int]=sort(diag(s_int),'descend');
        wopt_int=B_int(:,ind_int(1:Q));
        wopt_int=normc(real(wopt_int));
    end
    
    rho=trace(WGnew'*Ruucompressed*WGnew)/trace(WGnew'*denom_mat(Ruucompressed,Rvvcompressed,denom_sum)*WGnew);
    rho_track=[rho_track,rho];
    
    W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes,neighbors,Nu);
    
    w=form_mat(W,w,Q,nbsensnode,nbnodes);

    i=i+1;
    
    if(debug==1)
        for l=1:Q
            if sum(sum((wopt_int(:,l)-w(:,l)).^2))>sum(sum((-wopt_int(:,l)-w(:,l)).^2))
                w(:,l)=-w(:,l);
            end
        end
    
        plot(wopt_int(:,1),'r')
        hold on
        plot(w(:,1),'b')
        ylim([1.2*min(real(wopt_int(:,1))) 1.2*max(real(wopt_int(:,1)))]);
        hold off
        drawnow
    end
    
    if i>1
        norm_track=[norm_track,norm(w-w_old,'fro').^2/numel(w)];
    end
    
    if nargin>5

        for l=1:Q
            if sum(sum((W_star(:,l)-w(:,l)).^2))>sum(sum((-W_star(:,l)-w(:,l)).^2))
                w(:,l)=-w(:,l);
            end
        end

        norm_star_track=[norm_star_track,norm(w-W_star,'fro')^2/numel(W_star)];
    end

end

end

function [neighbors,path]=find_path(q,adj)
    [dist,path]=shortest_path(adj,q);
    neighbors=find(cell2mat(cellfun(@(c) length(c), path, 'uniform', false))==2);
    neighbors=sort(neighbors);
end

function Nu=constr_Nu(neighbors,path)

    nb_neighbors=length(neighbors);
    Nu=cell(nb_neighbors,1);
    for k=1:nb_neighbors
        Nu{k}=find(cell2mat(cellfun(@(c) ismember(neighbors(k),c), path, 'uniform', false))==1);
    end
    
end


function C_q=constr_C(W,Q,q,nbsensnode,nbnodes,neighbors,Nu)
   
    nb_neighbors=length(neighbors);

    ind=0:nb_neighbors-1;
    
    C_q=zeros(sum(nbsensnode),nbsensnode(q)+nb_neighbors*Q);
    C_q(:,1:nbsensnode(q))=[zeros(sum(nbsensnode(1:q-1)),nbsensnode(q));...
        eye(nbsensnode(q)); zeros(sum(nbsensnode(q+1:nbnodes)),nbsensnode(q))];
    for k=1:nb_neighbors
        ind_k=ind(k);
        for n=1:length(Nu{k})
            Nu_k=Nu{k};
            l=Nu_k(n);
            C_q(sum(nbsensnode(1:l-1))+1:sum(nbsensnode(1:l)),...
                nbsensnode(q)+ind_k*Q+1:nbsensnode(q)+ind_k*Q+Q)=W{l};
        end
    end
    
end

function [Ruucompressed,Rvvcompressed,Compressor]=compress(Ruu,Rvv,C_q)

    Ruucompressed=C_q'*Ruu*C_q;
    Rvvcompressed=C_q'*Rvv*C_q;
    Ruucompressed=make_sym(Ruucompressed);
    Rvvcompressed=make_sym(Rvvcompressed);

    Compressor=C_q'*C_q;
    Compressor=make_sym(Compressor);

end

function WG=comp_WG(Ruucompressed,rho,Rvvcompressed,Compressor,Q,denom_sum)
    [B,s]=eig(Ruucompressed-rho*denom_mat(Ruucompressed,Rvvcompressed,denom_sum),Compressor);
    [~,ind]=sort(diag(s),'descend');
    WG=B(:,ind(1:Q));
end

function W=update_W_efficient(W,WGnew,Q,q,nbsensnode,nbnodes,neighbors,Nu)
 
    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end
    
    nb_neighbors=length(neighbors);
    ind=0:nb_neighbors-1;
    
    for l=1:q-1
        for k=1:nb_neighbors
            if ~isempty(find(Nu{k} == l))
                start_r=nbsensnode(q)+ind(k)*Q+1;
                stop_r=nbsensnode(q)+ind(k)*Q+Q;
            end
        end
        G=diag(sign(diag(WGnew(start_r:stop_r,:))));
        W{l}=W{l}*G;
    end
    for l=q+1:nbnodes
        for k=1:nb_neighbors
            if ~isempty(find(Nu{k} == l))
                start_r=nbsensnode(q)+ind(k)*Q+1;
                stop_r=nbsensnode(q)+ind(k)*Q+Q;
            end
        end
        G=diag(sign(diag(WGnew(start_r:stop_r,:))));
        W{l}=W{l}*G;
    end
    
    
end

function W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes,neighbors,Nu)

    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end
    
    nb_neighbors=length(neighbors);
    ind=0:nb_neighbors-1;
    
    for l=1:q-1
        for k=1:nb_neighbors
            if ~isempty(find(Nu{k} == l))
                start_r=nbsensnode(q)+ind(k)*Q+1;
                stop_r=nbsensnode(q)+ind(k)*Q+Q;
            end
        end
        W{l}=W{l}*WGnew(start_r:stop_r,:);
    end
    for l=q+1:nbnodes
        for k=1:nb_neighbors
            if ~isempty(find(Nu{k} == l))
                start_r=nbsensnode(q)+ind(k)*Q+1;
                stop_r=nbsensnode(q)+ind(k)*Q+Q;
            end
        end
        W{l}=W{l}*WGnew(start_r:stop_r,:);
    end

end

function w=form_mat(W,w,Q,nbsensnode,nbnodes)

    teller=1;
    for l=1:nbnodes
        w(teller:teller+nbsensnode(l)-1,:)=W{l}(:,1:Q);
        teller=teller+nbsensnode(l);
    end

end

function R_denom=denom_mat(Ruu,Rvv,denom_sum)

    if denom_sum==1
        R_denom=Ruu+Rvv;
    else
        R_denom=Rvv;
    end
end

