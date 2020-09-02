function [w,rho_track,norm_track,norm_star_track]=distrib_trace_func_tree(params,data,tree,conv,debug,W_star)

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

leafs=setdiff([1:nbnodes],tree);
path_tree=comp_path(tree,leafs);

path=1:nbnodes;
rand_path=path(randperm(length(path)));

while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    w_old=w;
    
    if params.follow_path==0
        q=rem(i,nbnodes)+1;
        %q=1;
    elseif params.follow_path==1
        q=path_tree(rem(i,length(path_tree))+1);
    else
        q=rand_path(rem(i,length(path_tree))+1);
    end
    
    [Nu,descendents_q]=constr_Nu(q,tree,nbnodes);
    
    C_q=constr_C(W,Q,q,nbsensnode,nbnodes,tree,Nu,descendents_q);
    
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
    
    W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes,tree,Nu,descendents_q);
    
    w=form_mat(W,w,Q,nbsensnode,nbnodes);
    
    %rho_test=trace(w'*Ruu*w)/trace(w'*(Ruu+Rvv)*w);
    %rho_track_test=[rho_track_test,rho_test];

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

function [Nu,descendents_q]=constr_Nu(q,tree,nbnodes)

    children_q=[find(tree==q)];
    nb_children=length(children_q);
    nb_neighbors=nb_children+1;
    ind=1:nb_children;
    Nu=cell(nb_children,1);
    descendents_q=[];
    for l=q+1:nbnodes
        parent=l;
        while ((~ismember(parent,children_q)) && parent~=0)
            parent=tree(parent);
        end
        if parent~=0
            Nu{ind(children_q==parent)}=[Nu{ind(children_q==parent)},l];
            descendents_q=[descendents_q,l];
        end
    end
    
end


function C_q=constr_C(W,Q,q,nbsensnode,nbnodes,tree,Nu,descendents_q)
   
    children_q=[find(tree==q)];
    nb_children=length(children_q);
    nb_neighbors=nb_children+1;
    if q==1
       nb_neighbors=nb_neighbors-1; 
    end
    ind=1:nb_children;
    
    all_nodes=1:nbnodes;
    nodes_left=setdiff(all_nodes,[q,descendents_q]);
    
    C_q=zeros(sum(nbsensnode),nbsensnode(q)+nb_neighbors*Q);
    C_q(:,1:nbsensnode(q))=[zeros(sum(nbsensnode(1:q-1)),nbsensnode(q));...
        eye(nbsensnode(q)); zeros(sum(nbsensnode(q+1:nbnodes)),nbsensnode(q))];
    for k=1:nb_children
        ind_k=ind(k);
        if isempty(nodes_left)
            ind_k=ind_k-1;
        end
        for n=1:length(Nu{k})
            Nu_k=Nu{k};
            l=Nu_k(n);
            C_q(sum(nbsensnode(1:l-1))+1:sum(nbsensnode(1:l)),...
                nbsensnode(q)+ind_k*Q+1:nbsensnode(q)+ind_k*Q+Q)=W{l};
        end
    end
    
    for n=1:length(nodes_left)
        l=nodes_left(n);
        C_q(sum(nbsensnode(1:l-1))+1:sum(nbsensnode(1:l)),...
                nbsensnode(q)+1:nbsensnode(q)+Q)=W{l};
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

function W=update_W_efficient(W,WGnew,Q,q,nbsensnode,nbnodes,tree,Nu,descendents_q)
    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end

    children_q=[find(tree==q)];
    nb_children=length(children_q);
    ind=1:nb_children;
    
    for l=1:q-1
        if ~ismember(l,descendents_q)
            start_r=nbsensnode(q)+1;
            stop_r=nbsensnode(q)+Q;
        else
            for k=1:nb_children
                if ~isempty(find(Nu{k} == l))
                    start_r=nbsensnode(q)+ind(k)*Q+1;
                    stop_r=nbsensnode(q)+ind(k)*Q+Q;
                end
            end
        end
        G=diag(sign(diag(WGnew(start_r:stop_r,:))));
        W{l}=W{l}*G;
    end
    for l=q+1:nbnodes
        if ~ismember(l,descendents_q)
            start_r=nbsensnode(q)+1;
            stop_r=nbsensnode(q)+Q;
        else
            for k=1:nb_children
                ind_k=ind(k);
                if q==1
                   ind_k=ind_k-1; 
                end
                if ~isempty(find(Nu{k} == l))
                    start_r=nbsensnode(q)+ind_k*Q+1;
                    stop_r=nbsensnode(q)+ind_k*Q+Q;
                end
            end
        end
        G=diag(sign(diag(WGnew(start_r:stop_r,:))));
        W{l}=W{l}*G;
    end
    
    
end

function W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes,tree,Nu,descendents_q)

    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end
    
    children_q=[find(tree==q)];
    nb_children=length(children_q);
    ind=1:nb_children;
    
    for l=1:q-1
        if ~ismember(l,descendents_q)
            start_r=nbsensnode(q)+1;
            stop_r=nbsensnode(q)+Q;
        else
            for k=1:nb_children
                if ~isempty(find(Nu{k} == l))
                    start_r=nbsensnode(q)+ind(k)*Q+1;
                    stop_r=nbsensnode(q)+ind(k)*Q+Q;
                end
            end
        end
        W{l}=W{l}*WGnew(start_r:stop_r,:);
    end
    for l=q+1:nbnodes
        if ~ismember(l,descendents_q)
            start_r=nbsensnode(q)+1;
            stop_r=nbsensnode(q)+Q;
        else
            for k=1:nb_children
                ind_k=ind(k);
                if q==1
                   ind_k=ind_k-1; 
                end
                if ~isempty(find(Nu{k} == l))
                    start_r=nbsensnode(q)+ind_k*Q+1;
                    stop_r=nbsensnode(q)+ind_k*Q+Q;
                end
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

