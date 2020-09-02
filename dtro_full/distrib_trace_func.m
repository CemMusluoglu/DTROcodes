function [w,rho_track,norm_track,norm_star_track]=distrib_trace_func(params,data,conv,debug,W_star)

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
norm_track=[];
norm_star_track=[];

while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    w_old=w;
    
    q=rem(i,nbnodes)+1;
    
    C_k=constr_C(W,Q,q,nbsensnode,nbnodes);
    
    [Ruucompressed,Rvvcompressed,Compressor]=compress(Ruu,Rvv,C_k);
    
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
    
    W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes);
    
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
        norm_track=[norm_track,norm(w-w_old,'fro')^2/numel(w)];
    end
    
    if nargin>4
        
        for l=1:Q
            if sum(sum((W_star(:,l)-w(:,l)).^2))>sum(sum((-W_star(:,l)-w(:,l)).^2))
                w(:,l)=-w(:,l);
            end
        end
        
        norm_star_track=[norm_star_track,norm(w-W_star,'fro')^2/numel(W_star)];
    end

end

end


function C=constr_C(W,Q,q,nbsensnode,nbnodes)
    B_less=blkdiag(W{1:q-1});
    B_more=blkdiag(W{q+1:nbnodes});
    C=[zeros(sum(nbsensnode(1:q-1)),nbsensnode(q)),B_less,zeros(sum(nbsensnode(1:q-1)),...
            Q*(nbnodes-q)) ; eye(nbsensnode(q)),zeros(nbsensnode(q),Q*(nbnodes-1)) ;...
           zeros(sum(nbsensnode(q+1:nbnodes)),nbsensnode(q)+Q*(q-1)),B_more];
end

function [Ruucompressed,Rvvcompressed,Compressor]=compress(Ruu,Rvv,C_k)

    Ruucompressed=C_k'*Ruu*C_k;
    Rvvcompressed=C_k'*Rvv*C_k;
    Ruucompressed=make_sym(Ruucompressed);
    Rvvcompressed=make_sym(Rvvcompressed);

    Compressor=C_k'*C_k;
    Compressor=make_sym(Compressor);

end

function WG=comp_WG(Ruucompressed,rho,Rvvcompressed,Compressor,Q,denom_sum)
    [B,s]=eig(Ruucompressed-rho*denom_mat(Ruucompressed,Rvvcompressed,denom_sum),Compressor);
    [~,ind]=sort(diag(s),'descend');
    WG=B(:,ind(1:Q));
end

function W=update_W_efficient(W,WGnew,Q,q,nbsensnode,nbnodes)
    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end
    
    teller=1;
    for l=1:q-1
        G=diag(sign(diag(WGnew(nbsensnode(q)+teller:nbsensnode(q)+teller+Q-1,:))));
        W{l}=W{l}*G;
        teller=teller+Q;
    end
    for l=q+1:nbnodes
        G=diag(sign(diag(WGnew(nbsensnode(q)+teller:nbsensnode(q)+teller+Q-1,:))));
        W{l}=W{l}*G;
        teller=teller+Q;
    end
end

function W=update_W(W,WGnew,Q,q,nbsensnode,nbnodes)

    Wqold=W{q};
    W{q}=WGnew(1:nbsensnode(q),:);
    
    for l=1:Q
        if sum(sum((Wqold(:,l)-W{q}(:,l)).^2))>sum(sum((-Wqold(:,l)-W{q}(:,l)).^2))
            W{q}(:,l)=-W{q}(:,l);
            WGnew(:,l)=-WGnew(:,l);
        end
    end
    
    teller=1;
    for l=1:q-1
        W{l}=W{l}*WGnew(nbsensnode(q)+teller:nbsensnode(q)+teller+Q-1,:);
        teller=teller+Q;
    end
    for l=q+1:nbnodes
        W{l}=W{l}*WGnew(nbsensnode(q)+teller:nbsensnode(q)+teller+Q-1,:);
        teller=teller+Q;
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


