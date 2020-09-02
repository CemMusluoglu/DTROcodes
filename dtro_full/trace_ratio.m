function [W,rho_track,norm_track,norm_star_track]=trace_ratio(params,data,conv,W_star)

Q=params.Q;
nbsens=params.nbsens;
denom_sum=params.denom_sum;

Ruu=data.R_first;
Rvv=data.R_second;

tol_rho=conv.tol_rho;
nbiter=conv.nbiter;
i=0;
winit=randn(nbsens,Q);
winit=normc(winit);
rho=trace(winit'*Ruu*winit)/trace(winit'*denom_mat(Ruu,Rvv,denom_sum)*winit);
rho_old=rho+1;
rho_track=[];
norm_track=[];
norm_star_track=[];
W=winit;
while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    Wold=W;
    
    [B_int,s_int]=eig(Ruu-rho*denom_mat(Ruu,Rvv,denom_sum));
    [~,ind_int]=sort(diag(s_int),'descend');
    
    W=B_int(:,ind_int(1:Q));
    rho_old=rho;
    rho = trace(W'*Ruu*W)/trace(W'*denom_mat(Ruu,Rvv,denom_sum)*W);
    rho_track=[rho_track,rho];
    i=i+1;
    if i>1
        norm_track=[norm_track,norm(W-Wold,'fro')^2/numel(W)];
    end
    if nargin>3
        norm_star_track=[norm_star_track,norm(W-W_star,'fro')^2/numel(W_star)];
    end
end


end

function R_denom=denom_mat(Ruu,Rvv,denom_sum)

    if denom_sum==1
        R_denom=Ruu+Rvv;
    else
        R_denom=Rvv;
    end
end
