function [W,rho]=trace_ratio(Q,M,Ruu,Rvv,)
% Implementation of the centralized TRO
% INPUT: 
% M: Total
% Q: Number of filters
% Ruu, Rvv (M x M): covariance matrices
%
% OUTPUTS:
% W (M x Q): Estimation of a solution of the TRO problem
% rho: Evaluation of the TRO at W

i=0;
winit=randn(M,Q);
winit=normc(winit);
rho=trace(winit'*Ruu*winit)/trace(winit'*(Rvv+Ruu)*winit);
rho_old=rho+1;
tol_rho=1e-6;
nbiter=-1;
W=winit;
while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    [B_int,s_int]=eig(Ruu-rho*(Rvv+Ruu));
    [~,ind_int]=sort(diag(s_int),'descend');
    
    W=B_int(:,ind_int(1:Q));
    rho_old=rho;
    rho = trace(W'*Ruu*W)/trace(W'*(Rvv+Ruu)*W);

    i=i+1;
   
end

end

