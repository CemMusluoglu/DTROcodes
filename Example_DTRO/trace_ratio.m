function [W,rho]=trace_ratio(Q,M,Ruu,Rvv,denom_sum)
% Implementation of the centralized TRO
% INPUT: 
% M: Total
% Q: Number of filters
% Ruu, Rvv (M x M): covariance matrices
% denom_sum: Indicates whether to sum or not the denominator. If 1 then the
% matrix in the denominator is Ruu+Rvv otherwise it is Rvv.
%
% OUTPUTS:
% W (M x Q): Estimation of a solution of the TRO problem
% rho: Evaluation of the TRO at W

    i=0;
    winit=randn(M,Q);
    winit=normc(winit);
    rho=trace(winit'*Ruu*winit)/trace(winit'*denom_mat(Ruu,Rvv,denom_sum)*winit);
    rho_old=rho+1;
    tol_rho=1e-6;
    nbiter=-1;
    W=winit;
    while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)

        [B_int,s_int]=eig(Ruu-rho*denom_mat(Ruu,Rvv,denom_sum));
        [~,ind_int]=sort(diag(s_int),'descend');

        W=B_int(:,ind_int(1:Q));
        rho_old=rho;
        rho = trace(W'*Ruu*W)/trace(W'*denom_mat(Ruu,Rvv,denom_sum)*W);

        i=i+1;

    end

end

function R_denom=denom_mat(Ruu,Rvv,denom_sum)

    if denom_sum==1
        R_denom=Ruu+Rvv;
    else
        R_denom=Rvv;
    end
end

