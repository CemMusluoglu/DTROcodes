function [W,rho]=trace_ratio_simple(Ruu,Rvv,Q,M)

i=0;
winit=randn(M,Q);
winit=normc(winit);
rho=trace(winit'*Ruu*winit)/trace(winit'*(Rvv)*winit);
rho_old=rho+1;
tol_rho=1e-6;
nbiter=-1;
W=winit;
while (tol_rho>0 && abs(rho-rho_old)>tol_rho) || (i<nbiter)
    
    
    [B_int,s_int]=eig(Ruu-rho*(Rvv));
    [~,ind_int]=sort(diag(s_int),'descend');
    
    W=B_int(:,ind_int(1:Q));
    rho_old=rho;
    rho = trace(W'*Ruu*W)/trace(W'*(Rvv)*W);

    i=i+1;
   
end


end
