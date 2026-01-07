function [qh,uh] = evaluate_qh_uh(lam1,lam2,coeflocT,dim_W, dim_V)
% Evaluates the discrete solution at the baricentric coordinate (lam1,lam2)

%global dim_W dim_V 
alphaT = coeflocT.alpha;
betaT = coeflocT.beta;


qh=zeros(2,length(lam1));
uh=zeros(1,length(lam1));
for j = 1 : dim_V
    qh = qh + alphaT(j).*basis_v(j,lam1',lam2').';
end

for j = 1 : dim_W
    uh = uh + betaT(j).*basis_w(j,lam1',lam2')';
end

