function [qh,uh] = evaluate_qh_uh_gp(coeflocT,dim_W, dim_V,col_v_trian,col_w_trian,lam1)
% Evaluates the discrete solution at the quadrature points

%global dim_W dim_V 
alphaT = coeflocT.alpha;
betaT = coeflocT.beta;


qh=zeros(2,length(lam1));
uh=zeros(1,length(lam1));
for j = 1 : dim_V
    aux = [col_v_trian(j,:,1); col_v_trian(j,:,2)];
    qh = qh + alphaT(j).*aux;
end

for j = 1 : dim_W
    uh = uh + betaT(j).*col_w_trian(j,:);
end

