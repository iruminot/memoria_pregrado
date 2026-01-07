function lambdah = evaluate_lamh(lam1,coeflocT,dim_M)
% Evaluates the discrete solution lambda_h at the quadrature pont lam1

%global dim_M


gammaF1 = coeflocT.gammaF1;
gammaF2 = coeflocT.gammaF2;
gammaF3 = coeflocT.gammaF3;

lambdah = zeros(3,1);
for j = 1 : dim_M
    basis_m(j,lam1)'
    lambdah(1) = lambdah(1) + gammaF1(j).*basis_m(j,lam1)';
    lambdah(2) = lambdah(2) + gammaF2(j).*basis_m(j,lam1)';
    lambdah(3) = lambdah(3) + gammaF3(j).*basis_m(j,lam1)';
end