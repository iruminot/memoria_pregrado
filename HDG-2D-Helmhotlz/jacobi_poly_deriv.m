function jp_d =jacobi_poly_deriv(alpha,beta,n,x)
%------------------------------------------------------
% Description: This function returns the derivative of the 
% ***********  Jacobi polynomial P^{alpha, beta}_{n} evaluated at x.
%------------------------------------------------------

if n == 0
    jp_d = 0;
else
    jp_d = 0.5.*(n+alpha+beta+1.).*jacobi_poly(alpha+1,beta+1,n-1,x);
end