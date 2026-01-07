 function bt = basis_triangle(j,lam1,lam2)
%-----------------------------------------------------------------
%     Description: This function returns the j-th orthogonal basis function
%     ***********  on the reference triangle.

%     The vector (lam1, lam2, 1-lam1-lam2) is the vector of
%     barycentric coordinates in the reference triangle.
%-----------------------------------------------------------------

r = 2.*lam1-1;
s = 2.*lam2-1;

bb = (1+s)./2;
aa = 2.*lam1./(1.-lam2)-1;
x = -1.5+.5.*sqrt(1.+8.*j);
nm = ceil(x);
nm1=(nm+1).*(nm+2);
n = nm1./2-j;
m = nm-n;
theta_j = sqrt(2.0.*(2.0*n+1).*(m+n+1.0));

qr = zeros(length(lam1),1);
for ii = 1 : length(lam1)
    if lam2(ii) == 1      
        for i = 0 : n
            qr(ii) = qr(ii) + (factorial(n)./(factorial(n-i).*factorial(i))).^2.*(2.*lam1(ii)).^n;
        end
    else
        qr(ii) = jacobi_poly(0,0,n,aa(ii)).*(1.-bb(ii)).^n;
    end
    
end
bt = jacobi_poly(2*n+1,0,m,s).*qr.*theta_j;
