function bt_d2 = basis_triangle_lam2(j,lam1,lam2)
%-----------------------------------------------------------------
%Description: This function returns the derivative w.r.t. lam2 of the j-th orthogonal basis function
%     ***********  on the reference triangle.

%     The vector (lam1, lam2, 1-lam1-lam2) is the vector of
%     barycentric coordinates in the reference triangle.
%-----------------------------------------------------------------

r = 2.*lam1-1;
s = 2.*lam2-1;

bb = (1+s)./2;
aa = 2.*lam1/(1.-lam2)-1;
x = -1.5+.5*sqrt(1.+8.*j);
nm = ceil(x);
nm1=(nm+1)*(nm+2);
n = nm1/2-j;
m = nm-n;
theta_j = sqrt(2.0*(2.0*n+1)*(m+n+1.0));

AAA1=jacobi_poly(2*n+1,0,m,s).*(1.-bb).^n.*jacobi_poly_deriv(0,0,n,aa).*2.*lam1.*(1-lam2).^(-2);
AAA2=jacobi_poly_deriv(2*n+1,0,m,s).*(1.-bb).^n.*jacobi_poly(0,0,n,aa).*2;
AAA3=jacobi_poly(2*n+1,0,m,s).*jacobi_poly(0,0,n,aa).*n.*(1.-bb).^(n-1.).*(-1.);
bt_d2 = (AAA1+AAA2+AAA3).*theta_j;
