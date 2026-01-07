function bt_d1 = basis_triangle_lam1(j,lam1,lam2)
%Description: This function returns the derivative w.r.t. lam1 of the j-th orthogonal basis function
%***********  on the reference triangle.
%     The vector (lam1, lam2, 1-lam1-lam2) is the vector of
%     barycentric coordinates in the reference triangle.
%-----------------------------------------------------------------
r = 2.*lam1-1;
s = 2.*lam2-1;

bb = (1+s)./2;
aa = 2.*lam1./(1.-lam2)-1;
x = -1.5+.5*sqrt(1.+8.*j);
nm = ceil(x);
nm1=(nm+1)*(nm+2);
n = nm1/2-j;
m = nm-n;
theta_j = sqrt(2.0*(2.0*n+1)*(m+n+1.0)); 

bt_d1 = jacobi_poly_deriv(0,0,n,aa).*jacobi_poly(2*n+1,0,m,s).*(1.-bb).^n.*theta_j.*2./(1.-lam2);
