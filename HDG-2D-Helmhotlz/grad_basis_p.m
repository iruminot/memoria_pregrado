function gb_p = grad_basis_p(j,lam1,lam2)
% Description: This function returns the gradient of j-th basis function
% ***********   
%
% The vector (lam1,lam2,1-lam1-lam2) is the vector of 
% barycentric coordinates in the reference triangle. This means
% that, for any point x in the triangle, we can write
%         x  =  lam1 * v1  +  lam2 * v2  + lam3 * v3
% where v1,v2, and v3 are the vertices of the reference triangle

gb_p(1) = basis_triangle_lam1(j,lam1,lam2);
gb_p(2) = basis_triangle_lam2(j,lam1,lam2);

