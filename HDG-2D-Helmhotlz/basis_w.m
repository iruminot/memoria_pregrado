function bw = basis_w(j,lam1,lam2)
%-------------------------------------------------------------------
% Description: This function returns the j-th basis function
% ***********  for the pressure inside the reference triangle.

% The vector (lam1, lam2, 1-lam1-lam2) is the vector of
% barycentric coordinates in the reference triangle. This means
% that, for any point x in the triangle, we can write
%         x  =  lam1 * v1  +  lam2 * v2  + lam3 * v3
%where v1,v2, and v3 are the vertices of the reference triangle

bw = basis_triangle(j,lam1,lam2);
