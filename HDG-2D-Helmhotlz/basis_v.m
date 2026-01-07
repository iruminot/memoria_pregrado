function bv = basis_v(j,lam1,lam2)
%-------------------------------------------------------------------
% Description: This function returns the j-th basis function
% ***********  for the velocity inside the reference 
% The vector (lam1,lam2,1-lam1-lam2) is the vector of 
% barycentric coordinates in the reference triangle. This means
% that, for any point x in the triangle, we can write
%         x  =  lam1 * v1  +  lam2 * v2  + lam3 * v3
%
% where v1,v2, and v3 are the vertices of the reference triangle
bv = sparse(length(lam1),2);
if mod(j,2)== 0
   bv(:,2) = basis_triangle(j/2,lam1,lam2);
else
    bv(:,1) = basis_triangle((j+1)/2,lam1,lam2);
end