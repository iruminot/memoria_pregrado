function J = basis_jacv(j,lam1,lam2)
% Description: This function returns the j-th orthogonal basis function
% ***********  for the velocity jacobian inside the reference triangle.
% The vector (lam1,lam2,1-lam1-lam2) is the vector of 
% barycentric coordinates in the reference triangle.

if mod(j,2)==1
    J(1,1) = basis_triangle_lam1((j+1)/2,lam1,lam2);
	J(1,2) = 0;
	J(2,1) = basis_triangle_lam2((j+1)/2,lam1,lam2);
	J(2,2) = 0;
else
     J(1,1) = 0;
	J(1,2) = basis_triangle_lam1(j/2,lam1,lam2);
	J(2,1) = 0;
	J(2,2) = basis_triangle_lam2(j/2,lam1,lam2);
end