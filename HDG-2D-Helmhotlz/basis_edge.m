function be = basis_edge(j,lam)
% Description: This function returns the j-th orthogonal basis function
% ***********  on an edge of the reference triangle.
%
% The vector (lam, 1-lam) is the vector of barycentric 
% coordinates on the edge.
%-----------------------------------------------------------------
 
%be = legendre(j-1,2*lam-1).*sqrt(.5*(2.*j-1));
be = legendre(j-1,2*lam-1)/sqrt(1/(2*(j-1)+1));