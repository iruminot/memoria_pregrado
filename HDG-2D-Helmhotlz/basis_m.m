function bm = basis_m(j,par)
%-----------------------------------------------------------------
% Description: This function returns the j-th basis function
% ***********  for the pressure on the reference edge.
% The parameter "par", parametrizes the edge as follows. Any 
% point x of the edge is given by the expression
%       x  =  (1-par) * v1 + par * v2 
% Note that (1-par,par) is the vector of barycentric
% coordinates of the point x.

%bm = sqrt(2).*basis_edge(j,par);
bm = basis_edge(j,par);
