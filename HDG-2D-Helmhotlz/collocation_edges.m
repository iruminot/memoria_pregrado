function [col_m_edg col_w_edg col_v_edg n_quad_edge x_edge w_edge] = collocation_edges(dim_V, dim_W, dim_M, deg)
%-----------------------------------------------------------------------------------
% Description: Evaluates the basis functions of the edges at the 
%              quadrature points
%-----------------------------------------------------------------------------------

%global  dim_V dim_W dim_M deg

%getting the quadrature rule
%---------------------------
n_quad_edge = ceil(0.5*(deg+1));
[x,w_edge] = lgwt(n_quad_edge,0,1);

x_edge = [1-x x]';

%evaluating "col_m_edg"
%----------------------
col_m_edg =  zeros(dim_M,n_quad_edge);
for i = 1 : dim_M %loop on the basis functions
    col_m_edg(i,:) = basis_m(i,x_edge(1,:));
end 

%evaluating "col_w_edg"
%----------------------
col_w_edg =  zeros(dim_W,n_quad_edge,3);
for i = 1 : dim_W %loop on the basis functions
    col_w_edg(i,:,1) = basis_w(i,zeros(n_quad_edge,1),x_edge(2,:)');
    col_w_edg(i,:,2) = basis_w(i,x_edge(1,:)',zeros(n_quad_edge,1));
    col_w_edg(i,:,3) = basis_w(i,x_edge(2,:)',x_edge(1,:)');
end

%
%evaluating "col_v_edg"
%----------------------
col_v_edg =  cell(1,1);
for i = 1 : dim_V %loop on the basis functions
    col_v_edg{1}(i,:,1:2) = full(basis_v(i,zeros(n_quad_edge,1),x_edge(2,:)'));
    col_v_edg{2}(i,:,1:2) = full(basis_v(i,x_edge(1,:)',zeros(n_quad_edge,1)));
    col_v_edg{3}(i,:,1:2) = full(basis_v(i,x_edge(2,:)',x_edge(1,:)'));
end