function [col_jacv, col_w_trian, col_v_trian, n_quad_triangle, x_t, w_t] = collocation_triangles(dim_V, dim_W ,deg)
%*****************************************************************************************************************
%---------------------------------------------------------------
%Description: This subroutine evaluates all basis functions 
%***********  and the Jacobian matrix at the quadrature 
% 		   points on triangles. 
% 
%     The quadrature rule is exact for polynomials of degree "deg". 
%---------------------------------------------------------------

%global  dim_V dim_W deg

%n_quad_triangle = n_quad_rule(deg,2); % 2 for triangles = 2-simplexes
%n_qpuad_triangle_k1 = n_quad_rule(deg+2,2); % 2 for triangles = 2-simplexes

rule = ceil(0.5*(deg-1));

[ w_t, x_t ] = gm_rule_set(rule,2);
x_t = x_t';
n_quad_triangle = size(x_t,1); 

%evaluating "col_w_trian"
%-------------------------
col_w_trian = zeros(dim_W,n_quad_triangle);
for i = 1 : dim_W %loop on the basis functions
    col_w_trian(i,:) =  basis_w(i,x_t(:,1),x_t(:,2)).';
end

%evaluating "col_v_trian" and "col_jacv"
%---------------------------------------
%col_v_trian =  cell(dim_V,n_quad_triangle);
col_v_trian =  zeros(dim_V,n_quad_triangle,2);
col_jacv =  cell(dim_V,n_quad_triangle);

for i = 1 : dim_V %loop on the basis functions
    aux =  basis_v(i,x_t(:,1),x_t(:,2)).';
    col_v_trian(i,:,1) = aux(1,:);
    col_v_trian(i,:,2) = aux(2,:);
    for j = 1 : n_quad_triangle %loop on the quadrature points
        %col_v_trian{i,j} = basis_v(i,x_t(j,1),x_t(j,2));
        aux = basis_jacv(i,x_t(j,1),x_t(j,2));
        col_jacv{i,j} = aux.'; 
    end
end


% %evaluating "col_grad_p"
% %----------------------
% for i =1 : dim_W_k1% loop on the basis functions
%     for j =1 : n_qpoints_triangle_k1 %loop on the quadrature points
%         col_grad_p(i,j) = grad_basis_p(i,x_t(j,1),x_t(j,2));
%     end %loop on the quadrature points
% end %loop on the basis functions