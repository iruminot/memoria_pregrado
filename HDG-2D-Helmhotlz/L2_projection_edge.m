function coef = L2_projection_edge(bedg,func,dim_M, NGP_edge, col_m_edg, w_edge, x_edge, ex)
% Computes the coeficients of the L2 projection of a function "func" 
% into the edge "bedg". The matrix is |e|*I.

%global dim_M NGP_edge col_m_edg w_edge x_edge ex

if bedg.edge_or == 1 
    v0 = bedg.vert(1,:);
    v1 = bedg.vert(2,:);
else
    v0 = bedg.vert(2,:);
    v1 = bedg.vert(1,:);
end

coef = zeros(dim_M,1);
for ii = 1 : dim_M
    coef(ii) = 0;
    for gp = 1 : NGP_edge
        vx = v0 + x_edge(1,gp)*(v1-v0);
        coef(ii) = coef(ii) +  w_edge(gp)*col_m_edg(ii,gp)*func(vx(1),vx(2),ex);
    end
end