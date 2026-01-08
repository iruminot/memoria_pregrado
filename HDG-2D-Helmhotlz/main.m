function [hmax n_elem error_u error_q error_lam ] = main(k,ref,geom,ex)

%--------------------------------------------------------------------------
% Global variables
%global mesh
%global elem bedg edges

%global method
%global  geom ex dim_V dim_W dim_M  deg

%global quadratures
%global NGP_edge NGP_trian x_t w_t x_edge w_edge

% global collocation
%global col_jacv col_w_trian col_v_trian col_m_edg col_w_edg col_v_edg



%dimension of local spaces
dim_V = (k+1)*(k+2);
dim_W = (k+1)*(k+2)/2;
dim_M = k + 1;

%exactness quadrature rule
deg = 2*k;

% defining the geometry
[elem,bedg,edges,n_vert,n_edge,n_bedg,hmax]  = mesh_generator(geom,ref);

n_elem = size(elem,2);

% evaluating basis funcitons at quadrature points
[col_jacv, col_w_trian, col_v_trian, NGP_trian, x_t, w_t] = collocation_triangles(dim_V, dim_W ,deg);
[col_m_edg col_w_edg col_v_edg, NGP_edge, x_edge, w_edge] = collocation_edges(dim_V, dim_W, dim_M, deg);


% building local matrices
[A_local,b_local,Au_store,Ap_store,Bu_store,Bp_store] = local_solvers(n_elem,elem, NGP_edge, NGP_trian, w_t, x_t,w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian);


%Assembling global system
A_global = sparse(n_edge*dim_M,n_edge*dim_M);
b_global = sparse(n_edge*dim_M,1);
d = k+1;
r = 1:d;
for i = 1 : n_elem
    block1 = (elem(i).edge_gn(1)-1).*dim_M+1: elem(i).edge_gn(1).*dim_M;
    block2 = (elem(i).edge_gn(2)-1).*dim_M+1: elem(i).edge_gn(2).*dim_M;
    block3 = (elem(i).edge_gn(3)-1).*dim_M+1: elem(i).edge_gn(3).*dim_M;
    block = [block1';block2';block3'];
    b_global(block) = b_global(block) + b_local{i};
    A_global(block,block) = A_global(block,block) + A_local{i};
end


int_d = [];

%Imposing boundary conditions
for j = 1 : n_bedg % loop over boundary edges
    if bedg(j).edge_or == 1 
        v0 = bedg(j).vert(1,:);
        v1 = bedg(j).vert(2,:);
    else
        v0 = bedg(j).vert(2,:);
        v1 = bedg(j).vert(1,:);
    end
    
    int_d = L2_projection_edge(bedg(j),@dirichlet_bc,dim_M, NGP_edge, col_m_edg, w_edge, x_edge, ex);
    
    A_global(r+(j-1)*d,:) = sparse(length(r+(j-1)*d),n_edge*dim_M);
    A_global(r+(j-1)*d,r+(j-1)*d) = speye(k+1);
    b_global(r+(j-1)*d) = int_d;
end


% Solving the lienar system
gamma = A_global\b_global;

%compute local coefficients
coefloc = coeflocal(Au_store, Ap_store, Bu_store, Bp_store, gamma,elem,k);


%compute errors
error_u = cell(n_elem,1);
error_q = cell(n_elem,1);
error_lam = cell(n_elem,1);
for j = 1 : n_elem
    
    %Computing error in q and u 
    lam = x_t';
    BT = elem(j).BT;
    bT = elem(j).bT;
    
    
    %[qh,uh] = evaluate_qh_uh(lam(1,:),lam(2,:),coefloc(j),dim_W, dim_V);
    [qh,uh] = evaluate_qh_uh_gp(coefloc(j),dim_W, dim_V,col_v_trian,col_w_trian,lam(1,:));
    x = BT*x_t' + repmat(bT,1,NGP_trian);
    
    u_ex = u_exact(x(1,:),x(2,:),ex);
    q_ex = q_exact(x(1,:),x(2,:),ex);
    
    error_u{j,1} =  elem(j).area*(u_ex - uh).^2*w_t';
    error_q{j,1}  = elem(j).area*sum((q_ex - qh).^2,1)*w_t';
    
    %Computing error of the projection of lambda (numerical trace of u)
    loc_error = 0;
    for ee = 1 : 3
        edge_gn = elem(j).edge_gn(ee);
        u_proj = L2_projection_edge(edges(edge_gn),@u_exact,dim_M, NGP_edge, col_m_edg, w_edge, x_edge, ex);
        loc_error = loc_error + sum((u_proj-coefloc(j).gammaF(:,ee)).^2).*elem(j).edge_l(ee);
    end
    error_lam{j,1}  = loc_error.*max(elem(j).edge_l);
end
error_lam = sqrt(sum([error_lam{:,1}]));
error_u = sqrt(sum([error_u{:,1}]));
error_q = sqrt(sum([error_q{:,1}]));

%Ploting solution
figure
hold on
for j = 1 : n_elem
    x = elem(j).vert(:,1)';
    y = elem(j).vert(:,2)';
    BT = elem(j).BT;
    bT = elem(j).bT;

    lam = BT\(([x;y])-repmat(bT,1,length(x)));
    
    [qh,uh] = evaluate_qh_uh(lam(1,:),lam(2,:),coefloc(j),dim_W,dim_V);
    trisurf([1 2 3],x,y,uh);
end

hold off