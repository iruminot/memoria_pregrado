function [KK_store,BB_store,Au_store,Ap_store,Bu_store,Bp_store] = local_solvers(n_elem,elem, NGP_edge, NGP_trian, w_t, x_t, w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian)
%----------------------------------------------------------------------------------
%     Description: This module compute for each element (iele)
%                  the local matrices and the local load vector
%                  associated to the body forces, for the
%                  subsequent global assembling.
%----------------------------------------------------------------------------------
%     Programmers : Manuel Solano
%     Date        : August, 2015
%----------------------------------------------------------------------------------



    % local stiffness matrix & RHS stored for all element, storing structure KK_store (:,:,n_elem)
    KK_store = cell(n_elem);% dimension of a element in the cell: (3*dim_M,3*dim_M)
    BB_store = cell(n_elem);% dimension of a element in the cell: (3*dim_M,3*dim_M)
    Au_store = cell(n_elem);% dimension of a element in the cell: (dim_V,3*dim_M)
    Ap_store = cell(n_elem);% dimension of a element in the cell: (dim_W,3*dim_M)
    Bu_store = cell(n_elem);% dimension of a element in the cell: (dim_V,1)
    Bp_store = cell(n_elem);% dimension of a element in the cell: (dim_V,1)


    %parfor iele = 1 : n_elem
    for iele = 1 : n_elem
        [Au,Ap,Bu,Bp,mu_Ph_hat, w_Ph_hat, vn_Ph_hat] = subr_solver(iele,elem, NGP_edge, NGP_trian, w_t, x_t, w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian);

        mu_uh_n = vn_Ph_hat.';
        mu_Ph = w_Ph_hat.';

        KK = mu_uh_n*Au + mu_Ph*Ap - mu_Ph_hat;
        BB = mu_uh_n*Bu + mu_Ph*Bp;

        KK_store{iele,1} = KK;
        BB_store{iele,1} = -BB(:,:);
        Au_store{iele,1}=  Au;
        Ap_store{iele,1}=  Ap;
        Bu_store{iele,1}=  Bu;
        Bp_store{iele,1}=  Bp;
    end
end


%////////////////////////////////////////////////////////////////////////////////
function  [Au,Ap,Bu,Bp,mu_Ph_hat, w_Ph_hat, vn_Ph_hat] = subr_solver(iele,elem, NGP_edge, NGP_trian, w_t, x_t, w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian)
%**************************************************************************
%     Description: This subroutine obtain the relationship between local dof and
%                  multiplier dof. It also obtain the relationship of local dof
%                  and the source term f

%---------------------------------------------------------------------------------
%   Global variables
%    global dim_V dim_W dim_M
%---------------------------------------------------------------------------------

      [Ke, Fp, Fb,mu_Ph_hat, w_Ph_hat, vn_Ph_hat] = subr_assembly(iele,elem, NGP_edge, NGP_trian, w_t, x_t,w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian);
%---------------------------------------------------------------------------------
%     begin matrix multiplications 
%---------------------------------------------------------------------------------
      Ae = Ke\Fp;
      Be = Ke\Fb;
      
%---------------------------------------------------------------------------------
%     begin extracting corresponding matrices for u and p 
%---------------------------------------------------------------------------------
      Au = Ae(1 : dim_V,1 : 3*dim_M);
      Bu = Be(1 : dim_V,1);
      Ap = Ae(dim_V+1 : dim_V+dim_W, 1 : 3*dim_M);
      Bp = Be(dim_V+1 : dim_V+dim_W,1);
end 
%////////////////////////////////////////////////////////////////////////////////


%////////////////////////////////////////////////////////////////////////////////
function [Ke, Fp, Fb,mu_Ph_hat, w_Ph_hat, vn_Ph_hat] = ...
          subr_assembly(iele,elem, NGP_edge, NGP_trian, w_t, x_t,w_edge, ex, ...
          dim_V, dim_W, dim_M, col_m_edg,col_w_edg,col_v_edg,col_jacv, col_w_trian, ...
          col_v_trian)
%  ****************************************
%     Description: This subroutine assemble the local solver stiffness matrix and
%                  the RHS
%---------------------------------------------------------------------------------
%---------------------------------------------------------------------------------
%   Global variables
%---------------------------------------------------------------------------------
%     call functions

      mu_Ph_hat = subr_muPhhat(iele,NGP_edge, ex, dim_M, elem, w_edge,col_m_edg);
      w_Ph_hat = subr_wPh_hat(iele,NGP_edge, dim_W, dim_M, elem, w_edge, ex,col_m_edg, col_w_edg);
      vn_Ph_hat = subr_vnPh_hat(iele,NGP_edge, dim_V, dim_M, elem,w_edge,col_m_edg, col_v_edg);
	  Cv_uh = subr_Cvuh(iele,NGP_trian, dim_V, elem, w_t, x_t, ex,col_v_trian);
	  w_Ph = subr_wPh(iele,dim_W, ex, elem, w_edge,col_w_edg);
	  div_v_Ph =  subr_divv_Ph(iele,NGP_trian, dim_V, elem, w_t,col_jacv, col_w_trian);
	  w_f = subr_wf(iele,elem, w_t, x_t, ex,col_w_trian, NGP_trian);
      uv = subr_uhvh(iele,NGP_trian, dim_W, elem, w_t, x_t, ex,col_w_trian);
%---------------------------------------------------------------------------------

      Ke = sparse(dim_V+dim_W,dim_V+dim_W);     
      Ke(1 : dim_V,1 : dim_V) = Cv_uh(1 : dim_V,1 : dim_V);
      Ke(1 : dim_V,dim_V+1 : dim_V+dim_W) = -div_v_Ph(:,1 : dim_W);

	  div_v_Ph_TRANS = div_v_Ph.';
      Ke(dim_V+1 : dim_V+dim_W,1 : dim_V) =  div_v_Ph_TRANS(1 : dim_W,1 : dim_V);
      Ke(dim_V+1 : dim_V+dim_W,dim_V+1 : dim_V+dim_W) =  w_Ph(1 : dim_W,1 : dim_W) - uv;
      
%---------------------------------------------------------------------------------

      Fp = sparse(dim_V+dim_W,3*dim_M);
      Fp(1 : dim_V,1 : 3*dim_M) = -vn_Ph_hat;
      Fp(dim_V+1 : dim_V+dim_W,1 : 3*dim_M) = w_Ph_hat;

%---------------------------------------------------------------------------------

      Fb = sparse(dim_V+dim_W,1);
      Fb( dim_V+1 : dim_V+dim_W,1) = w_f;
      
end
%////////////////////////////////////////////////////////////////////////////////


%/////////////////////////////////////////////////////////////////////////////////
function mu_Ph_hat = subr_muPhhat(iele,NGP_edge, ex, dim_M, elem, w_edge,col_m_edg)
%  ********************************************************************************
%     Description: This subroutine computes the matrix tau<mu, Ph_hat>
%---------------------------------------------------------------------------------

    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_edge ex dim_M elem w_edge
    %    global col_m_edg
    %---------------------------------------------------------------------------------

    t = tau(ex);

    %	 begin integration of tau< mu, Ph_hat >
    for iedge = 1 : 3							% loop over edges
        orien = elem(iele).edge_or(iedge);      % orientation of iedge
        len = elem(iele).edge_l(iedge); 		% length of iedge    
        gp_vec = [1:NGP_edge].*(orien+1)./2+(NGP_edge+1-[1:NGP_edge]).*(1-orien)./2;
        component = (col_m_edg(:,gp_vec).*repmat(w_edge',dim_M,1))*col_m_edg(:,gp_vec).';
        if iedge == 1
            mu_Ph_hat(1:dim_M,1:dim_M)   = t(iedge)*component*len;
        elseif iedge == 2
            mu_Ph_hat(dim_M+1 : 2*dim_M,dim_M+1 : 2*dim_M) = t(iedge)*component*len;
        else
            mu_Ph_hat(2*dim_M+1 : 3*dim_M,2*dim_M+1 : 3*dim_M) = t(iedge)*component*len;
        end   
    end
end



%///////////////////////////////////////////////////////////////////////////////////
function aux = subr_Cvuh(iele,NGP_trian, dim_V, elem, w_t, x_t, ex,col_v_trian)
%******************************************************************************
%     Description: This subroutine computes (C v, uh)
%---------------------------------------------------------------

    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_trian dim_V elem w_t x_t ex
    %    global col_v_trian
    %---------------------------------------------------------------------------------

%	 begin integration of (C v, uh)

    BT = elem(iele).BT;
    bT = elem(iele).bT;
    
    x = BT*x_t' + repmat(bT,1,NGP_trian);
    c = tensor_C(x,ex); 
    
    rep = repmat(c.'.*w_t,dim_V,1);
    aux = col_v_trian(:,:,1)*(col_v_trian(:,:,1).*rep).' + col_v_trian(:,:,2)*(col_v_trian(:,:,2).*rep).';
    aux = aux.*elem(iele).area;
end

%////////////////////////////////////////////////////////////////////////////////
function  aux = subr_wPh(iele,dim_W, ex, elem, w_edge,col_w_edg)
%***************************************************************
%     Description: This subroutine computes  tau < w , Ph >
%--------------------------------------------------------------------------------

    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_edge dim_W ex elem w_edge
    %    global col_w_edg 
    %---------------------------------------------------------------------------------

    t = tau(ex);

    aux = sparse(dim_W,dim_W);                         % initialize w_Ph for each element
    
    for iedge = 1 : 3					% loop over iedge
        len = elem(iele).edge_l(iedge); % length of each edge
        component = (col_w_edg(:,:,iedge).*repmat(w_edge',dim_W,1))*col_w_edg(:,:,iedge).';
        aux = aux +  t(iedge).*component.*len;
    end
end


%////////////////////////////////////////////////////////////////////////////////
function  w_Ph_hat =  subr_wPh_hat(iele,NGP_edge, dim_W, dim_M, elem, w_edge, ex,col_m_edg, col_w_edg)
%  ***************************************************************************************************
%     Description: This subroutine computes tau < w , Ph_hat >
%---------------------------------------------------------------


    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_edge dim_W dim_M elem w_edge ex
    %    global col_m_edg col_w_edg
    %---------------------------------------------------------------------------------

    t = tau(ex);
    %   begin integration of  tau < w , Ph_hat >

    for iedge = 1 : 3								% loop over edges
        orien = elem(iele).edge_or(iedge);				% determine the orien of each edge
        len = elem(iele).edge_l(iedge);				% length of each edge
        
        gp_vec = [1:NGP_edge].*(orien+1)./2+(NGP_edge+1-[1:NGP_edge]).*(1-orien)./2;
        component = (col_m_edg(:,gp_vec).*repmat(w_edge',dim_M,1))*col_w_edg(:,:,iedge).';
        
        if iedge == 1
            w_Ph_hat(1:dim_W, 1:dim_M) = t(iedge).*component.'.*len;
        elseif iedge == 2
            w_Ph_hat(1:dim_W, dim_M+1:2*dim_M) = t(iedge).*component.'.*len;
        else
            w_Ph_hat(1:dim_W, 2*dim_M+1:3*dim_M) = t(iedge).*component.'.*len;
        end
    end
    ss=1;

end


%////////////////////////////////////////////////////////////////////////////////
function  vn_Ph_hat = subr_vnPh_hat(iele,NGP_edge, dim_V, dim_M, elem,w_edge,col_m_edg, col_v_edg)
%*************************************************************************************************
%     Description: This subroutine computes < v.n , Ph_hat >
%---------------------------------------------------------------


    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_edge dim_V dim_M elem w_edge
    %    global col_m_edg col_v_edg
    %---------------------------------------------------------------------------------


    for iedge = 1 : 3									% loop over edges
        normals = elem(iele).norm(iedge,:);				% determine the normals vector on each edge
        len = elem(iele).edge_l(iedge);					% length of each edge;
        orien = elem(iele).edge_or(iedge);					% determine the orien of each edge

        gp_vec = [1:NGP_edge].*(orien+1)./2+(NGP_edge+1-[1:NGP_edge]).*(1-orien)./2;
        mult = col_v_edg{iedge}(:,:,1).*normals(1)+col_v_edg{iedge}(:,:,2).*normals(2);
        component = (col_m_edg(:,gp_vec).*repmat(w_edge',dim_M,1))*mult';        

        
        if iedge == 1
            vn_Ph_hat(1:dim_V,1:dim_M) = component.';
        elseif iedge == 2
            vn_Ph_hat(1:dim_V,dim_M+1:2*dim_M) = component.';
        else
            vn_Ph_hat(1:dim_V,2*dim_M+1:3*dim_M) = component.';
        end             

    end
    dd=1;

end

%////////////////////////////////////////////////////////////////////////////////
function aux = subr_divv_Ph(iele,NGP_trian, dim_V, elem, w_t,col_jacv, col_w_trian)
%**********************************************************************************
%     Description: This subroutine computes ( div v , Ph )
%---------------------------------------------------------------


    %---------------------------------------------------------------------------------
    %   Global variables
    %    global NGP_trian dim_V dim_W elem w_t
    %    global col_jacv col_w_trian
    %---------------------------------------------------------------------------------

    div = r2_div(iele,NGP_trian, dim_V, col_jacv,elem);

    %----------------------------------------------------------------------------------
    %	 begin integration of  ( div v , Ph )

    aux = (elem(iele).area.*div.*repmat(w_t,dim_V,1))*col_w_trian';
end

%////////////////////////////////////////////////////////////////////////////////
function aux = subr_wf(iele,elem, w_t, x_t, ex,col_w_trian, NGP_trian)
%*********************************************************************
%     Description: This subroutine computes ( w , f )
%---------------------------------------------------------------

    %---------------------------------------------------------------------------------
    %   Global variables
    %    global elem w_t x_t ex
    %    global col_w_trian NGP_trian
    %---------------------------------------------------------------------------------

 
    % begin integration of  ( w , f )

    BT = elem(iele).BT;
    bT = elem(iele).bT;
    
    x = BT*x_t' + repmat(bT,1,NGP_trian);
    f_val = f(x(1,:),x(2,:),ex);
    aux = col_w_trian*(f_val.*w_t).'.*elem(iele).area;
end

function aux = subr_uhvh(iele,NGP_trian, dim_W, elem, w_t, x_t, ex,col_w_trian)
%******************************************************************************
%     Description: This subroutine computes omega^2(uh, vh)
%---------------------------------------------------------------


%	 begin integration of omega.^2(v, uh)   
    om = omega(ex);
    rep = repmat(w_t,dim_W,1);
    aux = col_w_trian*(col_w_trian.*rep).' ;
    aux = om.^2.*aux.*elem(iele).area;
end
