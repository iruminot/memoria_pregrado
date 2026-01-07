function div = r2_div(i_elem,n_qpoints_triangle, dim_V, col_jacv,elem)
%*********************************************************************
%
%-------------------------------------------------------------------------------
%     Description: This subroutine evaluates the divergence of the velocity 
%     ***********  basis functions at the quadrature points of the  "i_elem"-th 
%     triangle. 
%-------------------------------------------------------------------------------
%     getting transformation gradients
%     --------------------------------
      %global elem

      grad_lam_1 = (-1/(elem(i_elem).area.*2)).*elem(i_elem).norm(1,:);
      grad_lam_2 = (-1/(elem(i_elem).area.*2)).*elem(i_elem).norm(2,:);

  
%     computing the divergence
%     ------------------------
      
      for i = 1 : dim_V %loop on the velocity basis functions
          for j = 1 : n_qpoints_triangle %loop on the quadrature points divergence through chain rule
              div(i,j) =  grad_lam_1*col_jacv{i,j}(:,1) + grad_lam_2*col_jacv{i,j}(:,2);
          end
      end