function [elem,bedg,edges,n_vert,n_edge,n_bedg,hmax] = mesh_generator(geom,ref)
%-------------------------------------------------------------------------
% Description: Generates the mesh
% Programmer: Manuel Solano
% Date: Sept 6, 2013
%-------------------------------------------------------------------------

%ref = 2; %counter to refine the mesh

%ind = 1;
%while ind == 1
    hmax = 0;
    
    if geom==1    
        [vertices,elements,edges] = generate_rectangle(1,1,ref);  
        %edge = [];
        % load(strcat('mesh',num2str(ref),'.mat'));
        %elements = elem_back;
        %vertices = coord;        
        %edges = edge;
    elseif geom==2   
        [vertices,elements,edges] = generate_Lshaped(ref);
    elseif geom == 3
        edge = [];
        load borrar.mat
        vertices = coord;
        elements = elem_back;
        edges = [];
        edges = edge;
    else
        error('Domain not-defined')
    end

    n_edge = size(edges,1); % number of edges
%     %------- plotting the mesh ----------% 
%     figure
%     draw_mesh(elements(:,2:4),vertices(:,2:3))
%     %------------------------------------

    n_elem = size(elements,1); % number of elements
    n_vert = size(vertices,1); % number of vertices

    % Storing the elements in a structure
    for i = 1 : n_elem

    %       storing the coordinates of the vertices 
    %       ---------------------------------------
            elem(i).vert_gn(1) = elements(i,2); % local vertex 1 global number
            elem(i).vert_gn(2) = elements(i,3); % local vertex 2 global number
            elem(i).vert_gn(3) = elements(i,4); % local vertex 3 global number
    %
    %       verifying that the elements are not `flat'
    %       and that they are correctly oriented
    %       ------------------------------------------
            i1=elem(i).vert_gn(1); % local vertex 1 global number
            i2=elem(i).vert_gn(2); % local vertex 2 global number
            i3=elem(i).vert_gn(3); % local vertex 3 global number

    %       computing the area of the triangle
            area = area_of_triangle(vertices(i1,2:3),vertices(i2,2:3),vertices(i3,2:3));

    %       checking
            if (area == 0)  %flat element
               error('Element has zero area. This mesh cannot be used!!')

            elseif (area < 0)  %incorrect orientation
               disp('Element has the wrong orientation. This has been fixed')

    %          correcting the original orientation of the triangle
    %          interchanging the local vertices `2' and `3'

                elem(i).vert_gn(2) = i3;
                elem(i).vert_gn(3) = i2;

                i2 = elem(i).vert_gn(2); %setting the new i2
                i3 = elem(i).vert_gn(3); %setting the new i3

                area = -area; %now the area is positive
            end

    %       computation of the vertices
            elem(i).vert(1,1:2) = vertices(i1,2:3);
            elem(i).vert(2,1:2) = vertices(i2,2:3);
            elem(i).vert(3,1:2) = vertices(i3,2:3);
            
    %       computation of the matrix and vector of affine mapping x = BT*x_hat+bT
            elem(i).BT = [vertices(i1,2)-vertices(i3,2) vertices(i2,2)-vertices(i3,2);...
                          vertices(i1,3)-vertices(i3,3) vertices(i2,3)-vertices(i3,3)];
            elem(i).bT = [vertices(i3,2) vertices(i3,3)]';

    %       computation of the normals
            elem(i).norm(1,1:2) = [-(vertices(i2,3)-vertices(i3,3));  vertices(i2,2)-vertices(i3,2)];
            elem(i).norm(2,1:2) = [-(vertices(i3,3)-vertices(i1,3));  vertices(i3,2)-vertices(i1,2)];
            elem(i).norm(3,1:2) = [-(vertices(i1,3)-vertices(i2,3));  vertices(i1,2)-vertices(i2,2)];

    %       computation of the length of the edges
            elem(i).edge_l(1) = norm(elem(i).norm(1,:));
            elem(i).edge_l(2) = norm(elem(i).norm(2,:));
            elem(i).edge_l(3) = norm(elem(i).norm(3,:));

    %       ... the area
            elem(i).area=area;
            

    end %loop on the elements

    % Computation of the edges (faces)


    %computation of the orientation of the edges and the auxiliary variable "isize"
    %-----------------------------------------------------------------------------
    isize = zeros(n_vert,1); % number of elements sharing a vertex
    for i = 1 : n_elem % loop on the elements
        i1=elem(i).vert_gn(1);
        i2=elem(i).vert_gn(2);
        i3=elem(i).vert_gn(3);

        isize(i1) = isize(i1) + 1;
        isize(i2) = isize(i2) + 1;
        isize(i3) = isize(i3) + 1;

        elem(i).edge_or(1) = 1;
        elem(i).edge_or(2) = 1;
        elem(i).edge_or(3) = 1;

        if(i2>i3) 
            elem(i).edge_or(1)=-1;
        end
        if(i3>i1) 
            elem(i).edge_or(2)=-1;
        end
        if(i1>i2) 
            elem(i).edge_or(3)=-1;
        end

    end %loop on the elements

    % computation of the auxiliary variable "vert2elem". This rountine 
    % finds, for each vertex, the elements where it belongs
    
    icount = zeros(n_vert,1);
    for i = 1 : n_elem % loop on the elements

        i1 = elem(i).vert_gn(1);
        i2 = elem(i).vert_gn(2);
        i3 = elem(i).vert_gn(3);

        icount(i1) = icount(i1) + 1;                        
        icount(i2) = icount(i2) + 1;
        icount(i3) = icount(i3) + 1;

        vert2elem(i1).elem_gn(icount(i1)) = i;
        vert2elem(i2).elem_gn(icount(i2)) = i;
        vert2elem(i3).elem_gn(icount(i3)) = i;

        vert2elem(i1).vert_ln(icount(i1)) = 1;
        vert2elem(i2).vert_ln(icount(i2)) = 2;
        vert2elem(i3).vert_ln(icount(i3)) = 3;

    end

    % computation of the global numbers of the edges
    % ----------------------------------------------
    % and of the list of boundary edges
    % ---------------------------------

    % setting the number of boundary edges
      n_bedg = 2*n_vert-n_elem; % This is an upper bound to be updated pre
    
    % initializing the auxiliary counters
      iedg_count = n_bedg; % interior edge counter
      bedg_count = 0;      % boundary edge counter
      edges = [];
      perm = [1 2 3 1 2 ];
      for i = 1 : n_elem % loop on the elements
          for  l = 1 :3 % loop on the local sides
               i2 = elem(i).vert_gn(perm(l+1));
               i3 = elem(i).vert_gn(perm(l+2));
    
    %           calculation of the number of the neighboring elements "nghb"
    %           and the local number of the corresponding edge "ls"
    
    %           initialization
    %           if these numbers do not change, the edge is a boundary edge
    %           otherwise it is an interior edge and nghb, ls > 0
                nghb = 0;
                ls = 0;
    %
    %           computation of the actual values of "nghb" and "ls"
                for j = 1 : isize(i2) % loop on the number of triangles sharing the vertex "i2"
                   nj = vert2elem(i2).elem_gn(j);
                   for k = 1 : isize(i3) % loop on the number of triangles sharing the vertex "i3"
                       nk = vert2elem(i3).elem_gn(k);
                       if (nj==nk) & (nj~=i) 
                           nghb = nj;
                           ls=perm(vert2elem(i2).vert_ln(j)+1);
                       end
                   end
                end 
    
    %           calculation of the global number of the edge
                if nghb==0 % this is a boundary edge
                   bedg_count = bedg_count + 1;
                   elem(i).edge_gn(l) = bedg_count;
    
                   bedg(bedg_count).vert(1,:) = elem(i).vert(perm(l+1),:);
                   bedg(bedg_count).vert(2,:) = elem(i).vert(perm(l+2),:);
                   bedg(bedg_count).edge_l = norm(bedg(bedg_count).vert(2)-bedg(bedg_count).vert(1));
                   bedg(bedg_count).edge_gn = bedg_count;
                   bedg(bedg_count).edge_or = elem(i).edge_or(l);
                   bedg(bedg_count).vert_gn(1) = i2;
                   bedg(bedg_count).vert_gn(2) = i3;
                   bedg(bedg_count).elem_gn = i;
                
                   edges(bedg_count).vert(1,:) = elem(i).vert(perm(l+1),:);
                   edges(bedg_count).vert(2,:) = elem(i).vert(perm(l+2),:);
                   edges(bedg_count).edge_l = norm(bedg(bedg_count).vert(2,:)-bedg(bedg_count).vert(1,:));
                   edges(bedg_count).edge_gn = bedg_count;
                   edges(bedg_count).edge_or = elem(i).edge_or(l);
                   edges(bedg_count).vert_gn(1) = i2;
                   edges(bedg_count).vert_gn(2) = i3;
                   edges(bedg_count).elem_gn = i; 
                   hmax = max(hmax,edges(bedg_count).edge_l);
                elseif (elem(i).edge_or(l)==1) % the edge has positive orientation
                   iedg_count = iedg_count + 1;
                   
                   edges(iedg_count).vert(1,:) = elem(i).vert(perm(l+1),:);
                   edges(iedg_count).vert(2,:) = elem(i).vert(perm(l+2),:);
                   edges(iedg_count).edge_l = norm(edges(iedg_count).vert(2)-edges(iedg_count).vert(1));
                   edges(iedg_count).edge_gn = iedg_count;
                   edges(iedg_count).edge_or = elem(i).edge_or(l);
                   edges(iedg_count).vert_gn(1) = i2;
                   edges(iedg_count).vert_gn(2) = i3;
                   edges(iedg_count).elem_gn = i;
                   
                   elem(i).edge_gn(l) = iedg_count;
                   elem(nghb).edge_gn(ls) = iedg_count; % taking care of the neighboring element
                   hmax = max(hmax,edges(iedg_count).edge_l);
                end
                
          end     % loop on the edges
%     end % loop on the elements
    %ind = input('Refinar (1: Si, 0: No)--> '); 
    %ref = ref + 1;
      end


n_edge = iedg_count - n_bedg;
n_bedg0 = n_bedg;
n_bedg = bedg_count;
n_edge = n_edge + n_bedg;

aa(1:n_bedg)=edges(1:n_bedg);
aa(n_bedg+1:n_edge) = edges(n_bedg0+1:end);
edges = aa;
% Update number of edges
for i = 1 : n_elem % loop on the elements
    for  l = 1 :3 % loop on the local sides
        if elem(i).edge_gn(l) >  n_bedg
            elem(i).edge_gn(l) = elem(i).edge_gn(l) - n_bedg0 + n_bedg;
        end
    end
end






