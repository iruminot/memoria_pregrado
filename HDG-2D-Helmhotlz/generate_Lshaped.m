function  [coord,elem,edge] = generate_Lshaped(ref)
%-------------------------------------------------------------------------
%Description: Generates the mesh for the rectangle [0,a]x[0,b]
%Programmer: Manuel Solano
%Date: Sept 6, 2013
%-------------------------------------------------------------------------

% coords = [-1 -1; 1 -1; 1 0];
% aux_x = [linspace(0.1,0,2^ref)' zeros(2^ref,1) ];
% aux_y = [zeros(2^ref,1) linspace(0.1,0,2^ref)' ];
% aux_y(end,:) = [];
% coords = [coords;aux_x;aux_y;[0 1;-1 1]];
% 
% sc =size(coords,1);
% coords = [[1:sc]' coords];
% 
% conn = [[1:sc-1]' [2:sc]'];
% conn = [conn;[sc 1]];
% sc = size(conn,1);
% conn = [[1:sc]' conn ones(sc,1)];



coords = [1 -1 -1;
          2 1 -1;
          3 1 0;
          4 0 0;
          5 0 1;
          6 -1 1];
      
conn = [1 1 2 1;
        2 2 3 1;
        3 3 4 1;
        4 4 5 1;
        5 5 6 1;
        6 6 1 1];
% 

    sc = size(coords,1);
    se = size(conn,1);
    %%%%%%%%%%%  saving data on a file to be read by triangle %%%%%%%%%%% 
    %file = fopen(strcat(name,'.poly'),'w');
    file = fopen('aux_file.poly','w');
    fprintf(file, '%4.1i', [sc 2 0 0]);
    fprintf(file, '\n');
    fprintf(file, '%4.1i %2.16f %12.16f\n', coords');
    fprintf(file, '\n');
    fprintf(file, '%4.1i', [se 1 ]);
    fprintf(file, '\n');
    fprintf(file, '%4.1i %4.1i %4.1i %4.1i\n', conn');
    fprintf(file, '\n');
    fprintf(file, '%4.1i',0);
    fprintf(file, '\n');
    fclose(file);

    area = (0.25)^(ref);
    system(['./triangle -pqea',num2str(area,'%1.19f'), ' aux_file.poly']);
    %system(['triangle -pqea',num2str(area,'%1.19f'), ' aux_file.poly']);
%     
% if ref == 0
%     %%%%%%%%%%%  calling triangle %%%%%%%%%%% 
%     !triangle -pqea0.25 aux_file.poly
%     %!./triangle -pqea0.125 aux_file.poly
%     %!triangle -pqea0.25 aux_file.poly
% 
% elseif ref == 1
%        %%%%%%%%%%%  calling triangle %%%%%%%%%%% 
%     !triangle -pqea0.0625 aux_file.poly 
% elseif ref == 2
%        %%%%%%%%%%%  calling triangle %%%%%%%%%%% 
%     !triangle -pqea0.015625 aux_file.poly
% else
%     !triangle -pqea0.00390625 aux_file.poly
% end


%%%%%%%%%%%  reading mesh from triangle %%%%%%%%%%% 

%%% reading coordinates %%%
file = fopen('aux_file.1.node');
num=fscanf(file ,'%i',4); nverf=num(1);
coord=fscanf(file ,'%f',[4,nverf]);
coord = coord';
coord(:,4)=[];
fclose(file);

%%% reading elements %%%
file = fopen('aux_file.1.ele');
num=fscanf(file,'%i',3); nelf=num(1);
elem=fscanf(file,'%i',[4,nelf]);
elem = elem';
fclose(file);

%%% reading edges %%%
file = fopen('aux_file.1.edge');
num=fscanf(file,'%i',2); nelf=num(1);
edge=fscanf(file,'%i',[4,nelf]);
edge = edge';
fclose(file);