function  [coord,elem,edge] = generate_rectangle(a,b,ref)
%-------------------------------------------------------------------------
%Description: Generates the mesh for the rectangle [0,a]x[0,b]
%Programmer: Manuel Solano
%Date: Sept 6, 2013
%-------------------------------------------------------------------------
    coords = [1 0 0;
              2 a 0;
              3 a b;
              4 0 b];

    conn = [1 1 2 1;
            2 2 3 1;
            3 3 4 2;
            4 4 1 2];

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