function  [coord elem edge] = generate_Lshape(a,areaT,m)
%------------------------------------------------------------------------------
%Description: Generates the mesh for the L-shape [-a,a]^2 \ ]0,a[^2
%Programmer: Sergio Caucao
%Date: October 5, 2014
%------------------------------------------------------------------------------

if m == 1
    coords = [1 -a -a;
              2 a -a;
              3 a 0;
              4 0 0;
              5 0 a;
              6 -a a];
      
    conn = [1 1 2 1;
            2 2 3 1;
            3 3 4 2;
            4 4 5 2;
            5 5 6 2;
            6 6 1 1];

    sc = size(coords,1);
    se = size(conn,1);
%%%%%%%%%%%  saving data on a file to be read by triangle %%%%%%%%%%% 
    file = fopen('Lshape.poly','w');
    fprintf(file, '%4i', [sc 2 0 1]);
    fprintf(file, '\n');
    fprintf(file, '%4i %2.16f %12.16f\n', coords');
    fprintf(file, '\n');
    fprintf(file, '%4i', [se 1 1]);
    fprintf(file, '\n');
    fprintf(file, '%4i %4i %4i %4i\n', conn');
    fprintf(file, '\n');
    fprintf(file, '%4i\n',0);
    fprintf(file, '\n');
    fprintf(file, '%4i',1);
    fprintf(file, '\n');
    fprintf(file, '%4i', [1 0 0]);
    fclose(file);

%%%%%%%%%%%  calling triangle %%%%%%%%%%% 
    aux=['./triangle -pqea',num2str(areaT),' Lshape.poly'];	
    system(aux)

    filename = ['Lshape.1'];
    [coord elem edge] = geometry(filename);
else
    for j = (m-1):m
        filename = ['Lshape.',num2str(j)];
    	[coord_m elem_m edge_m] = geometry(filename);
    	coord_m(:,1) = [];
    	elem_m(:,1)  = [];
	    n_elem = size(elem_m,1);
        
        ax = coord_m(elem_m(:,1),1) - coord_m(elem_m(:,3),1);  % x(2)-x(1)
        ay = coord_m(elem_m(:,1),2) - coord_m(elem_m(:,3),2);  % y(2)-y(1)
        bx = coord_m(elem_m(:,2),1) - coord_m(elem_m(:,3),1);  % x(3)-x(1) 
        by = coord_m(elem_m(:,2),2) - coord_m(elem_m(:,3),2);  % y(3)-y(1)
        detB = ax.*by - ay.*bx;
	    
	    y1 = 1:n_elem;
	    y2 = abs(detB)'/4;
	    y  = [y1;y2];

        filearea = ['Lshape.',num2str(j),'.area'];	
	    fid      = fopen(filearea,'wt');
	    fprintf(fid, '%6d\n', n_elem);
	    fprintf(fid, '%6d %6.14f\n', y);
	    fclose(fid);

        aux = ['./triangle -raeq27ApQ ',' Lshape.',num2str(j)];
	    unix(aux);	    
    end
    aux = ['rm Lshape.',num2str(m),'*'];
    unix(aux);
    aux = ['mv Lshape.',num2str(m+1),'.node Lshape.',num2str(m),'.node'];
    unix(aux);
    aux = ['mv Lshape.',num2str(m+1),'.ele Lshape.',num2str(m),'.ele'];
    unix(aux);
    aux = ['mv Lshape.',num2str(m+1),'.edge Lshape.',num2str(m),'.edge'];
    unix(aux);
    aux = ['mv Lshape.',num2str(m+1),'.poly Lshape.',num2str(m),'.poly'];
    unix(aux);
    filename = ['Lshape.',num2str(m)];
    [coord elem edge] = geometry(filename);
end

%%%%%%%%%%%  reading mesh from triangle %%%%%%%%%%% 

function  [coord_m elem_m edge_m] = geometry(filename)

%%% reading coordinates %%%
file  = fopen([filename,'.node']);
num   = fscanf(file ,'%i',4); nverf=num(1);
coord_m = fscanf(file ,'%f',[4,nverf]);
coord_m = coord_m';
coord_m(:,4) = [];
fclose(file);

%%% reading elements %%%
file  = fopen([filename,'.ele']);
num  = fscanf(file,'%i',3); nelf=num(1);
elem_m = fscanf(file,'%f',[4,nelf]);
elem_m = elem_m';
fclose(file);


%%% reading edges %%%
file  = fopen([filename,'.edge']);
num  = fscanf(file,'%i',2); nelf=num(1);
edge_m = fscanf(file,'%i',[4,nelf]);
edge_m = edge_m';
fclose(file);
