function draw_mesh(element,coordinate)
%-------------------------------------------------------------------------
% Description: Draws a mesh given the elements and coordinates
% Programmer: Manuel Solano
% Date: Sept 6, 2013
%-------------------------------------------------------------------------

hold on
for j=1:size(element,1)
    trisurf([1 2 3],coordinate(element(j,:),1),coordinate(element(j,:),2),ones(3,1),'edgealpha',1, 'facecolor','w' );
end
view(0,90);
hold off
