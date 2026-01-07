function area = area_of_triangle(v1,v2,v3)
%-------------------------------------------------------------------
% Description: Computes the area of the triangle determined by
% ***********  the points v1,v2,v3. 
%                 
% Note that, if the points are oriented in a negative way, the 
% area is nonpositive.
%-------------------------------------------------------------------

ax = v2(1)-v1(1);
ay = v2(2)-v1(2);

bx = v3(1)-v1(1);
by = v3(2)-v1(2);

area = (ax*by-ay*bx)/2;