function p = u_exact(x,y,ex)


if ex == 1
    p = sin(x).*sin(y);
elseif ex == 2    
    p = x;
elseif ex == 3
    a = 0.05;
    b = a;
    p = 1./((x-a).^2+(y-b).^2);  
elseif ex == 4    
    p = x.^2; 
else
    error('Example not defined')
end

