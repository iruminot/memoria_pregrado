function u = q_exact(x,y,ex)


if ex == 1    
    u = -[cos(x).*sin(y);sin(x).*cos(y)];
elseif ex == 2   
    u = -[1+0.*x;0.*y];
elseif ex == 3;
    a = 0.05;
    b = a;
    u1 = -2.*(a-x)./((a - x).^2 + (b - y).^2).^2;
    u2 = -2.*(b-y)./((a - x).^2 + (b - y).^2).^2;
    u = [u1; u2];
elseif ex == 4;
    u1 = -2.*x;
    u2 = 0.*x;
    u = [u1; u2];
else
    error('Example not defined')
end

