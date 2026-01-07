function z = f(x,y,ex)


if ex == 1
    z = 2.*sin(x).*sin(y);
elseif ex == 2
    z = 0.*x;
elseif ex == 3
    a = 0.05;
    b = a;
    z = -4./((a - x).^2 + (b - y).^2).^2;
elseif ex == 4
    z = -2+0.*x;
else
    error('Example not defined')
end

