function c = tensor_C(p,ex)

if ex == 1
    c = ones(length(p(1,:)),1);
elseif ex == 2
    c = ones(length(p(1,:)),1);
elseif ex == 3
    c = ones(length(p(1,:)),1);  
elseif ex == 4
    c = ones(length(p(1,:)),1);   
else
    error('Permeability tensor not defined')
end