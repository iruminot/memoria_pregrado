function leg = legendre(n,x)
%-----------------------------------------------------------------
%     Description: This function returns the Legendre polynomial
%     ***********  L_{n}(x).
%-----------------------------------------------------------------

pmm = 1;
if n == 0
    leg = pmm;
    return
else
    pmmp1 = x*pmm;
    if n == 1
        leg = pmmp1;
		return
    else
        for i = 2 : n
			pnn = (x.*(2.*i-1.).*pmmp1-(i-1).*pmm)./i;
			pmm = pmmp1;
			pmmp1 = pnn;
        end
		leg = pnn;
		return
    end
end