function jp = jacobi_poly(alpha,beta,n,x)
%-----------------------------------------------------------------
%     Description: This function returns the Jacobi polynomial
%     ***********  P^{alpha, beta}_{n} evaluated at x.
%-----------------------------------------------------------------

pmm = ones(length(x),1);
if n==0
    jp = pmm;
	return
else
    pmmp1 = 0.5.*(alpha-beta)+0.5.*(alpha+beta+2).*x;
    if n == 1
		jp = pmmp1;
		return
    else
        for i = 2 : n
			c0 = -2.*(i-1.+alpha).*(i-1.+beta).*(2.*(i-1)+alpha+beta+2);
      		c2 = 2.*(i-1.+1.).*(i-1.+alpha+beta+1.).*(2.*(i-1.)+alpha+beta);
            c11=(2.*(i-1.)+alpha+beta+1.).*(alpha^2-beta.^2);
      		c12=(2.*(i-1.)+alpha+beta+2.).*(2.*(i-1.)+alpha+beta+1.).*(2.*(i-1.)+alpha+beta);
			pnn=((c11+c12.*x).*pmmp1+c0.*pmm)./c2;
			pmm = pmmp1;
			pmmp1 = pnn;
        end
		jp = pnn;
		return
    end
end
