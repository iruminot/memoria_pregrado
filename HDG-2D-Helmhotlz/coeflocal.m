function coefloc = coeflocal(Au_store, Ap_store, Bu_store, Bp_store, gamma,elem,k)
d = k+1;
r = 1:d;
for i = 1:size(elem,2)
    Kalpha = Au_store{i};
    Kbeta = Ap_store{i};
    btildealpha = Bu_store{i};
    btildebeta = Bp_store{i};
    e_gn= elem(i).edge_gn;
    coefloc(i).gammaF(:,1) = gamma(r+(e_gn(1)-1)*d);
    coefloc(i).gammaF(:,2) = gamma(r+(e_gn(2)-1)*d);
    coefloc(i).gammaF(:,3) = gamma(r+(e_gn(3)-1)*d);
    coefloc(i).alpha=Kalpha*[coefloc(i).gammaF(:,1);coefloc(i).gammaF(:,2);coefloc(i).gammaF(:,3)]+btildealpha;
    coefloc(i).beta=Kbeta*[coefloc(i).gammaF(:,1);coefloc(i).gammaF(:,2);coefloc(i).gammaF(:,3)]+btildebeta;
end
    
    
    
    