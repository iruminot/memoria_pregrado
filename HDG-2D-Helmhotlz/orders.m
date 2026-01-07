global  ex geom



geom = 1; %Rectangle
% geom = 2; %L-shaped
% geom = 3; %ring for testing
ex =  input('Choose example: ');
k = input('Choose polinomial degree: '); 
param = [1:5];

error_q = [];
error_u = [];
tic
for j = 1 : length(param)
    tic
    [h(j) N(j) error_u(j) error_q(j) error_lam(j)] = main(k,param(j),geom,ex);
    toc
end
tic

r_q = 2.*log(error_q(2:end)./error_q(1:end-1))./log(N(1:end-1)./N(2:end))
r_u = 2.*log(error_u(2:end)./error_u(1:end-1))./log(N(1:end-1)./N(2:end))
r_lam = 2.*log(error_lam(2:end)./error_lam(1:end-1))./log(N(1:end-1)./N(