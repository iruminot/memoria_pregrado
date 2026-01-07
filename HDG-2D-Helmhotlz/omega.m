function val = omega(ex)
% Frequency omega - \nabla u - omega^2 u = f


if ex == 1
    val = 0;
elseif ex == 2
    val = 0;
else
    error('coefficient alpha not defined for this problem');
end