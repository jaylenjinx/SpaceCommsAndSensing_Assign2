function q = qfunc_local(x)
% QFUNC_LOCAL  Gaussian tail probability Q(x) = 0.5*erfc(x/sqrt(2)).
%   Defined locally so the code runs without the Communications Toolbox.
    q = 0.5 * erfc(x / sqrt(2));
end
