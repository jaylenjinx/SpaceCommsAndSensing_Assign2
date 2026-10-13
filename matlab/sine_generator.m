function i_n = sine_generator(L, A, w0, n0)
% SINE_GENERATOR  Sinusoidal interference i(n) = A sin(w0 n) (Project 12.8).
%   i_n = sine_generator(L, A, w0) returns L samples starting at n = 0.
%   i_n = sine_generator(L, A, w0, n0) starts at n = n0, so a long run can
%   be generated in blocks without breaking the phase of the jammer.
    if nargin < 4, n0 = 0; end
    n   = n0 + (0:L-1);
    i_n = A * sin(w0 * n);
end
