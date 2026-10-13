function w = noise_generator(L, a)
% NOISE_GENERATOR  Uniform additive noise (Project 12.7, block 3).
%   w = noise_generator(L, a) returns L i.i.d. samples uniformly
%   distributed on (-a, a). The noise power is Pn = a^2/3.
    w = a * (2*rand(1, L) - 1);
end
