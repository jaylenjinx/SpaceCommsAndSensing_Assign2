function p = pn_generator(L)
% PN_GENERATOR  Binary (+/-1) pseudo-noise sequence (Project 12.8).
%   p = pn_generator(L) returns L equiprobable +/-1 chips generated with a
%   random number generator. The receiver uses the identical sequence.
    p = 2*double(rand(1, L) >= 0.5) - 1;
end
