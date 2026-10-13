function b = bit_generator(N)
% BIT_GENERATOR  Binary data generator (Project 12.7, block 1).
%   b = bit_generator(N) returns a 1xN row vector of independent,
%   equiprobable binary digits (0 or 1).
    b = double(rand(1, N) >= 0.5);
end
