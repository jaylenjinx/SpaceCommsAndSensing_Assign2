function s = modulator(b, M)
% MODULATOR  Rectangular-pulse binary modulator (Project 12.7, block 2).
%   s = modulator(b, M) maps every bit 1 to M consecutive +1's and every
%   bit 0 to M consecutive -1's, i.e. a sampled antipodal rectangular pulse.
%   Output length is M*length(b).
    a = 2*b(:).' - 1;                 % 0 -> -1, 1 -> +1
    s = reshape(repmat(a, M, 1), 1, []);
end
