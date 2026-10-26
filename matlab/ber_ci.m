function [lo, hi] = ber_ci(k, n, conf)
% BER_CI  Clopper-Pearson (exact binomial) confidence interval for a BER.
%   [lo, hi] = ber_ci(k, n) returns the 95% interval for the error
%   probability after k errors in n bits; [lo, hi] = ber_ci(k, n, conf)
%   uses confidence level conf. Works element-wise on vectors.
    if nargin < 3, conf = 0.95; end
    alpha = 1 - conf;
    lo = zeros(size(k)); hi = ones(size(k));
    for i = 1:numel(k)
        if k(i) > 0,    lo(i) = betaincinv(alpha/2, k(i), n(i) - k(i) + 1); end
        if k(i) < n(i), hi(i) = betaincinv(1 - alpha/2, k(i) + 1, n(i) - k(i)); end
    end
end
