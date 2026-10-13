function pe = exact_ber_uniform(M, a)
% EXACT_BER_UNIFORM  Exact bit error probability for Project 12.7.
%   pe = exact_ber_uniform(M, a) returns P(W > M), where W is the sum of M
%   i.i.d. Uniform(-a,a) noise samples. By symmetry this is the bit error
%   probability of the integrate-and-dump receiver for either bit value.
%   The pdf of W is obtained by numerically convolving M uniform pdfs.
%   Works for scalar M and vector a.
    pe = zeros(size(a));
    for k = 1:numel(a)
        ak = a(k);
        if ak <= 1                       % |W| <= M*a <= M: no error possible
            pe(k) = 0; continue;
        end
        dx = ak / 400;                   % grid spacing
        x1 = -ak:dx:ak;
        f1 = ones(size(x1)) / (2*ak);    % Uniform(-a,a) pdf
        f  = f1;
        for m = 2:M
            f = conv(f, f1) * dx;
        end
        x  = (0:numel(f)-1) * dx - M*ak; % support of the M-fold sum
        f  = f / (sum(f) * dx);          % remove discretisation drift
        pe(k) = sum(f(x > M)) * dx;
    end
end
