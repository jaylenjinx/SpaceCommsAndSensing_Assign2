function pe = exact_ber_uniform(M, a)
% EXACT_BER_UNIFORM  Exact bit error probability for Project 12.7.
%   pe = exact_ber_uniform(M, a) returns P(W > M), where W is the sum of M
%   i.i.d. Uniform(-a,a) noise samples. By symmetry this is the bit error
%   probability of the integrate-and-dump receiver for either bit value.
%
%   W = 2a*S - M*a with S ~ Irwin-Hall(M) (sum of M Uniform(0,1)), so
%       P(W > M) = P(S > x) = P(S < M - x),   x = M(1 + a)/(2a),
%   and the Irwin-Hall CDF  F(y) = (1/M!) sum_{k=0}^{floor(y)} (-1)^k
%   C(M,k) (y-k)^M  is evaluated on the short (y <= M/2) side, which keeps
%   the alternating sum well conditioned in double precision for M <= 30.
%   Works for scalar M and vector a.
    if M > 30
        error('exact_ber_uniform:M', 'Closed form is only well conditioned for M <= 30.');
    end
    pe = zeros(size(a));
    for i = 1:numel(a)
        y = M - M*(1 + a(i)) / (2*a(i));   % = M(a-1)/(2a)
        if y <= 0                           % a <= 1: |W| <= M, no errors
            pe(i) = 0; continue;
        end
        k = 0:floor(y);
        logt = gammaln(M+1) - gammaln(k+1) - gammaln(M-k+1) ...  % log C(M,k)
             + M*log(y - k) - gammaln(M+1);                       % (y-k)^M / M!
        logt(y - k <= 0) = -Inf;
        pe(i) = sum((-1).^k .* exp(logt));
    end
end
