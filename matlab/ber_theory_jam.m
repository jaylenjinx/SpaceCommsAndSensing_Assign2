function [ber_nopn, ber_pn] = ber_theory_jam(M, A, w0, Pn)
% BER_THEORY_JAM  Analytical BER with sinusoidal jamming (Project 12.8).
%   [ber_nopn, ber_pn] = ber_theory_jam(M, A, w0, Pn) averages over the
%   jammer phase phi seen at the start of a bit (uniform on [0, 2pi)).
%   The summed noise is treated as Gaussian with variance M*Pn (CLT).
%
%   Without PN the jammer adds a deterministic offset to the bit sum,
%       J(phi) = A * sum_{n=0}^{M-1} sin(w0 n + phi),
%   so  Pe(phi) = [Q((M+J)/s) + Q((M-J)/s)] / 2,  s = sqrt(M Pn).
%
%   With PN the jammer is multiplied by i.i.d. +/-1 chips, so its sum is
%   zero mean with variance A^2 * sum sin^2(w0 n + phi), giving
%       Pe(phi) = Q( M / sqrt(M Pn + A^2 sum sin^2(w0 n + phi)) ).
%   Scalar M, w0, Pn; A may be a vector.
    phi = 2*pi*(0:359)/360;
    n   = (0:M-1).';
    S   = sin(w0*n + phi);                 % M x 360
    J1  = sum(S, 1);                       % sum sin      per phase
    J2  = sum(S.^2, 1);                    % sum sin^2    per phase
    sig = sqrt(M*Pn);
    ber_nopn = zeros(size(A)); ber_pn = ber_nopn;
    for k = 1:numel(A)
        J = A(k) * J1;
        ber_nopn(k) = mean(0.5*(qfunc_local((M + J)/sig) + qfunc_local((M - J)/sig)));
        ber_pn(k)   = mean(qfunc_local(M ./ sqrt(M*Pn + A(k)^2 * J2)));
    end
end
