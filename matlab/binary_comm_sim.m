function [ber, nerr, nbits] = binary_comm_sim(M, a, opts)
% BINARY_COMM_SIM  Monte-Carlo run of the systems in Figs 12.16 and 12.18.
%   [ber, nerr, nbits] = binary_comm_sim(M, a, opts)
%     M     samples per bit
%     a     uniform noise half-width (noise on (-a,a), Pn = a^2/3)
%     opts  struct with optional fields
%       .use_pn    true to spread with a +/-1 PN sequence   (default false)
%       .A         sinusoidal interference amplitude         (default 0)
%       .w0        interference frequency in rad/sample      (default 0)
%       .min_err   stop once this many errors are seen       (default 200)
%       .max_bits  hard limit on simulated bits              (default 1e6)
%       .block     bits simulated per block                  (default 2e4)
%   Bits are simulated in blocks so large M does not exhaust memory; the
%   jammer phase runs continuously across blocks.
    if nargin < 3, opts = struct(); end
    use_pn   = getopt(opts, 'use_pn', false);
    A        = getopt(opts, 'A', 0);
    w0       = getopt(opts, 'w0', 0);
    min_err  = getopt(opts, 'min_err', 200);
    max_bits = getopt(opts, 'max_bits', 1e6);
    block    = getopt(opts, 'block', 2e4);
    % keep each block below ~2e7 samples
    block = max(1, min(block, floor(2e7 / M)));

    nerr = 0; nbits = 0; n0 = 0;
    while nerr < min_err && nbits < max_bits
        Nb = min(block, max_bits - nbits);
        L  = Nb * M;
        b  = bit_generator(Nb);                 % binary data generator
        s  = modulator(b, M);                   % modulator
        if use_pn
            p = pn_generator(L);                % PN sequence generator
            s = s .* p;                         % spreading
        end
        r = s + noise_generator(L, a);          % channel: wideband noise
        if A ~= 0
            r = r + sine_generator(L, A, w0, n0);   % channel: sine jammer
        end
        if use_pn
            r = r .* p;                         % despreading
        end
        y = demodulator(r, M);                  % demodulator
        [~, e] = detector(y, b);                % detector / error counter
        nerr = nerr + e; nbits = nbits + Nb; n0 = n0 + L;
    end
    ber = nerr / nbits;
end

function v = getopt(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end
