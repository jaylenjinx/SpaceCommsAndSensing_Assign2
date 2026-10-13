function y = demodulator(r, M)
% DEMODULATOR  Integrate-and-dump demodulator (Project 12.7, block 4).
%   y = demodulator(r, M) sums each block of M successive received
%   samples. The receiver is assumed to be bit-synchronised, so block k
%   is r((k-1)*M+1 : k*M). Returns one decision statistic per bit.
    y = sum(reshape(r, M, []), 1);
end
