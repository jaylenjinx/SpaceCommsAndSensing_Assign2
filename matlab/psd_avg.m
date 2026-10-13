function [P, w] = psd_avg(x, Nfft)
% PSD_AVG  Averaged, Hann-windowed periodogram (no Signal Toolbox needed).
%   [P, w] = psd_avg(x, Nfft) splits x into non-overlapping segments of
%   length Nfft, applies a Hann window, averages |FFT|^2 normalised so a
%   white sequence of power s2 gives P = s2, and returns the one-sided
%   estimate on w in [0, pi] rad/sample.
    K = floor(numel(x) / Nfft);
    X = reshape(x(1:K*Nfft), Nfft, K);
    h = 0.5 - 0.5*cos(2*pi*(0:Nfft-1).'/Nfft);   % Hann window
    X = X .* repmat(h, 1, K);
    P = mean(abs(fft(X, Nfft, 1)).^2, 2).' / sum(h.^2);
    P = P(1:Nfft/2+1);
    w = 2*pi*(0:Nfft/2) / Nfft;
end
