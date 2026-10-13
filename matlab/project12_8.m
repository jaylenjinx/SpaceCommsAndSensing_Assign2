%% Project 12.8: Binary spread-spectrum communications (Ingle & Proakis 4e)
% Simulates the system of Fig. 12.18: the Project 12.7 system with the
% modulator output multiplied by a +/-1 PN sequence p(n), a channel that
% adds wideband noise w(n) and a sinusoidal jammer i(n) = A sin(w0 n),
% and a receiver that multiplies by the same p(n) before demodulating.
%
% Experiments
%   1. Brief: BER with and without PN for M = 50, 100, 500, 1000, A = M,
%      for an in-band (low-frequency) and an out-of-band jammer.
%   2. BER against jammer amplitude A (jammer-to-signal ratio) per M.
%   3. BER against jammer frequency w0 for fixed M.
%   4. Power spectra before and after despreading.
clear; close all; clc;
rng(128);
if exist('OCTAVE_VERSION', 'builtin'), graphics_toolkit('gnuplot'); end

a  = sqrt(3);  Pn = a^2/3;                 % wideband noise, Pn = 1
M_list = [50 100 500 1000];
w_in   = 0.002;                            % rad/sample, inside 2*pi/M for all M
w_out  = 1;                                % rad/sample (~0.318*pi)
min_err = Inf;                             % experiment 1: fixed bit budget
bits_for = @(M) max(2e4, round(2e7/M));    % ~2e7 samples per point

%% 1. Brief experiment: A = M
fprintf('Experiment 1: A = M, Pn = %.2f\n', Pn);
w_set = [w_in w_out];
T = zeros(numel(M_list), 2*numel(w_set));  % simulated BERs
Tth = T;                                   % theory
for i = 1:numel(M_list)
    M = M_list(i); A = M;
    for j = 1:numel(w_set)
        o = struct('A', A, 'w0', w_set(j), 'min_err', min_err, 'max_bits', bits_for(M));
        o.use_pn = false; ber0 = binary_comm_sim(M, a, o);
        o.use_pn = true;  ber1 = binary_comm_sim(M, a, o);
        [t0, t1] = ber_theory_jam(M, A, w_set(j), Pn);
        T(i, 2*j-1:2*j)   = [ber0 ber1];
        Tth(i, 2*j-1:2*j) = [t0 t1];
        fprintf('  M = %4d  w0 = %.3f  no PN: %.4f (th %.4f)   PN: %.4f (th %.4f)\n', ...
            M, w_set(j), ber0, t0, ber1, t1);
    end
end

f = figure('visible', 'off');
semilogx(M_list, T(:,1), 'ro-', M_list, T(:,2), 'bo-', ...
         M_list, T(:,3), 'rs--', M_list, T(:,4), 'bs--', 'LineWidth', 1.2, 'MarkerSize', 7);
grid on; ylim([0 0.6]); xlim([40 1200]);
xlabel('Samples per bit M (with A = M)'); ylabel('Bit error rate');
title('Project 12.8: BER with jammer amplitude A = M, P_n = 1');
legend(sprintf('No PN, \\omega_0 = %.3f', w_in), sprintf('PN, \\omega_0 = %.3f', w_in), ...
       sprintf('No PN, \\omega_0 = %g', w_out), sprintf('PN, \\omega_0 = %g', w_out), ...
       'Location', 'southwest');
save_fig(f, 'p128_ber_A_eq_M');

%% 2. BER against jammer amplitude (in-band jammer)
A_list = logspace(-1, 3, 25);
jsr_db = 10*log10(A_list.^2/2);            % jammer power / signal power per sample
B0 = nan(numel(M_list), numel(A_list)); B1 = B0; Th0 = B0; Th1 = B0;
for i = 1:numel(M_list)
    M = M_list(i);
    [Th0(i,:), Th1(i,:)] = ber_theory_jam(M, A_list, w_in, Pn);
    for k = 1:numel(A_list)
        o = struct('A', A_list(k), 'w0', w_in, 'min_err', 1000, ...
                   'max_bits', round(bits_for(M)/4));
        o.use_pn = false; B0(i,k) = binary_comm_sim(M, a, o);
        o.use_pn = true;  B1(i,k) = binary_comm_sim(M, a, o);
    end
    fprintf('Experiment 2: M = %d done\n', M);
end

f = figure('visible', 'off');
cols = {'b', 'r', [0 0.6 0], 'm'};
h = zeros(1, numel(M_list));
for i = 1:numel(M_list)
    h(i) = semilogy(jsr_db, max(B1(i,:), 1e-6), 'o', 'Color', cols{i}, 'MarkerSize', 5); hold on;
    semilogy(jsr_db, max(Th1(i,:), 1e-6), '-', 'Color', cols{i}, 'LineWidth', 1.1);
    semilogy(jsr_db, max(B0(i,:), 1e-6), 'x', 'Color', cols{i}, 'MarkerSize', 6);
    semilogy(jsr_db, max(Th0(i,:), 1e-6), '--', 'Color', cols{i}, 'LineWidth', 0.8);
    jA = 10*log10(M_list(i)^2/2);              % A = M for this M
    semilogy([jA jA], [1e-6 1], '-.', 'Color', cols{i}, 'LineWidth', 0.6);
end
grid on; ylim([1e-5 1]); xlim([min(jsr_db) max(jsr_db)]);
xlabel('Jammer-to-signal power ratio A^2/2 (dB)'); ylabel('Bit error rate');
title(sprintf('BER vs jammer power, \\omega_0 = %.3f rad/sample (o PN, x no PN)', w_in));
legend(h, arrayfun(@(m) sprintf('M = %d', m), M_list, 'UniformOutput', false), ...
       'Location', 'southeast');
save_fig(f, 'p128_ber_vs_jsr');

%% 3. BER against jammer frequency
M  = 100;
w_list = pi * linspace(0.002, 1, 120);     % simulated points
w_fine = pi * linspace(0.001, 1, 2000);    % theory curve
A_pair = [M 5];
W0 = zeros(2, numel(w_list)); W1 = W0;
Wt0 = zeros(2, numel(w_fine)); Wt1 = Wt0;
for q = 1:2
    for k = 1:numel(w_list)
        o = struct('A', A_pair(q), 'w0', w_list(k), 'min_err', 500, 'max_bits', 4e4);
        o.use_pn = false; W0(q,k) = binary_comm_sim(M, a, o);
        o.use_pn = true;  W1(q,k) = binary_comm_sim(M, a, o);
    end
    for k = 1:numel(w_fine)
        [Wt0(q,k), Wt1(q,k)] = ber_theory_jam(M, A_pair(q), w_fine(k), Pn);
    end
end
fprintf('Experiment 3 done\n');

f = figure('visible', 'off');
floor_ = 1e-5;
for q = 1:2
    subplot(2,1,q);
    h1 = semilogy(w_fine/pi, max(Wt0(q,:), floor_), 'r-', 'LineWidth', 0.8); hold on;
    h2 = semilogy(w_list/pi, max(W0(q,:), floor_), 'rx', 'MarkerSize', 5);
    h3 = semilogy(w_fine/pi, max(Wt1(q,:), floor_), 'b-', 'LineWidth', 1.4);
    h4 = semilogy(w_list/pi, max(W1(q,:), floor_), 'bo', 'MarkerSize', 4);
    grid on; ylim([floor_ 1]);
    title(sprintf('M = %d, A = %g, P_n = 1', M, A_pair(q)));
    ylabel('BER');
    if q == 1
        legend([h1 h2 h3 h4], 'No PN (theory)', 'No PN (sim)', 'PN (theory)', 'PN (sim)', ...
               'Location', 'southwest');
    end
end
xlabel('\omega_0 / \pi');
save_fig(f, 'p128_ber_vs_w0');

%% 4. Spectra before and after despreading
M = 100; A = 10; w0 = 0.2*pi; Nb = 400; L = Nb*M; Nfft = 1024;
b = bit_generator(Nb);
s = modulator(b, M);
p = pn_generator(L);
jam = sine_generator(L, A, w0);
w = noise_generator(L, a);
r = s.*p + jam + w;                         % received (spread) signal
z = r .* p;                                 % after despreading
[Psp, om] = psd_avg(s.*p, Nfft);            % transmitted spread signal
Pj  = psd_avg(jam, Nfft);                   % jammer in the channel
Pw  = psd_avg(w, Nfft);                     % noise in the channel
Ps  = psd_avg(s, Nfft);                     % desired signal after despreading
Pjd = psd_avg(jam.*p, Nfft);                % jammer after despreading
Pwd = psd_avg(w.*p, Nfft);                  % noise after despreading
Pz  = psd_avg(z, Nfft);                     % total after despreading
dB  = @(P) 10*log10(P + 1e-6);

f = figure('visible', 'off');
subplot(2,1,1);
h1 = plot(om/pi, dB(Psp), 'b'); hold on;
h2 = plot(om/pi, dB(Pj), 'r');
h3 = plot(om/pi, dB(Pw), 'Color', [0.5 0.5 0.5]);
grid on; ylim([-40 60]); ylabel('PSD (dB)');
title('In the channel (before despreading)');
legend([h1 h2 h3], 's(n)p(n)', 'i(n)', 'w(n)', 'Location', 'northeast');
subplot(2,1,2);
h1 = plot(om/pi, dB(Ps), 'b'); hold on;
h2 = plot(om/pi, dB(Pjd), 'r');
h3 = plot(om/pi, dB(Pwd), 'Color', [0.5 0.5 0.5]);
grid on; ylim([-40 60]); ylabel('PSD (dB)'); xlabel('\omega / \pi');
title('After despreading (input to the demodulator)');
legend([h1 h2 h3], 's(n)', 'i(n)p(n)', 'w(n)p(n)', 'Location', 'northeast');
save_fig(f, 'p128_spectra');

%% Save numeric results
res128 = struct('M_list', M_list, 'w_in', w_in, 'w_out', w_out, 'Pn', Pn, ...
    'T', T, 'Tth', Tth, 'A_list', A_list, 'jsr_db', jsr_db, 'B0', B0, 'B1', B1, ...
    'Th0', Th0, 'Th1', Th1, 'w_list', w_list, 'A_pair', A_pair, ...
    'W0', W0, 'W1', W1, 'Wt0', Wt0, 'w_fine', w_fine, 'Wt1', Wt1, ...
    'psd_om', om, 'Psp', Psp, 'Pj', Pj, 'Pw', Pw, 'Ps', Ps, 'Pjd', Pjd, 'Pwd', Pwd, ...
    'psd_M', M, 'psd_A', A, 'psd_w0', w0);
save(fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'res128.mat'), 'res128', '-v7');
write_table_128(res128);
