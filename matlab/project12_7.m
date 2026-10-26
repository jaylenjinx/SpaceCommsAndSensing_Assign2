%% Project 12.7: Binary data communications system (Ingle & Proakis 4e)
% Simulates the system of Fig. 12.16:
%   bit generator -> modulator (M samples/bit) -> + uniform noise ->
%   integrate-and-dump demodulator -> detector & error counter
% and measures the bit error rate (BER) against SNR = M / Pn, Pn = a^2/3.
%
% Two sweeps are run, as suggested in the brief:
%   (a) Pn fixed (= 1) while M is changed,
%   (b) M fixed (= 10) while Pn (i.e. a) is changed.
% Results are compared with the Gaussian (CLT) approximation
% Pe ~ Q(sqrt(M/Pn)) and with the exact result from exact_ber_uniform.
clear; close all; clc;
rng(12);                                   % reproducible results
if exist('OCTAVE_VERSION', 'builtin'), graphics_toolkit('gnuplot'); end

min_err  = 200;                            % stop after this many errors
max_bits = 2e7;                            % ... or this many bits

%% Illustration of the signals through each block
M = 10; a = sqrt(3);                       % Pn = 1, SNR = 10 dB
b = bit_generator(8);
s = modulator(b, M);
w = noise_generator(numel(s), a);
r = s + w;
y = demodulator(r, M);
b_hat = detector(y, b);
n = 0:numel(s)-1;

f = figure('visible', 'off');
subplot(3,1,1); stairs(n, s, 'LineWidth', 1.2); ylim([-1.5 1.5]); grid on;
title(sprintf('Modulator output s(n), bits = [%s]', num2str(b)));
ylabel('s(n)');
subplot(3,1,2); plot(n, r, '.-'); grid on;
title(sprintf('Received r(n) = s(n) + w(n), a = %.2f (P_n = %.2f)', a, a^2/3));
ylabel('r(n)');
subplot(3,1,3); stem(0:numel(y)-1, y, 'filled'); hold on;
plot([-0.5 numel(y)-0.5], [0 0], 'k--'); grid on;
title(sprintf('Demodulator output y(k), detected bits = [%s]', num2str(b_hat)));
xlabel('bit index k'); ylabel('y(k)');
save_fig(f, 'p127_waveforms');

%% (a) Pn fixed, M varied
a_fix = sqrt(3);  Pn_fix = a_fix^2/3;      % Pn = 1
M_list = [1 2 3 4 6 8 10 12 14 16 18 20];
ber_a = zeros(size(M_list)); nb_a = ber_a; ne_a = ber_a;
for k = 1:numel(M_list)
    [ber_a(k), ne_a(k), nb_a(k)] = binary_comm_sim(M_list(k), a_fix, ...
        struct('min_err', min_err, 'max_bits', max_bits));
    fprintf('(a) M = %2d  SNR = %5.2f dB  BER = %.3e  (%d bits)\n', ...
        M_list(k), 10*log10(M_list(k)/Pn_fix), ber_a(k), nb_a(k));
end
snr_a   = M_list / Pn_fix;
gauss_a = qfunc_local(sqrt(snr_a));
exact_a = zeros(size(M_list));
for k = 1:numel(M_list), exact_a(k) = exact_ber_uniform(M_list(k), a_fix); end

%% (b) M fixed, Pn varied
M_fix  = 10;
snr_db = 0:1:13;
Pn_b   = M_fix ./ 10.^(snr_db/10);
a_b    = sqrt(3*Pn_b);
ber_b  = zeros(size(snr_db)); nb_b = ber_b; ne_b = ber_b;
for k = 1:numel(snr_db)
    [ber_b(k), ne_b(k), nb_b(k)] = binary_comm_sim(M_fix, a_b(k), ...
        struct('min_err', min_err, 'max_bits', max_bits));
    fprintf('(b) a = %5.2f  SNR = %5.2f dB  BER = %.3e  (%d bits)\n', ...
        a_b(k), snr_db(k), ber_b(k), nb_b(k));
end
gauss_b = qfunc_local(sqrt(10.^(snr_db/10)));
exact_b = exact_ber_uniform(M_fix, a_b);

%% Plot BER vs SNR for both sweeps
snr_fine = linspace(-1, 14, 200);
f = figure('visible', 'off');
semilogy(10*log10(snr_a), ber_a, 'o', 'MarkerSize', 7, 'LineWidth', 1.2); hold on;
semilogy(snr_db, ber_b, 's', 'MarkerSize', 7, 'LineWidth', 1.2);
semilogy(snr_fine, qfunc_local(sqrt(10.^(snr_fine/10))), 'k-', 'LineWidth', 1.2);
semilogy(10*log10(snr_a), exact_a, 'b--', 'LineWidth', 1);
semilogy(snr_db, exact_b, 'r:', 'LineWidth', 1.4);
grid on; xlim([-1 14]); ylim([1e-8 1]);
xlabel('SNR = M / P_n (dB)'); ylabel('Bit error rate');
title('Project 12.7: BER of the binary system with uniform noise');
legend('Sim. (a): P_n = 1, M varied', 'Sim. (b): M = 10, P_n varied', ...
       'Gaussian approx. Q(\surd SNR)', 'Exact (a)', 'Exact (b)', ...
       'Location', 'southwest');
save_fig(f, 'p127_ber');

%% Save numeric results
res127 = struct('M_list', M_list, 'snr_a_db', 10*log10(snr_a), 'ber_a', ber_a, ...
    'nbits_a', nb_a, 'nerr_a', ne_a, 'exact_a', exact_a, 'gauss_a', gauss_a, ...
    'M_fix', M_fix, 'snr_b_db', snr_db, 'a_b', a_b, 'ber_b', ber_b, ...
    'nbits_b', nb_b, 'nerr_b', ne_b, 'exact_b', exact_b, 'gauss_b', gauss_b, ...
    'wave_b', b, 'wave_s', s, 'wave_r', r, 'wave_y', y, 'wave_bhat', b_hat, 'wave_M', 10);
save(fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'res127.mat'), 'res127', '-v7');
write_table_127(res127);
