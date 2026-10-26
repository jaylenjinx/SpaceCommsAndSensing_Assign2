function write_table_127(r)
% WRITE_TABLE_127  Write the sweep (a) results as a LaTeX table body,
%   including the error count and a 95% Clopper-Pearson interval.
    outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
    fid = fopen(fullfile(outdir, 'table_127.tex'), 'w');
    [lo, hi] = ber_ci(r.nerr_a, r.nbits_a);
    for k = 1:numel(r.M_list)
        fprintf(fid, '%d & %.2f & %d & %s & %s & [%s, %s] & %s & %s \\\\\n', ...
            r.M_list(k), r.snr_a_db(k), r.nerr_a(k), bits(r.nbits_a(k)), ...
            sci(r.ber_a(k)), sci(lo(k)), sci(hi(k)), sci(r.exact_a(k)), ...
            sci(r.gauss_a(k)));
    end
    fclose(fid);
end

function s = sci(x)
    if x == 0, s = '$0$'; return; end
    e = floor(log10(x)); m = x / 10^e;
    if e >= -2, s = sprintf('$%.4f$', x);
    else, s = sprintf('$%.2f\\times10^{%d}$', m, e); end
end

function s = bits(n)
    e = floor(log10(n)); m = n / 10^e;
    s = sprintf('$%.1f\\times10^{%d}$', m, e);
end
