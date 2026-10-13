function write_table_128(r)
% WRITE_TABLE_128  Write the A = M results as a LaTeX table body.
    outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
    fid = fopen(fullfile(outdir, 'table_128.tex'), 'w');
    for i = 1:numel(r.M_list)
        fprintf(fid, '%d & %d', r.M_list(i), r.M_list(i));
        for c = 1:4
            fprintf(fid, ' & %.4f\\,(%.4f)', r.T(i,c), r.Tth(i,c));
        end
        fprintf(fid, ' \\\\\n');
    end
    fclose(fid);
end
