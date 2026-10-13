function save_fig(f, name)
% SAVE_FIG  Write figure f to ../figures/<name>.pdf and .png.
%   Uses exportgraphics in MATLAB R2020a+ and print otherwise (Octave).
    outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    base = fullfile(outdir, name);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
        exportgraphics(f, [base '.png'], 'Resolution', 200);
    else
        set(f, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.4 4.2], ...
               'PaperSize', [6.4 4.2]);
        print(f, '-dpdf', [base '.pdf']);
        print(f, '-dpng', '-r200', [base '.png']);
    end
end
