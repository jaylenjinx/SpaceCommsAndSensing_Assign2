function save_fig(f, name, h)
% SAVE_FIG  Write figure f to ../figures/<name>.pdf and .png.
%   h is the page height in inches for print (default 4.2).
%   Uses exportgraphics in MATLAB R2020a+ and print otherwise (Octave).
    if nargin < 3, h = 4.2; end
    outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    base = fullfile(outdir, name);
    if exist('exportgraphics', 'file') > 0
        exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
        exportgraphics(f, [base '.png'], 'Resolution', 200);
    else
        set(f, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.4 h], ...
               'PaperSize', [6.4 h]);
        print(f, '-dpdf', [base '.pdf']);
        print(f, '-dpng', '-r200', [base '.png']);
    end
end
