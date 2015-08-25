
function savefig(fig_handle, fname, varargin)
    % Saves the figure in the provided figure handle or the current figure to
    % a) a png file and b) a .fig file. Since the default dpi of 100 might
    % be changed during this process to increase the quality of the .png file
    % it will set all "unit" properties of axes/annotations etc in the
    % figure to "normalized" before calling print.
    P.dpi = 300;
    P.fig = 1;
    P.png = 1;
    P.eps = 0;
    P.emf = 0;
    P.ai = 0;
    P.special = 0;
%     if ~ishandle(fig_handle)
%         varargin = [fname varargin];
%         fname = fig_handle;
%         fig_handle = gcf;
%     end
dev = {'win'
    'winc'
    'meta'
    'bitmap'
    'setup'
    'ps'
    'psc'
    'ps2'
    'ps2c'
    'psc2'
    'eps'
    'epsc'
    'eps2'
    'eps2c'
    'epsc2'
    'hpgl'
    'ill'
    'mfile'
    'tiff'
    'tiffnocompression'
    'bmp'
    'hdf'
    'png'
    'svg'
    'jpeg'
    'laserjet'
    'ljetplus'
    'ljet2p'
    'ljet3'
    'ljet4'
    'pxlmono'
    'cdjcolor'
    'cdjmono'
    'deskjet'
    'cdj550'
    'djet500'
    'cdj500'
    'paintjet'
    'pjetxl'
    'pjxl'
    'pjxl300'
    'dnj650c'
    'bj10e'
    'bj200'
    'bjc600'
    'bjc800'
    'epson'
    'epsonc'
    'eps9high'
    'ibmpro'
    'pcxmono'
    'pcxgray'
    'pcx16'
    'pcx256'
    'pcx24b'
    'bmpmono'
    'bmp16m'
    'bmp256'
    'pngmono'
    'pnggray'
    'png16m'
    'png256'
    'pbm'
    'pbmraw'
    'pgm'
    'pgmraw'
    'ppm'
    'ppmraw'
    'pkm'
    'pkmraw'
    'tifflzw'
    'tiffpack'
    'tiff24nc'
    'pdfwrite'};
    P = mysort.util.parseInputs(P, 'savefig', varargin);    
    global savefig__;
    if ~isempty(savefig__) && savefig__ == 0
        return
    end
    if nargin == 1
        fname = fig_handle;
        fig_handle = gcf();
    end
    set(fig_handle, 'PaperPositionMode', 'auto');   % Use screen size
    visibility = get(fig_handle, 'visible');
    set(fig_handle, 'visible', 'on');
    mysort.plot.figureChildrenSet(fig_handle, 'Units', 'normalized');
    if P.png
        print(fig_handle, ['-r' num2str(P.dpi)], '-dpng', [fname '.png']);
    end
    if P.emf
        saveas(fig_handle, [fname '.emf'], 'emf');
    end
    if P.fig
        saveas(fig_handle, [fname '.fig'], 'fig');
    end
    if P.eps
        saveas(fig_handle, [fname '.eps'], 'eps');
    end
    if P.ai
        saveas(fig_handle, [fname '.ai'], 'ai');
    end    
    if P.special
        print(fig_handle, [fname '.pdf'], ['-d' 'pdfwrite'] )
    end
    set(fig_handle, 'visible', visibility);