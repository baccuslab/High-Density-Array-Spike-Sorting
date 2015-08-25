function outpath = figureOutPath()
    if ~isempty(strfind(computer, 'WIN'))
        outpath   = 'C:\LocalData\HarrisBuzsaki\BotmPaperFigures\';
    else
        outpath = '/net/bs-filesvr01/export/group/hierlemann/recordings/collaborations/Harris/BotmPaperFigures/';
    end

end