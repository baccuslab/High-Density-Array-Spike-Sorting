
function a = figureDate(str)
    if nargin == 0
        str = '';
    else
        str = [str ' - '];
    end
    cn = mysort.util.getCallerName();
%     if ~isempty(strcmp(cn, '(shell?)'))
% Try to get the name of the last run m file or the one that is open in
% matlab editor... but how?
% This does not work as it will always be figureDate
%         cn = mfilename;
%     end
    
    str = [str date() ' - ' cn];
    a = annotation('textbox',[.1 .1 .1 .1],...
        'Units','pixel',...
        'Position', [1 10 400 12],...
        'FitBoxToText','off',...
        'String',str,...
        'Interpreter', 'none',...
        'LineStyle','none');
    set(a, 'Units', 'normalized')
    p = get(a, 'Position');
    set(a, 'Position', [0 p(2) 1 p(4)]);
end