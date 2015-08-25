
function checkForForbiddenParameter(parameterCell, in)
    for i=1:2:length(in)
        for ii=1:length(parameterCell) 
            if strcmp(in{i}, parameterCell{ii})
                error('This Parameter is not allowed to be set by the user! (%s)', in{i});
            end
        end
    end