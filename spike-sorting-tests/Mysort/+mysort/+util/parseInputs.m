
function [P unexplained_vars]= parseInputs(P, in, treat_unexp_vars, dummy)
    % Parses the inputs to a function given in varargin. E.g. if a the 
    % function foo is defined as
    %
    % function f = foo(p1,p2,varargin)
    %   P.bar1 = -1;
    %   P.bar2 =  1;
    %   P = parseInputs(P, varargin);
    %   f = P.bar1 * P.bar2;
    % end
    %
    % then the field "bar1" and "bar2" of the struct P can be set by the 
    % caller of foo by providing the following two parameters to varargin:
    %
    % f = foo(p1,p2,'bar2',1, 'bar1', 5);
    %
    % This allows to define functions with default values to parameters in P
    % which can be overwritten by the caller in a human readable and
    % convinient way: The name of the parameter as string followed by the
    % actual value.
    %
    % Remark: 
    % If "otherP" is defined as input the next parameter will be handled
    % as a struct/class containing the parameter as fields/properties. This
    % allows to cascade a struct down through a caller hierarchy or to set
    % the properties of a class in a nice manner. Also, if the first
    % argument in varargin is not a string but a struct it is treated as a
    % "otherP".
    %
    % Inputs:
    %   P    - struct with default values as fields
    %   in   - varargin in of the calling function
    %   unexp_vars - can be either nothing (defaults to "split") or one of
    %                the following strings or corresponding numbers:
    %                "split" or 1 => not defined parameters are returned in
    %                a separate struct similar to P
    %                "merge" or 2 => undefined parameters are just attached
    %                to P as if they were defined
    %                "error" or -1 => an error is raised when an
    %                undefined parameter is given. good for debugging
    %
    % Output:
    %   P    - struct with the parameters. In "in" unspecified parameters
    %          remain on the default value.
    %   unexplained_vars - the struct containing the unexplained i.e.
    %                      by P not defined parameters
    %  Modified 16-Feb-12
    %  Author: Felix Franke
    %          felfranke@gmail.com
    unexplained_vars = struct();
    ERROR = -1;
    SPLIT = 1;
    MERGE = 2;
    
    % THIS IS FOR DOWNWARD COMPATIBILITY
    if ischar(in) && iscell(treat_unexp_vars)
        caller_name = in;
        in = treat_unexp_vars;
        if nargin < 4
            treat_unexp_vars = SPLIT;
        else
            treat_unexp_vars = dummy;
        end
    else
        caller_name = mysort.util.getCallerName();
    end
    % END: THIS IS FOR DOWNWARD COMPATIBILITY
    
    
    if nargin < 3; treat_unexp_vars = SPLIT; end
    switch lower(treat_unexp_vars)
        case 'error'
            treat_unexp_vars = ERROR;
        case 'split'
            treat_unexp_vars = SPLIT;
        case 'merge'
            treat_unexp_vars = MERGE;
        otherwise
            if ~any([ERROR SPLIT MERGE] == treat_unexp_vars)
                error('unknown code for treat_unexp_vars!')
            end
    end
    TUV = treat_unexp_vars;

    if isempty(in)
        return
    end
    
    if isstruct(in)
        % first and only argument is struct!
        [P unexplained_vars] = copyP(P,in);
        return
    elseif isstruct(in{1})
        % first argument is struct!
        [P unexplained_vars] = copyP(P,in{1});
        % remove from list
        in(1) = [];
    end
    fnames = fieldnames(P);
%     % Make all fieldnames with more than 2 latters lower case
%     for i=1:length(fnames)
%         if length(fnames{i}) > 2
%             fnames{i} = lower(fnames{i});
%         end
%     end
    for ii=1:2:length(in)
        if strcmp(in{ii}, 'otherP')
            [P unexplained_vars] = copyP(P,in{ii+1}, unexplained_vars);
        else
            str = in{ii};
%             if length(str) > 2
%                 str = lower(str);
%             end
            if ~isempty(find(strcmp(fnames, str),1)) || TUV==MERGE
                P.(str) = in{ii+1};
            elseif TUV == SPLIT
                unexplained_vars.(str) = in{ii+1};
            else
                error('Error! Parameter in %s unknown: %s', caller_name, str);
            end
        end
    end

    %------------------------------------------------------------------
    function [Pdest Punexp] = copyP(Pdest, Psource, Punexp)
        sourceNames = fieldnames(Psource);
        if nargin < 3; Punexp = []; end
        if isstruct(Pdest)
            for f=1:length(sourceNames)
                if isfield(Pdest, sourceNames{f}) || TUV==MERGE
                    if isstruct(Pdest.(sourceNames{f})) && isstruct(Psource.(sourceNames{f}))
                        [Pdest.(sourceNames{f}) un] = copyP(Pdest.(sourceNames{f}), Psource.(sourceNames{f}));  
                    else
                        Pdest.(sourceNames{f}) = Psource.(sourceNames{f});
                    end
%                 elseif TUV == MERGE
%                     Pdest.(sourceNames{f}) = Psource.(sourceNames{f});
                elseif TUV == SPLIT
                    Punexp.(sourceNames{f}) = Psource.(sourceNames{f});
                else
                    error('FieldName in %s unknown: %s', caller_name, sourceNames{f});
                end
            end
        else % is a class
            destNames = properties(Pdest);
            for f=1:length(sourceNames)
                if ~isempty(find(strcmp(destNames, sourceNames{f}),1)) || TUV==MERGE
                    try
                        Pdest.(sourceNames{f}) = Psource.(sourceNames{f});
                    catch ME
                        % Porperties can be private
                    end
                elseif TUV == SPLIT
                    Punexp.(sourceNames{f}) = Psource.(sourceNames{f});
                else
                    error('FieldName in %s unknown: %s', caller_name, sourceNames{f});                
                end
            end              
        end
    end
    %------------------------------------------------------------------
    function s = recursiveStructMerge(s, s2)
        sourceNames = fieldnames(s2);
       
    end
end