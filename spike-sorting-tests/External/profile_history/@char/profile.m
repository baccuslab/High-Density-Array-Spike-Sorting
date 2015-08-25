function varargout = profile(varargin)
    % wrapper function from the builtin profile() function, with extra -timing input arg

    % Unfortunately, the following fails since profile() is an m-file, not built-in:
    %varargout = builtin('profile','-timing',varargin{:});

    % Also unfortunately, @char/profile.m precedes the Matlab path (http://mathworks.com/help/matlab/matlab_prog/function-precedence-order.html)
    % So copy the original profile.m function here, rename it and then use it:
    curFolder = fileparts(mfilename('fullpath'));
    copyfile([matlabroot '/toolbox/matlab/codetools/profile.m'], [curFolder '/profile_orig.m']);
    if nargout
        if any(strcmpi(varargin,'on'))
            varargout{:} = profile_orig('-timestamp',varargin{:});
        else
            varargout{:} = profile_orig(varargin{:});
        end
    else  % no output args requested
        if any(strcmpi(varargin,'on'))
            profile_orig('-timestamp',varargin{:});
        else
            profile_orig(varargin{:});
        end
    end
end
