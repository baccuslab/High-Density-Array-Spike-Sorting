function varargout = parseopt(opt, varargin)
%PARSEOPT Get values from fields of option struct.
%   [V1, V2, ...] = parseopt(opt, F1, F2, ...) where opt is the
%   structure to be processed and F1, F2, ... are the fields to be
%   checked. The values of the fields F1, F2, ... are returned in
%   V1, V2, etc. If a field is not present in opt, the empty matrix
%   is returned in the corresponding output argument.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

narg = nargout;
varargout = cell(size(varargin));
for i = 1:narg
  if ~ischar(varargin{i})
    error('Field names must be strings.')
  elseif isfield(opt, varargin{i})
    varargout{i} = getfield(opt, varargin{i});
  else
    varargout{i} = [];
  end
end

  
