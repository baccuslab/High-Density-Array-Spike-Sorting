function [jj] = opt_findstr(H,s)
% OPT_FINDSTR - find string S in character array H, 
%   without case sensitivity.
% 
% Calling:
%  [jj] = opt_findstr(H,s)

jj=[];%if not found return empty
ji = 1;
for ii = 1:size(H,1)
  if findstr(lower(H(ii,:)), lower(s)),
    jj(ji) = ii;
    ji = ji+1;
  end
end
