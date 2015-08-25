function [opt1] = opt_fcn(opt_type,opt_args)
% OPT_FCN optical surface with general shape
% OPT_TYPE should be the function name, opt_args
% should be a string matrix, see README_OPT for specification.
% The specification and handling of OPT_ARGS are left to mercy of
% the person writing the OPT_FCN. See the function to call for
% details. A requirement on the optical function is that it should
% return 0 when called with OPT_FCN(R,'s',optargs) when R is a
% point on the surface and growing positive and negative values
% when R is not on the surface. When called with
% OPT_FCN(R,'n',opt_args) the function should return the surface
% normal.
% 
% Calling:
% [opt1] = opt_fcn(opt_type,opt_args)
% 
% See also OPT_LENS, OPT_SCREEN, OPT_GRID, OPT_PRISM, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

if nargin ~= 2,
  help opt_fcn;
  ok = 0;
  return;
end;

%opt_type
opt1 = opt_elem(opt_type);

opt1.fcn1 = inline([opt1.type,'(point_on_line(r_0,e_l,l),s_or_n,arglist)^2'],'l','r_0','e_l','s_or_n','arglist');
opt1.fcn2 = inline([opt1.type,'(r_int,s_or_n,arglist)'],'r_int','s_or_n','arglist');

ii = opt_findstr(opt_args,'glass');
opt1.glass = strtok(opt_args(ii,12:end));

opt_args = opt_args([1:ii-1 ii+1:end],:);

opt1.arglist = struct;
for ii = 1:size(opt_args,1)
  
  if ~isempty(deblank(opt_args(ii,:)))
    val = str2num(opt_args(ii,13:end));
    if isempty(val)
      val = fliplr(deblank(fliplr(deblank(opt_args(ii,13:end)))));
    end
    opt1.arglist = setfield(opt1.arglist,deblank(opt_args(ii,1:11)),val);
  end
  
end
