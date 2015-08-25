function opt_elems = opt_build(file)
% OPT_BUILD - Builds optical system specified in FILE
% 
% Calling:
% OPT_ELEMS = OPT_BUILD(FILE)
% 
% See also README OPT.EXMPL

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

fp = fopen(file,'r');

if fp == -1
  
  error(['Error, could not open file: ',file])
  
end

curr_line = fgetl(fp);
curr_type = strtok(curr_line(10:end));
curr_line = fgetl(fp);

while ~feof(fp) & ~strcmp(curr_type,'end')
  
  opt_args = '';
  % Each optical element starts with ``#type''
  % Read all parameters specifying that optical
  % element.
  while ~strcmp(curr_line(1:5),'#type') & ~strcmp(curr_type,'end')
    
    opt_args = str2mat(opt_args,curr_line);
    curr_line = fgetl(fp);
    
  end
  %opt_args
  % Build the current optical element.
  switch curr_type
   case 'file'
    opt_elem = opt_build(deblank(fliplr(deblank(fliplr(opt_args(end,10:end))))));
    %curr_type
   case 'aperture'
    opt_elem = opt_aperture(curr_type,opt_args);
   case 'grid'
    opt_elem = opt_grid(curr_type,opt_args);
   case 'lens'
    opt_elem = opt_lens(curr_type,opt_args);
   case 'prism'
    opt_elem = opt_prism(curr_type,opt_args);
   case 'screen'
    opt_elem = opt_screen(curr_type,opt_args);
   case 'slit'
    opt_elem = opt_slit(curr_type,opt_args);
   otherwise
    %curr_type
    %opt_args
    opt_elem = opt_fcn(curr_type,opt_args);
  end
  %keyboard
  if ~exist('opt_elems','var')
    
    opt_elems = [opt_elem];
    
  else
    
    opt_elems = [opt_elems opt_elem];
    
  end
  curr_type = strtok(curr_line(10:end));
  curr_line = fgetl(fp);
  %curr_type
  %curr_line
end

fclose(fp);
