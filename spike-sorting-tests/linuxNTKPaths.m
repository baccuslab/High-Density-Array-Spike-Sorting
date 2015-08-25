p = {'/usr/local/bsse/etc/matlab' 
     '/usr/local/hima/hidens/current/matlab/mex_ntkparser'
     '/usr/local/bsse/el6/etc/matlab'};
 
for i=1:length(p)
    addpath(p{i});
end