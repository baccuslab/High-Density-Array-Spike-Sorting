% Startup script, put this into your matlab startup dir or your matlabrc
% and set this path to where the matlab mysort directory is:
global MYSORT_PATH;
MYSORT_PATH = 'PUT PATH HERE';
% addpath('/net/bs-sw/sw-repo/hierlemann/hidens/14154/matlab/mex_ntkparser')
addpath('/home/frankef/bel.svn/cmosmea_external/ntk_tools/trunk/mex_ntkparser')

this_path = pwd;
disp(['Startup in ' this_path])

cd(MYSORT_PATH)
disp(['Changed path to ' pwd])
setup();