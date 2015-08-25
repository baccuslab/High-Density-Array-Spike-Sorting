% Calls all unit tests

addpath(fullfile(pwd, 'matlab_xunit_3.1','matlab_xunit', 'xunit'));
runtests(fullfile('tests', 'epoch'))
runtests(fullfile('tests', 'util'))
