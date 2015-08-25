%run on demo data

clc
clear
mesh_density=20;    % number of points generated 
n=100;               % length of the solution x
mesh_rho=20;         % NO points between 0 and 1 to consider for rho
mesh_delta=20;       % NO points between 0 and 1 to consider for delta 
tests_per_point=20;  % NUMBER OF TRIALS/ SIMULATIONS MADE

[h1]=Figure(mesh_density,n,mesh_rho, mesh_delta, tests_per_point)