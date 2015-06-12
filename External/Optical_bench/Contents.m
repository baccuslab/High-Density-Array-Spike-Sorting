%  OPT_BENCH - ray-tracer through optical systems
% Version 1.0 20100429
%   
% Ray tracer functions
%   opt_build        - Builds optical system specified in FILE
%   opt_trace        - ray tracing through optical systems.
%   opt_intersection - Determine the impact point of an optical
%   opt_refraction   - calculation of optical refraction - Snells law.
%
% Glass characterisation
%   opt_refrindx     - Calculate refractive index 
%   opt_absorption   - Calculate the absorption coefficient
%   
% Ploting
%   opt_plotoptics   - Plot the optical system.
%   
% Default structures
%   opt_ray          - RAY creator for OPT_TOOLS
%   opt_elem         - Default opt_elem structure
%   
% Optical elements
%   opt_aperture     - Circular aperture, iris
%   opt_flens        - Thick aproximation of thin lens f focal width
%   opt_grid         - Simple difractive grating, 
%   opt_lens         - Spherical thick lens, single glass.
%   opt_prism        - Wedge prism
%   opt_screen       - Screen - imaging detector.
%   opt_slit         - optical slit.
%   
% Example scripts
%   opt_exmpl_astigmatic   - astigmatic aberation example
%   opt_exmpl_chromatic_ab - chromatic aberation example
%   opt_exmpl_coma         - example showing coma
%   opt_exmpl_field_curv   - example showing field curvature
%   opt_exmpl_sphere_coma  - example showing spherical aberation and coma
