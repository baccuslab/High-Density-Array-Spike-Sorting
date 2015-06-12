function  [r_out] = point_on_line(r_0,e_l,l);
% POINT_ON_LINE  - calculates the vector to a point
% R_OUT that are L away from R_0 in the direction
% E_L.
% 
% Calling:
% [r_out] = poin_on_line(r_0,e_l,l);
% 

%       Bjorn Gustavsson
%	Copyright (c) 1997 by Bjorn Gustavsson
%	$Revision: 1.1 $  $Date: 1997/12/17 15:08:28 $

%
%       $Log: point_on_line.m,v $
%       Revision 1.1  1997/12/17 15:08:28  bjorn
%       Initial revision
%

r_out = r_0+l*e_l;
