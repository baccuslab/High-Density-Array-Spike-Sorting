function out_v = opt_rot(in_v,rot_v,theta)
% OPT_ROT - rotate vector IN_V THEDA radians around ROT_V
%   
% Calling:
% out_v = opt_rot(in_v,rot_v,theta)

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

transpose_after = 1;

if size(in_v)==[3 1]
  in_v = in_v';
  rot_v = rot_v';
  transpose_after = 0;
end
rotv = rot_v/norm(rot_v);
e1 = rotv;
out_v = in_v;

if norm(in_v)
  
  e2 = cross(e1,in_v/norm(in_v));
  
  if norm(e2)
    
    e2 = e2/norm(e2);
    e3 = cross(e1,e2);
    
    trmtr = [e1;e2;e3];
    
    rmatr = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
    
    out_v = trmtr'*rmatr*trmtr*in_v';
    
    if transpose_after
      
      out_v = out_v';
      
    end
    
  end
  
end
