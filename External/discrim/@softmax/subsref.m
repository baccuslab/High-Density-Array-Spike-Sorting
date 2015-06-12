function g = subsref(f, s)
%SOFTMAX/SUBSREF Access fields of SOFTMAX object or classify new data.

%   $Id: subsref.m,v 1.1 1999/06/04 18:50:50 michael Exp $
%   $Log: subsref.m,v $
%   Revision 1.1  1999/06/04 18:50:50  michael
%   Initial revision
%

switch s(1).type
  case '.'
   switch s(1).subs
    case 'weights'
     h = f.weights;
    case 'classifier'
     h = f.classifier;
    otherwise
     h = subsref(f.classifier, s(1));
   end
   
   if length(s) > 1
     g = subsref(h, s(2:end));
   else
     g = h;
   end
   
 case '{}'
  subsref(f.classifier, s)
 case '()'
  h = classify(f, s(1).subs{:});
  
  if length(s) > 1
    g = subsref(h, s(2:end));
  else
    g = h;
  end
end
