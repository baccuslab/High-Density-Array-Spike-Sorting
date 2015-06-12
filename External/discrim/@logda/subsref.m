function g = subsref(f, s)
%LOGDA/SUBSREF  Subscripted reference of LOGDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

switch s(1).type
  case '.'
   switch s(1).subs
    case 'coefs'
     h = f.coefs;
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
