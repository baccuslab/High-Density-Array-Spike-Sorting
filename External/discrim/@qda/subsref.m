function g = subsref(f, s)
%QDA/SUBSREF Subscripted reference of LDA object.

%   Copyright (c) 1999 Michael Kiefte. 

%   $Log$

switch s(1).type
 case '.'
  switch s(1).subs
   case 'means'
    h = f.means;
   case 'scale'
    h = f.scale;
   case 'ldet'
    h = f.ldet;
   case 'est'
    h = f.est;
   case 'nu'
    if f.est ~= 't'
      error(['Field ''nu'' only defined for t-estimator QDA' ...
	     ' objects.'])
    end
    h = f.nu;
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

