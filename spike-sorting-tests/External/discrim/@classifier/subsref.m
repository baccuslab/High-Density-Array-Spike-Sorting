function g = subsref(f, s)
%CLASSIFIER/SUBSREF Subscripted reference of CLASSIFIER object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

switch s(1).type
  case '.'
   switch s(1).subs
    case 'prior'
     if isempty(f.prior)
       n = f.counts;
       h = n/sum(n);
     elseif length(f.prior) == 1
       n = length(f.counts);
       h = repmat(1/n, 1, n);
     else
       h = f.prior;
     end
    case 'counts'
     h = f.counts;
    case 'nvar'
     h = size(f.range, 2);
    case 'nclass'
     h = length(f.counts);
    case 'range'
     h = f.range;
    otherwise
     error(sprintf('Reference to non-existent field ''%s''.', ...
		   s(1).subs))     
   end
   
   if length(s) > 1
     g = subsref(h, s(2:end));
   else
     g = h;
   end
 case '{}'
  error(['Cell array subscript reference not defined for CLASSIFIER' ...
	 ' objects.'])
 case '()'
  h = classify(f, s(1).subs{:});
  
  if length(s) > 1
    g = subsref(h, s(2:end));
  else
    g = h;
  end
end



