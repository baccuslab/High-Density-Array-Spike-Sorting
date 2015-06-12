function outcell=parseSize(sz,numargsout,dim)
%
%by Matt Jacobson, 
%Copyright, The University of Michigan, 2005
%
%parseSize - this function is an OOP aid for overloading the size() method in 
%user-defined classes. If the user provides the vector of dimensions of a hypothetical
%array object, parseSize will take care of the main input/output argument processing 
%in a way that imitates the conventions of MATLAB's usual SIZE function.
%
%EXAMPLE: consider the following class and notice how parseSize is
%used to implement the SIZE method in a single line,
%
%     classdef myClass
%         properties 
%             dimension
%         end % properties
%
%         methods
%             function varargout=size(obj,varargin)
%                 varargout=parseSize(obj.dimension,nargout,varargin{:});
%             end
%         end
%     end 
%
%
%In the next several examples, we see that all of the same calling syntaxes
%as for the usual MATLAB size() function are now enabled for this class:
%
%     >> obj=myClass; obj.dimension=[3,5,4,1];
% 
%     >> [m,n,p,q,r,s]=size(obj); [m,n,p,q,r,s]
% 
%      ans =
% 
%          3     5     4     1     1     1
% 
%     >> [m,n]=size(obj); [m,n]
% 
%      ans =
% 
%          3    20
% 
%     >> s=size(obj)
% 
%     s =
% 
%          3     5     4
% 
% 
%     >> s2=size(obj,2)
% 
%     s2 =
% 
%          5
% 
%
%SYNTAX:
%
%  argsout=parseSize(dimensions,numargsout,dim)
%
%in:
%
% dimensions:  A vector specifying the dimensions of a hypothetical array-like
%              object. Trailing ones are permitted.
%
% numargsout: The number of requested outputs in a call to the class'
%             size method, e.g., for [m,n,p]=size(obj) a value of 
%             numargsout=3 should be passed to parseSize.
%
% dim: if the calling syntax to the class's size method is size(obj,i), 
%      one would pass dim=i to parseSize.
%
%out:
%
% argsout: a cell array to be passed as the varargout of the size method. 




 idx=find(sz~=1,1,'last');
 if isempty(idx)
     sz=[1,1];
 elseif idx==1
     sz=[sz(1), 1];
 else
   sz=sz(1:idx);
 end
         
         

%dimSpecified=isvar('dim');
dimSpecified=logical(exist('dim','var'));


if dimSpecified,

  if ~isscalar(dim), error('DIM argument must be scalar.'); end
  
  if dim>length(sz), sz(end+1:dim)=1; end
  
  sz=sz(dim);

end


if numargsout<=1,

    outcell{1}=sz;
    return;

elseif dimSpecified,

  error('Too many output arguments.')

else

 
  if length(sz)>=numargsout

    sz(numargsout)=prod(sz(numargsout:end));
    sz(numargsout+1:end)=[];
    
  else
    
    sz(end+1:numargsout)=1;  
      
  end

  outcell=num2cell(sz);

end

 end
 