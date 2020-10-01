function [f,g]=uofg(x)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UOFG Evaluate objective function and possibly its gradient at x.
%
%   f=uofg(x) returns the value of the objective function in f.
%
%   [f,g]=uofg(x) also returns the gradient vector in g.
%
if nargout == 1
   f=utools('uofg',x);
else
   [f,g]=utools('uofg',x);
end
