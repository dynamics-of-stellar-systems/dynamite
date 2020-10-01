function [f,g]=cofg(x)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% COFG Evaluate objective function and possibly its gradient at x.
%
%   f=cofg(x) returns the value of the objective function in f.
%
%   [f,g]=cofg(x) also returns the gradient vector in g.
%
if nargout == 1
   f=ctools('cofg',x);
else
   [f,g]=ctools('cofg',x);
end
