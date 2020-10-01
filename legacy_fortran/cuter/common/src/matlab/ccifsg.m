function [ci,sgci]=ccifsg(x,i)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CCIFSG Evaluate an individual constraint function and possibly its
% gradient at x, in a sparse format.
%
%   ci=ccifsg(x,i) returns the vector of i-th constraint value at x in ci.
%
%   [ci,sgci]=ccfsg(x,i) returns the sparse constraint gradient
%   in sgci.
%
if nargout == 1 & nargin == 2
   ci=ctools('ccifsg',x,i);
elseif nargout == 2 & nargin == 2
   [ci,sgci]=ctools('ccifsg',x,i);
end
