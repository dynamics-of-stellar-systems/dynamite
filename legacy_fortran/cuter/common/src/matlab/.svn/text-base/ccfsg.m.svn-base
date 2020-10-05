function [c,cjac]=ccfsg(x,jtrans)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CCFSG Evaluate constraint functions and possibly their gradients at x.
%
%   c=ccfsg(x) returns the vector of constraint values in c.
%
%   [c,cjac]=ccfsg(x) returns the sparse constraint Jacobian matrix in cjac.
%   The i,j-th nonzero entry in cjac corresponds to the partial derivative
%   of the i-th constraint with respect to the j-th variable.
%
%   [c,cjac]=ccfsg(x,jtrans), with jtrans set to 1, gives the transpose
%   of the Jacobian, where the i,j-th component is the partial derivative
%   of the j-th constraint with respect to the i-th variable.
%   If jtrans is not given, it is assumed to be 0.
%
if nargout == 1 & ( nargin == 1 | nargin == 2 )
   c=ctools('ccfsg',x);
elseif nargout == 2 & nargin == 1
   [c,cjac]=ctools('ccfsg',x);
else
   [c,cjac]=ctools('ccfsg',x,jtrans);
end
