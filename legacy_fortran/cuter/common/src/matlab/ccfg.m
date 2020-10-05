function [c,cjac]=ccfg(x,jtrans)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CCFG Evaluate constraint functions and possibly their gradients at x.
%
%   c=ccfg(x) returns the vector of constraint values in C.
%
%   [c,cjac]=ccfg(x) also returns the constraint Jacobian matrix in cjac.
%   cjac(i,j) contains the partial derivative of the i-th constraint
%   with respect to the j-th variable.
%
%   [c,cjac]=ccfg(x,jtrans), with jtrans set to 1, gives the transpose
%   of the Jacobian, where the i,j-th component is the partial derivative
%   of the j-th constraint with respect to the i-th variable.
%   If jtrans is not given, it is assumed to be 0.
%
if nargout == 1 & ( nargin == 1 | nargin == 2 )
   c=ctools('ccfg',x);
elseif nargout == 2 & nargin == 1
   [c,cjac]=ctools('ccfg',x);
else
   [c,cjac]=ctools('ccfg',x,jtrans);
end
