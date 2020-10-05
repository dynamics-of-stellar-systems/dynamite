function [g,h]=ugrdh(x)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UGRDH Evaluate both gradient and Hessian of objective function at x.
%
%   [g,h]=ugrdh(x) returns the gradient vector in g and the Hessian matrix
%   in h.  The Hessian is stored as a full matrix.
%
[g,h]=utools('ugrdh',x);
