function [g,h]=ugrsh(x)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UGRSH Evaluate both gradient and Hessian of objective function at x.
%
%   [g,h]=ugrsh(x) returns gradient vector in g and Hessian matrix in h.
%   h is a sparse matrix containing only the upper triangular entries.
%
[g,h]=utools('ugrsh',x);
