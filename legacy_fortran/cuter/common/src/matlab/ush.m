function h=ush(x)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% USH Evaluate Hessian of objective function at x.
%
%   h=ush(x) returns the Hessian matrix in h, where h is a sparse matrix.
%   Since the Hessian is symmetric, h contains only the upper triangular
%   entries.
%
h=utools('ush',x);
