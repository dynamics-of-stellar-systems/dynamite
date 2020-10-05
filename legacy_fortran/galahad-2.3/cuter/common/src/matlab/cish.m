function h=cish(x,i)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CISH Evaluate Hessian of individual problem function at x.
%
%   h=cish(x,i) returns the (sparse) Hessian of the i-th problem
%   function in h. The function is selected using the index i.
%   If i=0 -> objective function,
%      i>0 -> i-th constraint function.
%
%   h is a sparse matrix containing only the upper triangular entries.
%
h=ctools('cish',x,i);
