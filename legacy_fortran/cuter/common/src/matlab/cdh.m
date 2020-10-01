function h=cdh(x,v)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CDH Evaluate Hessian of Lagrangian function at x and v.
%
%   cdh(x,v) returns the Hessian of the Lagrangian function,
%   stored as a full matrix.
%
h=ctools('cdh',x,v);
