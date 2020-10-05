function Hi=cidh(x,i)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CIDH Evaluate individual problem function Hessian x.
%   The value of i specifies the function:
%    i=0 -> objective function,
%    i>0 -> i-th constraint function.
%
%   Hi=cidh(x,i) returns the Hessian of the i-th function
%                in Hi, stored as a full matrix.
%
Hi=ctools('cidh',x,i);
