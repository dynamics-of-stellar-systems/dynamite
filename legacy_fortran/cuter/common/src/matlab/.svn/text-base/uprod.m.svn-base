function q=uprod(x,p,goth)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UPROD Evaluate product of vector p with Hessian at x.
%
%   uprod(x,p) evaluates the Hessian at x and returns the product of the 
%   Hessian with the vector p.  
%
%   uprod(x,p,goth), with goth=1, skips the evaluation of the Hessian.
%   If goth is not given, it is assumed to be 0.
%
%   Set goth to 1 if the Hessian has already been evaluated by a call to
%   udh, ush, ugrdh or ugrsh at the current point, or if a previous call,
%   with goth set to 0, has been made to uprod or ubandh at the current point.
%
if nargin == 2
   q=utools('uprod',x,p);
else
   q=utools('uprod',x,p,goth);
end
