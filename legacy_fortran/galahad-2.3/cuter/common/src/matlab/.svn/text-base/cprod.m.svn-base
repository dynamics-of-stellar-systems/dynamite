function q=cprod(x,v,p,goth)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CPROD Evaluate product of vector p with Hessian of Lagrangian at x and v.
%
%   cprod(x,v,p) evaluates the Hessian of the Lagrangian at x and v,
%   and returns the product of the Hessian with the vector p.
%
%   cprod(x,v,p,goth), with goth=1, skips the evaluation of the Hessian.
%   If goth is not given, it is assumed to be 0.
%
%   Set goth to 1 if the Hessian has already been evaluated by a call to
%   cdh, csh, cgrdh or csgrsh at the current point, or if a previous call,
%   with goth set to 0, has been made to cprod at the current point.
%
if nargin == 3
   q=ctools('cprod',x,v,p);
else
   q=ctools('cprod',x,v,p,goth);
end
