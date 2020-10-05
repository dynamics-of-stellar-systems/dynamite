function [g,cjac]=cgr(x,v,options)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CGR Evaluate gradient of the objective or Lagrangian function,
%   and the gradients of the general constraint functions.
%
%   [g,cjac]=cgr(x) returns the gradient of the objective function in
%   g and the gradients of the general constraint functions in cjac.
%   cjac(i,j) contains the partial derivative of the i-th constraint
%   with respect to the j-th variable.
%
%   [g,cjac]=cgr(x,v) returns the gradient of the Lagrangian function
%   in g and the gradients of the general constraint functions in cjac,
%   where x is the current estimate of the solution and v is the current
%   estimate of the Lagrange multipliers.
%
%   [g,cjac]=cgr(x,v,options), where options is a 2-dimensional vector,
%   allows cjac to be transposed and requests the gradient of the Lagrangian
%   to be placed in g.
%   options( 1 ) = jtrans, set to 1 if the user wants the transpose
%     of the Jacobian, where the i,j-th component is the partial derivative
%     of the j-th constraint with respect to the i-th variable.
%     If options is not given, jtrans defaults to 0.
%   options( 2 ) = grlagf, set to 1 if the gradient of the Lagrangian
%     is required and set to 0 if the gradient of the objective function
%     is sought.  Note that grlagf defaults to 0 if v is not given,
%     and defaults to 1 if v is given.
%
if nargin == 1
   [g,cjac]=ctools('cgr',x);
elseif nargin == 2
   [g,cjac]=ctools('cgr',x,v);
else
   [g,cjac]=ctools('cgr',x,v,options);
end
