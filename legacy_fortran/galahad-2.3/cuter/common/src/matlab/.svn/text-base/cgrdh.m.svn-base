function [g,cjac,h]=cgrdh(x,v,options)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CGRDH Evaluate gradient of objective or Lagrangian function, gradients
%   of general constraint functions, and Hessian of Lagrangian.
%
%   [g,cjac,h]=cgrdh(x,v) returns the gradient of the objective function in
%   g, the gradients of the general constraint functions in cjac, and the
%   Hessian of the Lagrangian in h.  The Hessian is stored as a full matrix.
%   cjac(i,j) contains the partial derivative of the i-th constraint
%   with respect to the j-th variable.
%
%   [g,cjac,h]=cgrdh(x,v,options), where options is a 2-dimensional vector,
%   allows cjac to be transposed and requests the gradient of the Lagrangian
%   to be placed in g.
%   options( 1 ) = jtrans, set to 1 if the user wants the transpose
%     of the Jacobian, where the i,j-th component is the partial derivative
%     of the j-th constraint with respect to the i-th variable.
%     If options is not given, jtrans defaults to 0.
%   options( 2 ) = grlagf, set to 1 if the gradient of the Lagrangian
%     is required and set to 0 if the gradient of the objective function
%     is sought.  If options is not given, grlagf defaults to 0.
%
if nargin == 2
   [g,cjac,h]=ctools('cgrdh',x,v);
else
   [g,cjac,h]=ctools('cgrdh',x,v,options);
end
