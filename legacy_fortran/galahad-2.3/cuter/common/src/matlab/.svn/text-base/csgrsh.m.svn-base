function [g,cjac,h]=csgrsh(x,v,options)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CSGRSH Evaluate gradient of objective or Lagrangian function, gradient
%   of general constraint functions, and Hessian of Lagrangian.
%
%   [g,cjac,h]=csgrsh(x,v) returns the sparse gradient of the objective
%   function in g, the sparse constraint Jacobian matrix in cjac, and the
%   sparse Hessian of the Lagrangian in h.
%
%   The i,j-th nonzero entry in cjac corresponds to the partial derivative
%   of the i-th constraint with respect to the j-th variable.
%
%   The i,j-th nonzero entry in h corresponds to the partial derivative
%   of the Lagrangian function with respect to the i-th and j-th variables.
%   Since the Hessian is symmetric, only the upper triangular entries are
%   returned.
%
%   [g,cjac,h]=csgrsh(x,v,options), where options is a 2-dimensional vector,
%   allows cjac to be transposed and requests the gradient of the Lagrangian
%   to be placed in g.
%   options( 1 ) = jtrans, set to 1 if the user wants the transpose
%     of the Jacobian, where the i,j-th nonzero entry is the partial derivative
%     of the j-th constraint with respect to the i-th variable.
%     If options is not given, jtrans defaults to 0.
%   options( 2 ) = grlagf, set to 1 if the gradient of the Lagrangian
%     is required and set to 0 if the gradient of the objective function
%     is sought.  If options is not given, grlagf defaults to 0.
%
if nargin == 2
   [g,cjac,h]=ctools('csgrsh',x,v);
else
   [g,cjac,h]=ctools('csgrsh',x,v,options);
   if options(1) ~= 0
      cjac=cjac';
   end
end
