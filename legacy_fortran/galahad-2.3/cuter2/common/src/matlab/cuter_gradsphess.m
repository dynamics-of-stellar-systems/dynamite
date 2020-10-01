function [varargout] = cuter_gradsphess( varargin )
% Return the Hessian of the Lagrangian, the Jacobian of the constraints
% and the gradient of either the objective function or the Lagrangian
% Usage:  [g,H] = cuter_gradsphess( x )   if the problem is unconstrained, or
%       [g,J,H] = cuter_gradsphess( x, v, gradf )  if it is constrained
%                 where gradf is either true or false;
%                  gradf = true : returns the gradient of the objective in g
%                          false: returns the gradient of the Lagrangian in g
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('gradsphess',varargin{:});
