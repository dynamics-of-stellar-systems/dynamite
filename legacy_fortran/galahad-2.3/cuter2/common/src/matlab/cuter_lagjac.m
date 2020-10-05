function [varargout] = cuter_lagjac( varargin )
% Return the gradient of the objective or Lagrangian and Jacobian
% [g,J] = cuter_lagjac(x)  or  [g,J] = cuter_lagjac(x,v)
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('lagjac',varargin{:});
