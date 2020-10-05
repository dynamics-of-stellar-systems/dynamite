function [varargout] = cuter_jtprod( varargin )
% Return the product of the transpose Jacobian at x with a vector p.
% Usage:  r = cuter_jtprod( x, p )  --> recomputes J(x)
%         r = cuter_jtprod( p )     --> assumes J(x) was computed previously
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('Jtprod',varargin{:});
