function [varargout] = cuter_jprod( varargin )
% Return the product of the Jacobian at x with a vector p.
% Usage:  r = cuter_jprod( x, p )  --> recomputes J(x)
%         r = cuter_jprod( p )     --> assumes J(x) was computed previously
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('Jprod',varargin{:});
