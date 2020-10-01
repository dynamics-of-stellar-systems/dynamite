function [varargout] = cuter_hprod( varargin )
% Return the matrix-vector product between the Hessian of the
% Lagrangian (or of the objective if problem is unconstrained) and a
% given vector p
% Usage:  r = cuter_hprod( x, v, p )   (Re)computes the Hessian at (x,v)
%         r = cuter_hprod( x, p )      Same, for unconstrained problems
%         r = cuter_hprod( p )         assumes H(x,v) was computed previously
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('hprod',varargin{:});
