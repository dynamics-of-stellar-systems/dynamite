function [varargout] = cuter_slagjac( varargin )
% Return the sparse Jacobian and gradient of either the objective
% function or the Lagrangian.
% Usage:  [J,g] = cuter_slagjac(x) or [J,g] = cuter_slagjac(x,v)
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('slagjac',varargin{:});
