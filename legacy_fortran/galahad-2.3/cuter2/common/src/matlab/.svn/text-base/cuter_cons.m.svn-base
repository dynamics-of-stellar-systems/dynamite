function [varargout] = cuter_cons( varargin )
% Return constraint bodies and Jacobian if requested.
% or return a single constraint value and its gradient if requested
% Usage:  c = cuter_cons(x)    or  [c,J]   = cuter_cons(x)
%        ci = cuter_cons(x,i)  or  [ci,gi] = cuter_cons(x,i)
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('cons',varargin{:});
