function [varargout] = cuter_objcons( varargin )
% Evaluate objective function value and constraint bodies.
% Usage:  [f,c] = cuter_objcons(x)
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('objcons',varargin{:});
