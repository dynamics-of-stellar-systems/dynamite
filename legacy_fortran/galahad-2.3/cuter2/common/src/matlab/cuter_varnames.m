function [varargout] = cuter_varnames( varargin )
% Return variable names.
% Usage: vnames = cuter_varnames()
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('varnames',varargin{:});
