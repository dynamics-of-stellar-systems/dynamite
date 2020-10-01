function [varargout] = cuter_setup( varargin )
% Set up main problem structure
% Usage:  prob = cuter_setup()
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('setup',varargin{:});
