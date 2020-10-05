function [varargout] = cuter_scons( varargin )
% Return constraint bodies and sparse Jacobian
% or return a single constraint value and its gradient in sparse format
% Usage: [c,J] = cuter_scons(x)
%        [ci, sgci] = cuter_scons( prob.x, i )
    varargout = cell(1,nargout);
    [varargout{:}] = mcuter('scons',varargin{:});
