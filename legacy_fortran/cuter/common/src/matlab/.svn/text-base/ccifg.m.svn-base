function [ci,gci]=ccifg(x,i)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CCIFG Evaluates a particular constraint function and possibly
%       its gradient at x.
%
%   [ci,gci]=ccifg(x,i) returns the value of constraint number i at x
%                       in ci and its gradient in gci.
%
%   ci=ccifg(x,i)  returns the value of constraint number i at x in ci.
%
if nargout == 1 & nargin == 2
     ci=ctools('ccifg',x,i);
elseif nargout == 2 & nargin == 2
     [ci,gci]=ctools('ccifg',x,i);
end
