function [pname,xnames]=unames(n)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UNAMES Get problem name and variable names.
%
%   [pname,xnames]=unames returns the problem name in pname and
%   all variable names in xnames.  The variable names are stored
%   row-wise in xnames.
%
%   [pname,xnames]=unames(n) gets only the first n variable names.
%
if nargin == 0
   [pname,xnames]=utools('unames');
else
   [pname,xnames]=utools('unames',n);
end
xnames=setstr(xnames);
