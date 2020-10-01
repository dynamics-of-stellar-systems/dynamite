function [pname,xnames,gnames]=cnames(n,m)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CNAMES Get problem name, variable names, and constraint names.
%
%   [pname,xnames,gnames]=cnames returns the problem name in pname,
%   all variable names in xnames, and all constraint names in gnames.
%   The names are stored row-wise in xnames and gnames.
%
%   [pname,xnames,gnames]=cnames(n,m) gets only the first n variable names
%   and only the first m constraint names.
%
if nargin == 0
   [pname,xnames,gnames]=ctools('cnames');
else
   [pname,xnames,gnames]=ctools('cnames',n,m);
end
xnames=setstr(xnames);
gnames=setstr(gnames);
