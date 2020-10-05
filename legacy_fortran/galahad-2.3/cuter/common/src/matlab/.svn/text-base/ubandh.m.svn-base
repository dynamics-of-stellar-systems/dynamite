function bandh=ubandh(x,nsemib,goth)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UBANDH Get elements of Hessian that lie within given semi-bandwidth.
%
%   bandh=ubandh(x,nsemib) evaluates the Hessian at x and returns in bandh
%   the Hessian elements that lie within nsemib of its diagonal.
%   bandh(1,i) contains the diagonal entry in column i of the Hessian.
%   bandh(j+1,i) contains the entry j places below diagonal in column i
%   of the Hessian.
%
%   bandh=ubandh(x,nsemib,goth), with goth=1, skips the evaluation of the
%   Hessian.  If goth is not given, it is assumed to be 0.
%
%   Set goth to 1 if the Hessian has already been evaluated by a call to
%   udh, ush, ugrdh or ugrsh at the current point, or if a previous call,
%   with goth set to 0, has been made to uprod or ubandh at the current point.
%
if nargin == 2
   bandh=utools('ubandh',x,nsemib);
else
   bandh=utools('ubandh',x,nsemib,goth);
end
