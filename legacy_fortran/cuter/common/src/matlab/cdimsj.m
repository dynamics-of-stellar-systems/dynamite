function nnzj=cdimsj
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CDIMSJ gets the number of nonzeros to store the Jacobian of the constraints
%
%   nnzj=cdimsj returns the number of nonzero elements required to store
%              the Jacobian matrix of the constraint functions.
%
nnzj=ctools('cdimsj');
