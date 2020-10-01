function nnz=cdimsh
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CDIMSH gets the number of nonzeros to store the Hessian of the Lagrangian
%
%   nnz=cdimsh returns the number of nonzero elements required to store
%              the Hessian matrix of the Lagrangian.
%
nnz=ctools('cdimsh');
