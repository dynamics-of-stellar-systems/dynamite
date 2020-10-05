function [x,bl,bu,v,cl,cu,equatn,linear]=csetup(options)
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% CSETUP Set up the data structures for subsequent computations.
%   CSETUP must be called before any other constrained tool.
%
%   [x,bl,bu,v,cl,cu]=csetup returns the starting point in x and the
%   lower and upper bounds on the variables in bl and bu, respectively.
%   v contains the initial estimates of the Lagrange multipliers,
%   and cl and cu contain the lower and upper bounds on the constraints,
%   respectively.
%
%   [x,bl,bu,v,cl,cu,equatn,linear]=csetup also returns equatn, an array
%   whose i-th component is 1 if the i-th constraint is an equation and
%   0 otherwise, and linear, an array whose i-th component is 1 if the
%   i-th constraint is linear and 0 otherwise.
%
%   [x,bl,bu,v,cl,cu]=csetup(options) or 
%   [x,bl,bu,v,cl,cu,equatn,linear]=csetup(options),
%   where options is a 3-dimensional vector,
%   permits the reordering of the constraints and variables.
%   options( 1 ) = efirst, set to 1 if the user wishes the general
%     equations to occur before the general inequalities.
%   options( 2 ) = lfirst, set to 1 if the user wishes the general
%     linear or affine constraints to occur before the general nonlinear
%     ones.  If both efirst and lfirst are set to 1, the linear constraints
%     will occur before the nonlinear ones.  The linear constraints are
%     ordered so that the linear equations occur before the linear
%     inequalities.  Likewise, the nonlinear equations appear before the
%     before the nonlinear inequalities.
%   options( 3 ) = nvfrst, set to 1 if the user wishes the nonlinear
%     variables to occur before the linear variables.  Any variable which
%     belongs to a nontrivial group or to a nonlinear element in a trivial
%     group is treated as a nonlinear variable.  If the number of variables
%     which appear nonlinearly in the objective function (say n_1) is
%     different from the number of variables which appear nonlinearly in
%     the constraints (say m_1), then the variables are ordered so that
%     the smaller set occurs first.  For example, if n_1 < m_1, the n_1
%     nonlinear objective variables occur first, followed by the nonlinear
%     Jacobian variables not belonging to the first n_1 variables, followed
%     by the linear variables.
%   If options is not given, the constraints and variables are not reordered.
%
if nargin == 0 & nargout == 6
   [x,bl,bu,v,cl,cu]=ctools('csetup');
elseif nargin == 1 & nargout == 6
   [x,bl,bu,v,cl,cu]=ctools('csetup',options);
elseif nargin == 0 & nargout == 8
   [x,bl,bu,v,cl,cu,equatn,linear]=ctools('csetup');
else
   [x,bl,bu,v,cl,cu,equatn,linear]=ctools('csetup',options);
end
