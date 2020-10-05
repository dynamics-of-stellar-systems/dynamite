% Test script to exercise the Matlab/CUTEr interface.
% D. Orban, 2011.

clear all;

fprintf('Reading problem dimensions...\n');
[nvar,ncon] = cuter_dims();

fprintf('Initializing problem...\n');
prob = cuter_setup();

assert(nvar == prob.n);
assert(ncon == prob.m);

fprintf('Reading variable names...\n');
vnames = cuter_varnames();

fprintf('Evaluating objective function...\n');
f = cuter_obj(prob.x);

fprintf('Evaluating objective and gradient...\n');
[f,g] = cuter_obj(prob.x);

if prob.m > 0
    fprintf('Reading constraint names...\n');
    cnames = cuter_connames();

    fprintf('Evaluating objective and constraints...\n');
    [f,c] = cuter_objcons(prob.x);

    fprintf('Evaluating dense Hessian of Lagrangian...\n');
    H = cuter_hess(prob.x,prob.v);

    fprintf('Evaluating gradient and Hessian of Lagrangian\n');
    fprintf('  and Jacobian of constraints...\n');
    [g,J,H] = cuter_gradhess(prob.x,prob.v,false,false);

    fprintf('Evaluating gradient and Hessian of Lagrangian\n');
    fprintf('  and transpose Jacobian of constraints...\n');
    [g,Jt,H] = cuter_gradhess(prob.x,prob.v,false,true);

    fprintf('Computing Hessian-vector product...\n');
    p = rand(prob.n,1);
    Hp = cuter_hprod(p);  % Assume Hessian was evaluated earlier.
    fprintf('  error = %7.1e\n', norm(Hp-H*p)/(1+norm(Hp)));

    fprintf('Computing dense Hessian of first constraint...\n');
    H1 = cuter_ihess(prob.x,1);

    %fprintf('Computing sparse Hessian of first constraint...\n');
    %H1sp = cuter_isphess(prob.x,1);

    fprintf('Computing Jacobian-vector product...\n');
    Jp = cuter_jprod(p);  % Assume J was evaluated earlier.
    fprintf('  error = %7.1e\n', norm(Jp-J*p)/(1+norm(Jp)));

    fprintf('Computing transpose Jacobian-vector product...\n');
    q = rand(prob.m,1);
    Jtq = cuter_jtprod(q);  % Assume J was evaluated earlier.
    fprintf('  error = %7.1e\n', norm(Jtq-Jt*q)/(1+norm(Jtq)));

    fprintf('Computing gradient of Lagrangian and dense constraint Jacobian...\n');
    [g,J] = cuter_lagjac(prob.x,prob.v);

    fprintf('Computing constraints and sparse Jacobian...\n');
    [c,J] = cuter_scons(prob.x);

    fprintf('Computing first constraint body and sparse gradient...\n');
    [c1,g1] = cuter_scons(prob.x,1);

    %fprintf('Computing gradient of Lagrangian and sparse constraint Jacobian...\n');
    %[sg,sJ] = cuter_slagjac(prob.x,prob.v);

    fprintf('Computing sparse Lagrangian Hessian...\n');
    Hsp = cuter_sphess(prob.x,prob.v);

else

    fprintf('Evaluating dense objective Hessian...\n');
    H = cuter_hess(prob.x);

    fprintf('Evaluating objective gradient and dense Hessian...\n');
    [g,H2] = cuter_gradhess(prob.x);

    fprintf('  error = %7.1e\n', norm(H-H2)/(1+norm(H)));

    fprintf('Computing Hessian-vector product...\n');
    p = rand(prob.n,1);
    Hp = cuter_hprod(p);  % Assume Hessian was evaluated earlier.
    fprintf('  error = %7.1e\n', norm(Hp-H*p)/(1+norm(Hp)));

    fprintf('Computing sparse objective Hessian...\n');
    Hsp = cuter_sphess(prob.x);

    fprintf('  error = %7.1e\n', norm(H-full(Hsp))/(1+norm(H)));
end
