% Table 1(a,b): Example 1, Euler and standalone AM1.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example1');
N = (10*2.^(0:7))';
E1 = zeros(numel(N),3); E2 = E1;
iterations = zeros(numel(N),3); residuals = iterations;
for j = 1:numel(N)
    E1(j,1) = euler(p.f,p.exact,p.a,p.b,p.u0,N(j));
    E1(j,2) = RBF_euler(p.f,p.exact,p.a,p.b,p.u0,N(j),1);
    E1(j,3) = RBF_euler(p.f,p.exact,p.a,p.b,p.u0,N(j),2);
    kernels = {'classical','MQ','MA'};
    for k = 1:3
        [E2(j,k),~,~,s] = RBF_am1(p.f,p.jacobian,p.exact,p.a,p.b,p.u0,N(j),kernels{k});
        iterations(j,k) = s.max_newton_iterations;
        residuals(j,k) = s.max_equation_residual;
    end
end
P1 = observed_orders(N,E1); P2 = observed_orders(N,E2);
check_published('table1a',E1,P1);
check_published('table1b',E2,P2);
assert(isequal(E2(:,2),E2(:,3)),'MQ/MA AM1 common-weight identity failed.');
assert(max(iterations(:)) <= 4,'A displayed Newton solve needed over four iterations.');
opt.fragment = true;
a = table_latex({'Euler','MQ-RBF Euler','MA-RBF Euler'},N,E1,P1,opt);
b = table_latex({'AM1','MQ-RBF AM1','MA-RBF AM1'},N,E2,P2,opt);
write_paper_table(fullfile(paths.tables,'table1.tex'),{a,b}, ...
    {'(a) Euler formulas','(b) Standalone Adams--Moulton formulas'}, ...
    ['Endpoint errors and observed rates for Example~1 at $b=1$: ', ...
     '(a) Euler and (b) standalone AM1 formulas. Rates use unrounded errors.'],'table1');
settings = struct('problem','example1','precision','double', ...
    'startup','Euler: initial value only; RBF Euler: exact v1; all AM1: exact v1,v2', ...
    'newton_increment_tolerance',1e-14,'newton_guess','forward Euler');
save_results(paths,'table1a',{'N','Euler_E','Euler_p','MQ_E','MQ_p','MA_E','MA_p'}, ...
    [N,E1(:,1),P1(:,1),E1(:,2),P1(:,2),E1(:,3),P1(:,3)],settings);
save_results(paths,'table1b',{'N','AM1_E','AM1_p','MQ_E','MQ_p','MA_E','MA_p', ...
    'AM1_Newton_max','MQ_Newton_max','MA_Newton_max','AM1_residual','MQ_residual','MA_residual'}, ...
    [N,E2(:,1),P2(:,1),E2(:,2),P2(:,2),E2(:,3),P2(:,3),iterations,residuals],settings);
end
