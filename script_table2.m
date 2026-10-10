% Table 2(a,b): Example 2, unsafeguarded AB2 and PC.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example2');
N = (10*2.^(0:5))'; E1 = zeros(numel(N),3); E2 = E1;
for j = 1:numel(N)
    E1(j,1) = ab2(p.f,p.exact,p.a,p.b,p.u0,N(j));
    E1(j,2) = RBF_adams(p.f,p.exact,p.a,p.b,p.u0,N(j),1);
    E1(j,3) = RBF_adams(p.f,p.exact,p.a,p.b,p.u0,N(j),3);
    E2(j,1) = ab2am1(p.f,p.exact,p.a,p.b,p.u0,N(j));
    E2(j,2) = RBF_adams(p.f,p.exact,p.a,p.b,p.u0,N(j),2);
    E2(j,3) = RBF_adams(p.f,p.exact,p.a,p.b,p.u0,N(j),4);
end
P1 = observed_orders(N,E1); P2 = observed_orders(N,E2);
check_published('table2a',E1,P1); check_published('table2b',E2,P2);
opt.fragment = true;
a = table_latex({'AB2','MQ-RBF AB2','MA-RBF AB2'},N,E1,P1,opt);
b = table_latex({'AB2--AM1','MQ-RBF AB2--AM1','MA-RBF AB2--AM1'},N,E2,P2,opt);
write_paper_table(fullfile(paths.tables,'table2.tex'),{a,b}, ...
    {'(a) Adams--Bashforth formulas','(b) Adams--Bashforth--Moulton predictor--corrector'}, ...
    ['Endpoint errors and observed rates for Example~2 at $b=3$: ', ...
     '(a) AB2 and (b) AB2--AM1 predictor--corrector formulas.'],'table2');
settings = struct('problem','example2','precision','double','safeguard','none', ...
    'startup','exact v1 for classical AB2/PC; exact v1,v2 for RBF methods');
headers = {'N','classical_E','classical_p','MQ_E','MQ_p','MA_E','MA_p'};
save_results(paths,'table2a',headers,[N,E1(:,1),P1(:,1),E1(:,2),P1(:,2),E1(:,3),P1(:,3)],settings);
save_results(paths,'table2b',headers,[N,E2(:,1),P2(:,1),E2(:,2),P2(:,2),E2(:,3),P2(:,3)],settings);
end
