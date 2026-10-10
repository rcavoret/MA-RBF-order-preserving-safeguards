% Section 5, Example 1: two RK4 starts versus exact starts for AB/AM/PC.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example1'); N = [320;640]; methods = {'AB','AM','PC'};
E = zeros(2,3); Exact = E;
for j = 1:2
    for k = 1:3
        [t,v,s] = MA_adams(p,N(j),methods{k},'none',0,'rk4');
        m = error_measures(p,t,v,s); E(j,k) = m.endpoint;
        [t,v,s] = MA_adams(p,N(j),methods{k},'none',0,'exact');
        m = error_measures(p,t,v,s); Exact(j,k) = m.endpoint;
    end
end
P = observed_orders(N,E); Pexact = observed_orders(N,Exact);
assert(all(abs(P(2,:)-[2.9994,2.9931,2.9811]) < 1e-4));
opt = struct('filename',fullfile(paths.tables,'text_startup.tex'));
table_latex({'MA-RBF AB2 (RK4)','MA-RBF AM1 (RK4)','MA-RBF PC (RK4)'},N,E,P,opt);
save_results(paths,'text_startup',{'N','AB_RK4_E','AM_RK4_E','PC_RK4_E', ...
    'AB_RK4_p','AM_RK4_p','PC_RK4_p','AB_exact_E','AM_exact_E','PC_exact_E', ...
    'AB_exact_p','AM_exact_p','PC_exact_p'},[N,E,P,Exact,Pexact], ...
    struct('problem','example1','precision','double','startup','two classical RK4 steps', ...
    'AM_Newton_increment_tolerance',1e-14));
end
