% Section 5, Example 1: same-order AB2 reference against MA-RBF Euler.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example1'); N = 1280;
E = [ab2(p.f,p.exact,p.a,p.b,p.u0,N), ...
     RBF_euler(p.f,p.exact,p.a,p.b,p.u0,N,2)];
assert(abs(E(1)-1.91e-7) <= 0.005e-7);
assert(abs(E(2)-2.67e-7) <= 0.005e-7);
opt = struct('filename',fullfile(paths.tables,'text_same_order.tex'));
table_latex({'AB2','MA-RBF Euler'},N,E,[NaN,NaN],opt);
save_results(paths,'text_same_order',{'N','AB2_endpoint','MA_Euler_endpoint'},[N,E], ...
    struct('problem','example1','precision','double','startup','common exact v1'));
end
