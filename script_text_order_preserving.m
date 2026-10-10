% Section 5: fixed-threshold order preservation and the AB3 reference; 80 digits.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
[one, restore] = start_precision(80); %#ok<NASGU>
p = manuscript_problem('example3',one); N = [1024;2048];
methods = {'AB','AM','PC'}; E = zeros(2,4); A = zeros(2,3); residuals = zeros(2,1);
for j = 1:2
    for k = 1:3
        fprintf('80 digits: order-preserving %s, N=%d\n',methods{k},N(j));
        [t,v,s] = MA_adams(p,N(j),methods{k},'order-preserving',one/10,'exact');
        m = error_measures(p,t,v,s); E(j,k) = m.maximum; A(j,k) = s.activations;
        if k == 2, residuals(j) = s.max_equation_residual; end
    end
    [~,ta,va] = ab3(p.f,p.exact,p.a,p.b,p.u0,N(j));
    m = error_measures(p,ta.',va.'); E(j,4) = m.maximum;
end
P = observed_orders(N,E);
check_published('order_preserving',E(end,:),P(end,:),[],true);
opt = struct('maximum_error',true,'error_digits',4, ...
    'filename',fullfile(paths.tables,'text_order_preserving.tex'));
table_latex({'OP AB','OP AM','OP PC','AB3'},N,E,P,opt);
save_results(paths,'text_order_preserving',{'N','AB_Emax','AM_Emax','PC_Emax','AB3_Emax', ...
    'AB_p','AM_p','PC_p','AB3_p','AB_A','AM_A','PC_A','AM_equation_residual'}, ...
    [N,E,P,A,residuals],struct('problem','example3','digits',80,'delta',0.1, ...
    'placement','order-preserving','startup','exact v1,v2','AM_root','cancellation-free quadratic'));
end
