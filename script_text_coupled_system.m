% Section 5: genuinely coupled system and orthogonal/translation equivariance; double precision.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('coupled'); N = [512;1024]; methods = {'AB','AM','PC'};
E = zeros(2,3); gaps = E; residuals = zeros(2,1);
U = [cos(0.7),-sin(0.7);sin(0.7),cos(0.7)]; shift = [0.3;-0.4];
rotated = transform_problem(p,U,shift);
for j = 1:2
    for k = 1:3
        [t,v,s] = MA_adams(p,N(j),methods{k},'order-preserving',0.1,'exact');
        m = error_measures(p,t,v,s); E(j,k) = m.maximum;
        if k == 2, residuals(j) = s.max_equation_residual; end
        [~,w,sw] = MA_adams(rotated,N(j),methods{k},'order-preserving',0.1,'exact');
        assert(sw.completed && s.completed);
        difference = w-bsxfun(@plus,U*v,shift);
        gaps(j,k) = max(sqrt(sum(difference.^2,1)));
        assert(gaps(j,k) < 1e-12,'Equivariance check exceeds the roundoff tolerance.');
    end
end
P = observed_orders(N,E);
% These manuscript mantissas are double-roundoff sensitive.
check_published('coupled',E(end,:),P(end,:),[],false);
if max(gaps(:)) >= 3e-14
    warning('MARBF:PlatformRoundoff', ...
        'Equivariance gap %.3g exceeds the paper''s machine-specific 3e-14, but is below 1e-12.',max(gaps(:)));
end
opt = struct('maximum_error',true,'error_digits',4, ...
    'filename',fullfile(paths.tables,'text_coupled_system.tex'));
table_latex({'OP AB','OP AM','OP PC'},N,E,P,opt);
save_results(paths,'text_coupled_system',{'N','AB_Emax','AM_Emax','PC_Emax','AB_p','AM_p','PC_p', ...
    'AB_equivariance_gap','AM_equivariance_gap','PC_equivariance_gap','AM_equation_residual'}, ...
    [N,E,P,gaps,residuals],struct('problem','coupled','precision','double','delta',0.1, ...
    'startup','exact v1,v2','norm','Euclidean','rotation',0.7,'translation',shift, ...
    'AM_solver','direct 2-by-2 linear solve'));
end
