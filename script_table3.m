% Table 3(a,b): Original safeguards, maximum errors, activations and even-mesh failure.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example3');
N = (128*2.^(0:6))'; deltas = [0.1,0.2]; methods = {'AB','PC'};
E = zeros(numel(N),2,2); A = E; failed_at = zeros(numel(N),2);
for j = 1:numel(N)
    for k = 1:2
        % An incomplete trajectory must NOT receive a full-interval error.
        [tu,vu,su] = MA_adams(p,N(j),methods{k},'none',0,'exact');
        mu = error_measures(p,tu,vu,su);
        assert(~su.completed && su.stop_n == N(j)/2 && isnan(mu.maximum));
        failed_at(j,k) = su.stop_n;
        for d = 1:2
            [t,v,s] = MA_adams(p,N(j),methods{k},'original',deltas(d),'exact');
            assert(s.completed,'A safeguarded trajectory failed.');
            m = error_measures(p,t,v,s);
            E(j,d,k) = m.maximum; A(j,d,k) = s.activations;
        end
    end
end
P1 = observed_orders(N,E(:,:,1)); P2 = observed_orders(N,E(:,:,2));
check_published('table3a',E(:,:,1),P1,A(:,:,1));
check_published('table3b',E(:,:,2),P2,A(:,:,2));
opt = struct('fragment',true,'maximum_error',true,'counts',A(:,:,1));
a = table_latex({'$\delta_f=0.1$','$\delta_f=0.2$'},N,E(:,:,1),P1,opt);
opt.counts = A(:,:,2);
b = table_latex({'$\delta_f=0.1$','$\delta_f=0.2$'},N,E(:,:,2),P2,opt);
write_paper_table(fullfile(paths.tables,'table3.tex'),{a,b}, ...
    {'(a) MA-RBF AB2','(b) MA-RBF AB2--AM1 predictor--corrector'}, ...
    ['Safeguard activation in Example~3 on $[-1/2,1/2]$ ', ...
     'for (a) MA-RBF AB2 and (b) MA-RBF AB2--AM1. ', ...
     'Here $\mathcal E_N$ is the maximum nodal error and $A_N$ the number of active updates. ', ...
     'Both unsafeguarded schemes encounter a zero denominator on every listed mesh.'],'table3');
settings = struct('problem','example3','precision','double','startup','exact v1,v2', ...
    'placement','original','deltas',deltas,'activation_indices','n=1,...,N-2');
headers = {'N','delta01_Emax','delta01_p','delta01_A','delta02_Emax','delta02_p','delta02_A','unsafeguarded_stop_n'};
save_results(paths,'table3a',headers,[N,E(:,1,1),P1(:,1),A(:,1,1),E(:,2,1),P1(:,2),A(:,2,1),failed_at(:,1)],settings);
save_results(paths,'table3b',headers,[N,E(:,1,2),P2(:,1),A(:,1,2),E(:,2,2),P2(:,2),A(:,2,2),failed_at(:,2)],settings);
end
