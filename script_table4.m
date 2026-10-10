% Table 4: Nondegenerate control, both thresholds and unsafeguarded equivalence.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('control'); N = [16;32;64;128;256];
methods = {'AB','PC'}; E = zeros(numel(N),2); counts = zeros(numel(N),1);
trajectory_difference = zeros(numel(N),4);
for j = 1:numel(N)
    for k = 1:2
        [t0,v0,s0] = MA_adams(p,N(j),methods{k},'none',0,'exact');
        assert(s0.completed);
        for d = 1:2
            delta = d/10;
            [t,v,s] = MA_adams(p,N(j),methods{k},'original',delta,'exact');
            assert(s.completed && s.activations == 0,'Control safeguard unexpectedly active.');
            counts(j) = s.activations;
            diff = max(abs(v(:)-v0(:)));
            trajectory_difference(j,2*(k-1)+d) = diff;
            assert(diff < 5e-13,'Inactive safeguard and unsafeguarded trajectory disagree.');
            m = error_measures(p,t,v,s);
            if d == 1, E(j,k) = m.maximum; end
        end
    end
end
P = observed_orders(N,E); check_published('table4',E,P,counts);
opt = struct('fragment',true,'maximum_error',true,'shared_count',counts);
a = table_latex({'MA-RBF AB2','MA-RBF AB2--AM1'},N,E,P,opt);
write_paper_table(fullfile(paths.tables,'table4.tex'),{a},{''}, ...
    ['Nondegenerate control for Example~3 on $[1/4,1/2]$. ', ...
     'For $\delta_f=0.1,0.2$, the safeguard is inactive and the safeguarded ', ...
     'and unsafeguarded results coincide to roundoff.'],'table4');
settings = struct('problem','control','precision','double','startup','exact v1,v2', ...
    'interval',[1/4,1/2],'step','h=1/(4N)','deltas',[0.1,0.2]);
save_results(paths,'table4',{'N','AB_Emax','AB_p','PC_Emax','PC_p','A', ...
    'AB_delta01_difference','AB_delta02_difference','PC_delta01_difference','PC_delta02_difference'}, ...
    [N,E(:,1),P(:,1),E(:,2),P(:,2),counts,trajectory_difference],settings);
end
