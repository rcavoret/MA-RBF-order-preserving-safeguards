% Independent 80-digit MATLAB checks of Tables 3/4, coupled mantissas and vector identities.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
[one, restore] = start_precision(80); %#ok<NASGU>
p = manuscript_problem('example3',one); N = (128*2.^(0:6))';
methods = {'AB','PC'}; values = [];
for k = 1:2
    E = zeros(numel(N),2); A = E;
    for j = 1:numel(N)
        for d = 1:2
            fprintf('80 digits: Table 3, %s, N=%d, delta=%.1f\n',methods{k},N(j),d/10);
            [t,v,s] = MA_adams(p,N(j),methods{k},'original',d*one/10,'exact');
            m = error_measures(p,t,v,s); E(j,d) = m.maximum; A(j,d) = s.activations;
        end
    end
    P = observed_orders(N,E);
    if k == 1, name = 'table3a'; else, name = 'table3b'; end
    check_published(name,E,P,A,true);
    values = [values;N,repmat(k,numel(N),1),E,P,A]; %#ok<AGROW>
end
save_results(paths,'verification_80digits_table3',{'N','method_1AB_2PC','delta01_Emax','delta02_Emax', ...
    'delta01_p','delta02_p','delta01_A','delta02_A'},values,struct('digits',80,'table',3));
p = manuscript_problem('control',one); N = [16;32;64;128;256]; E = zeros(5,2);
for j = 1:numel(N)
    for k = 1:2
        [t,v,s] = MA_adams(p,N(j),methods{k},'original',one/10,'exact');
        m = error_measures(p,t,v,s); E(j,k) = m.maximum;
        assert(s.activations == 0);
    end
end
P = observed_orders(N,E); check_published('table4',E,P,zeros(5,1),true);
save_results(paths,'verification_80digits_table4',{'N','AB_Emax','PC_Emax','AB_p','PC_p'}, ...
    [N,E,P],struct('digits',80,'table',4));
% The paper does not specify 80 digits for the coupled example. These high-
% precision values expose two roundoff-sensitive final mantissa digits.
p = manuscript_problem('coupled',one); N = [512;1024]; E = zeros(2,3);
methods = {'AB','AM','PC'};
for j = 1:2
    for k = 1:3
        [t,v,s] = MA_adams(p,N(j),methods{k},'order-preserving',one/10,'exact');
        m = error_measures(p,t,v,s); E(j,k) = m.maximum;
    end
end
P = observed_orders(N,E); check_published('coupled_80',E(end,:),P(end,:),[],true);
save_results(paths,'verification_80digits_coupled',{'N','AB_Emax','AM_Emax','PC_Emax','AB_p','AM_p','PC_p'}, ...
    [N,E,P],struct('digits',80,'note','See NUMERICAL_NOTES.md for two final-digit differences'));
% Direct 80-digit vector trajectories versus the stable residual evaluator.
values = []; N = 512;
for ip = 1:2
    names = {'vector_log','vector_cancel'}; p = manuscript_problem(names{ip},one);
    for k = 1:3
        [t,v,s] = MA_adams(p,N,methods{k},'order-preserving',one/N,'exact');
        direct = double((p.exact(t)-v)*N^3);
        res = forcing_residuals(N,names{ip},methods{k},'order-preserving',1/N);
        gap = max(abs(direct(:)-res.scaled_errors(:)));
        assert(s.completed && gap < 5e-10);
        values = [values;ip,k,gap]; %#ok<AGROW>
    end
end
save_results(paths,'verification_80digits_vector_identity', ...
    {'problem_1log_2cancel','method_1AB_2AM_3PC','max_scaled_identity_gap'},values, ...
    struct('digits',80,'N',512,'comparison','direct recurrence vs compensated scaled-residual sums'));
fprintf('All requested 80-digit cross-checks passed.\n');
end
