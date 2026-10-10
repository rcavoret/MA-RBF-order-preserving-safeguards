% Section 5, Example 3: even/odd mesh robustness; no error assigned to failed runs.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
p = manuscript_problem('example3'); evenN = (128*2.^(0:6))';
methods = {'AB','PC'}; values = zeros(numel(evenN)*4,8); row = 0;
for j = 1:numel(evenN)
    for odd = 0:1
        N = evenN(j)-odd;
        for k = 1:2
            [t,v,s] = MA_adams(p,N,methods{k},'none',0,'exact');
            m = error_measures(p,t,v,s);
            if odd
                assert(s.completed,'An odd companion run did not complete.');
            else
                assert(~s.completed && s.stop_n == N/2);
                assert(isnan(m.maximum) && isnan(m.endpoint));
            end
            row = row+1;
            values(row,:) = [N,k,s.completed,s.stop_n,s.completed_steps,m.maximum,m.endpoint,odd];
        end
    end
end
save_results(paths,'text_robustness',{'N','method_1AB_2PC','completed','stop_n', ...
    'completed_steps','maximum_error','endpoint_error','odd_mesh'},values, ...
    struct('problem','example3','precision','double','safeguard','none', ...
    'startup','exact v1,v2','NaN_error','incomplete run; NOT an error estimate'));
fprintf('All 14 even-mesh runs stop at n=N/2; all 14 odd-mesh runs complete.\n');
end
