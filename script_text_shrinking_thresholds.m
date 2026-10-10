% Section 5 / powersharp: shrinking thresholds alpha=1/4,1/2,3/4, both placements.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
% The paper specifies alpha, but no mesh list for this supplementary test.
% The fully fixed reproduction grid here is N=2^7,...,2^20.
N = (2.^(7:20))'; alphas = [1/4,1/2,3/4]; values = [];
for a = alphas
    e = zeros(numel(N),2); emax = e; scaled = e; opmax = e;
    for j = 1:numel(N)
        h = 1/N(j); delta = h^a;
        methods = {'AB','AM'};
        for k = 1:2
            o = forcing_residuals(N(j),'forcing',methods{k},'original',delta);
            e(j,k) = h^3*o.endpoint;
            emax(j,k) = h^3*o.maximum;
            scaled(j,k) = h^(1-a)*o.endpoint;
            op = forcing_residuals(N(j),'forcing',methods{k},'order-preserving',delta);
            opmax(j,k) = op.maximum;
        end
    end
    rates = observed_orders(N,abs(e));
    values = [values;N,repmat(a,numel(N),1),e,rates,scaled,emax,opmax]; %#ok<AGROW>
end
save_results(paths,'text_shrinking_thresholds',{'N','alpha','original_endpoint_AB','original_endpoint_AM_PC', ...
    'endpoint_rate_AB','endpoint_rate_AM_PC','endpoint_over_h2alpha_AB','endpoint_over_h2alpha_AM_PC', ...
    'original_maximum_AB','original_maximum_AM_PC','OP_maximum_over_h3_AB','OP_maximum_over_h3_AM_PC'}, ...
    values,struct('problem','F=t+t^2','precision','double with compensated scaled-residual sums', ...
    'delta','h^alpha','predicted_endpoint_constants',[10/9,-2/9], ...
    'startup','exact v1,v2','mesh_choice','N=2^7,...,2^20, explicitly supplied by this package'));
fprintf('Endpoint ratios tend to 10/9 and -2/9; convergence is slow, especially for alpha=3/4.\n');
end
