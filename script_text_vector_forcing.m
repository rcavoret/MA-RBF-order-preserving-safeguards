% Section 5: vector logarithmic coefficient and tangent-projection cancellation.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
N = (512*2.^(0:6))'; values = []; problems = {'vector_log','vector_cancel'};
for ip = 1:2
    previous = NaN(2,2);
    for j = 1:numel(N)
        row = [N(j),ip]; current = zeros(2,2);
        methods = {'AB','AM'};
        for k = 1:2
            o = forcing_residuals(N(j),problems{ip},methods{k},'order-preserving',1/N(j));
            current(:,k) = o.zero;
            D = (o.zero-previous(:,k))/log(2);
            row = [row,o.zero.',D.',o.endpoint.',o.maximum]; %#ok<AGROW>
            if j == numel(N) && k == 1
                if ip == 1
                    assert(max(abs(D-[3.415094;-0.001573])) < 5.01e-7);
                else
                    assert(abs(o.maximum-3.1680) < 5.01e-5);
                end
            end
        end
        previous = current;
        values = [values;row]; %#ok<AGROW>
    end
end
headers = {'N','problem_1log_2cancel','AB_Szero_1','AB_Szero_2','AB_D_1','AB_D_2', ...
    'AB_endpoint_scaled_1','AB_endpoint_scaled_2','AB_maximum_over_h3', ...
    'AM_PC_Szero_1','AM_PC_Szero_2','AM_PC_D_1','AM_PC_D_2', ...
    'AM_PC_endpoint_scaled_1','AM_PC_endpoint_scaled_2','AM_PC_maximum_over_h3'};
save_results(paths,'text_vector_forcing',headers,values, ...
    struct('precision','double with compensated scaled-residual sums','delta','h', ...
    'norm','Euclidean','startup','exact v1,v2','reciprocal','radial', ...
    'log_forcing','(t+t^2,t^2)','cancellation_forcing','(t,t^2)'));
end
