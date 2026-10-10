% Section 5: sharp logarithmic coefficients, endpoint cancellation, direct 80-digit identity.
% Run this script without editing parameters.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
[one, restore] = start_precision(80); %#ok<NASGU>
N = (512*2.^(0:8))';
Szero = zeros(numel(N),2); Send = Szero; OriginalEnd = Szero; D = NaN(size(Szero));
previous = [];
for j = 1:numel(N)
    sums = scalar_sharp_sums(N(j),one);
    Szero(j,:) = double(sums.zero); Send(j,:) = double(sums.endpoint);
    OriginalEnd(j,:) = double(sums.original_endpoint);
    if j > 1
        D(j,:) = double((sums.zero-previous)/log(2*one));
    end
    previous = sums.zero;
end
assert(max(abs(D(end,:)-[3.416235,-0.249950])) < 5.01e-7);
assert(max(abs(Send(end,:)-[-3.753485,0.274644])) < 5.01e-7);
% Independent evaluation of the actual computed maps, not just the residuals.
p = manuscript_problem('forcing',one); M = 512; theory = scalar_sharp_sums(M,one);
methods = {'AB','AM','PC'}; placements = {'original','order-preserving'};
identity = zeros(6,4); row = 0;
for g = 1:2
    for k = 1:3
        [t,v,s] = MA_adams(p,M,methods{k},placements{g},one/M,'exact');
        assert(s.completed);
        err = (p.exact(t)-v)*M^3;
        refzero = theory.zero(min(k,2));
        if g == 1, refend = theory.original_endpoint(min(k,2));
        else, refend = theory.endpoint(min(k,2)); end
        gap0 = abs(double(err(M/2+1)-refzero));
        gap1 = abs(double(err(end)-refend));
        assert(max(gap0,gap1) < 1e-60,'Direct recurrence/residual identity failed.');
        row = row+1; identity(row,:) = [g,k,gap0,gap1];
    end
end
settings = struct('problem','F=t+t^2','digits',80,'delta','h', ...
    'startup','exact v1,v2','signed_error','exact minus computed', ...
    'summation','exact residual sums evaluated through harmonic-number identities');
save_results(paths,'text_sharp_asymptotics',{'N','Szero_AB','Szero_AM_PC','D_AB','D_AM_PC', ...
    'OP_endpoint_scaled_AB','OP_endpoint_scaled_AM_PC','original_endpoint_scaled_AB','original_endpoint_scaled_AM_PC'}, ...
    [N,Szero,D,Send,OriginalEnd],settings);
save_results(paths,'text_sharp_identity',{'placement_1original_2OP','method_1AB_2AM_3PC', ...
    'zero_scaled_difference','endpoint_scaled_difference'},identity,settings);
text = sprintf(['%% Computed at 80 digits, N=%d.\n', ...
    '$D_N^{AB}=%.6f$, $D_N^{AM/PC}=%.6f$.\n', ...
    '$N^3 e_N^{AB}(1/2)=%.6f$, $N^3 e_N^{AM/PC}(1/2)=%.6f$.\n'], ...
    N(end),D(end,1),D(end,2),Send(end,1),Send(end,2));
write_text_file(fullfile(paths.tables,'text_sharp_asymptotics.tex'),text);
end
