% Reproduce ALL numbered/textual results and 80-digit checks; requires Symbolic Math Toolbox.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
% Fail clearly BEFORE the full run if the explicitly required toolbox is absent.
[~, restore] = start_precision(80); clear restore
check_invariants();
names = {'script_table1.m','script_table2.m','script_table3.m','script_table4.m', ...
    'script_stability_regions.m','script_text_same_order.m','script_text_startup.m', ...
    'script_text_robustness.m','script_text_order_preserving.m', ...
    'script_text_sharp_asymptotics.m','script_text_shrinking_thresholds.m', ...
    'script_text_vector_forcing.m','script_text_coupled_system.m','script_verify_80digits.m'};
run_script_list(paths,names,'full');
check_output_files(paths);
text = sprintf('Full MATLAB reproduction and verification passed on %s.\nMATLAB: %s\n',datestr(now,30),version);
write_text_file(fullfile(paths.checks,'verification_full.txt'),text);
fprintf('%s',text);
end
