% Run the native base-MATLAB verification suite and regenerate its checked outputs.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
check_invariants();
names = {'script_table1.m','script_table2.m','script_table3.m','script_table4.m', ...
    'script_stability_regions.m','script_text_same_order.m','script_text_startup.m', ...
    'script_text_robustness.m','script_text_shrinking_thresholds.m', ...
    'script_text_vector_forcing.m','script_text_coupled_system.m'};
run_script_list(paths,names,'verification_core');
check_output_files(paths);
text = sprintf(['Core MATLAB verification passed on %s.\nMATLAB: %s\n', ...
    'Includes all numbered outputs and all double-precision textual checks.\n', ...
    '80-digit checks require run_all.m or the separately identified precision scripts.\n'], ...
    datestr(now,30),version);
write_text_file(fullfile(paths.checks,'verification_core.txt'),text);
fprintf('%s',text);
end
