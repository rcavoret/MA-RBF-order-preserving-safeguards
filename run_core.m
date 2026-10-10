% Reproduce every numbered table and figure.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
check_invariants();
run_script_list(paths,{'script_table1.m','script_table2.m','script_table3.m','script_table4.m', ...
    'script_stability_regions.m'},'core');
check_output_files(paths);
end
