function run_script_list(paths, names, report_name)
%RUN_SCRIPT_LIST Execute each driver in a fresh function workspace and record timings.
elapsed = zeros(numel(names),1);
for k = 1:numel(names)
    fprintf('\n=== %s ===\n',names{k});
    timer = tic;
    execute_one(fullfile(paths.root,names{k}));
    elapsed(k) = toc(timer);
    fprintf('Completed in %.3f seconds.\n',elapsed(k));
end
text = sprintf('script,elapsed_seconds\n');
for k = 1:numel(names)
    text = [text,sprintf('%s,%.6f\n',names{k},elapsed(k))]; %#ok<AGROW>
end
write_text_file(fullfile(paths.checks,[report_name '_timings.csv']),text);
end

function execute_one(file)
run(file);
end
