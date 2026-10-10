function save_results(paths, name, headers, values, settings)
%SAVE_RESULTS Save computed data, full-precision CSV and execution metadata.
% Nothing is loaded from a previous run; generated files are never inputs.
if size(values,2) ~= numel(headers)
    error('MARBF:Dimensions','CSV header/data column count does not agree.');
end
metadata.matlab_version = version;
metadata.computer = computer;
metadata.created = datestr(now,30);
metadata.settings = settings;
file = fullfile(paths.data,[name '.csv']);
fid = fopen(file,'w');
if fid < 0, error('MARBF:OutputFile','Cannot write %s.',file); end
closer = onCleanup(@() fclose(fid));
for j = 1:numel(headers)
    if j > 1, fprintf(fid,','); end
    fprintf(fid,'%s',headers{j});
end
fprintf(fid,'\n');
for i = 1:size(values,1)
    for j = 1:size(values,2)
        if j > 1, fprintf(fid,','); end
        fprintf(fid,'%.17g',values(i,j));
    end
    fprintf(fid,'\n');
end
clear closer
save(fullfile(paths.data,[name '.mat']),'headers','values','metadata');
fprintf('Saved computed data: %s\n',file);
end
