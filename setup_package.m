function paths = setup_package()
%SETUP_PACKAGE Add only the package's source folders and create output folders.
% Paths are absolute and do not depend on the current working directory.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'Classical methods for IVP'));
addpath(fullfile(root, 'RBF methods for IVP'));
addpath(fullfile(root, 'Stability regions'));
addpath(fullfile(root, 'Reproducibility'));
paths.root = root;
paths.tables = fullfile(root, 'output', 'tables');
paths.figures = fullfile(root, 'output', 'figures');
paths.data = fullfile(root, 'output', 'data');
paths.checks = fullfile(root, 'output', 'checks');
names = {'tables', 'figures', 'data', 'checks'};
for j = 1:numel(names)
    if ~exist(paths.(names{j}), 'dir')
        [ok, message] = mkdir(paths.(names{j}));
        if ~ok, error('MARBF:OutputDirectory', '%s', message); end
    end
end
end
