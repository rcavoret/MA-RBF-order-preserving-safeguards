function tex = table_latex(met_name, N, met_err, met_ord, options)
%TABLE_LATEX Format COMPUTED errors/rates; compatible with the original 4 inputs.
% With no output argument, print the table as in the original repository.
% With an output argument, return its LaTeX source. OPTIONS may contain:
% fragment, maximum_error, counts, shared_count, error_digits, filename.
% Counts are measured by the integrator, never inferred from exact histories.
if nargin < 5, options = struct(); end
if ~isfield(options,'fragment'), options.fragment = false; end
if ~isfield(options,'maximum_error'), options.maximum_error = false; end
if ~isfield(options,'counts'), options.counts = []; end
if ~isfield(options,'shared_count'), options.shared_count = []; end
if ~isfield(options,'error_digits'), options.error_digits = 2; end
if ~isfield(options,'filename'), options.filename = ''; end
if iscell(met_err)
    for j = 1:numel(met_err), met_err{j} = met_err{j}(:); end
    met_err = [met_err{:}];
end
if iscell(met_ord)
    for j = 1:numel(met_ord), met_ord{j} = met_ord{j}(:); end
    met_ord = [met_ord{:}];
end
N = N(:); nm = numel(met_name);
if ~isequal(size(met_err),[numel(N),nm]) || ~isequal(size(met_ord),size(met_err))
    error('MARBF:TableDimensions','Table dimensions do not agree.');
end
nc = 2+~isempty(options.counts);
if ~isempty(options.counts) && ~isequal(size(options.counts),size(met_err))
    error('MARBF:TableCounts','Expected one activation count per error.');
end
if options.maximum_error
    ename = '$\mathcal E_N$'; pname = '$p_N^{\max}$';
else
    ename = '$E_N$'; pname = '$p_N$';
end
extra = '';
if ~isempty(options.shared_count), extra = ' r'; end
spec = ['@{}r',repmat(repmat(' r',1,nc),1,nm),extra,'@{}'];
tex = '';
if ~options.fragment
    tex = sprintf('\\begin{table}[!htbp]\n\\centering\n\\small\n\\begin{adjustbox}{max width=\\linewidth}\n');
end
tex = [tex,sprintf('\\begin{tabular}{%s}\n\\toprule\n',spec)];
for m = 1:nm
    tex = [tex,sprintf(' & \\multicolumn{%d}{c}{%s}',nc,met_name{m})]; %#ok<AGROW>
end
if ~isempty(options.shared_count), tex = [tex,' & ']; end
tex = [tex,sprintf(' \\\\\n')];
for m = 1:nm
    lo = 2+(m-1)*nc;
    tex = [tex,sprintf('\\cmidrule(lr){%d-%d}',lo,lo+nc-1)]; %#ok<AGROW>
end
tex = [tex,sprintf('\n$N$')];
for m = 1:nm
    tex = [tex,' & ',ename,' & ',pname]; %#ok<AGROW>
    if nc == 3, tex = [tex,' & $A_N$']; end %#ok<AGROW>
end
if ~isempty(options.shared_count), tex = [tex,' & $A_N$']; end
tex = [tex,sprintf(' \\\\\n\\midrule\n')];
for i = 1:numel(N)
    tex = [tex,sprintf('%d',N(i))]; %#ok<AGROW>
    for m = 1:nm
        e = met_err(i,m); p = met_ord(i,m);
        if isnan(e), es = '--'; else, es = sprintf('%.*E',options.error_digits,e); end
        if isnan(p), ps = '--'; else, ps = sprintf('%.4f',p); end
        tex = [tex,sprintf(' & %s & %s',es,ps)]; %#ok<AGROW>
        if nc == 3
            tex = [tex,sprintf(' & %d',options.counts(i,m))]; %#ok<AGROW>
        end
    end
    if ~isempty(options.shared_count)
        tex = [tex,sprintf(' & %d',options.shared_count(i))]; %#ok<AGROW>
    end
    tex = [tex,sprintf(' \\\\\n')]; %#ok<AGROW>
end
tex = [tex,sprintf('\\bottomrule\n\\end{tabular}\n')];
if ~options.fragment, tex = [tex,sprintf('\\end{adjustbox}\n\\end{table}\n')]; end
if ~isempty(options.filename), write_text_file(options.filename,tex); end
if nargout == 0, fprintf('%s',tex); end
end
