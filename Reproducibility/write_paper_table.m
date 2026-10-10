function write_paper_table(file, panels, titles, caption, label)
%WRITE_PAPER_TABLE Assemble one complete numbered table from computed panels.
% The resulting .tex needs booktabs and adjustbox, as does the manuscript.
np = numel(panels);
tex = sprintf('%% Generated from computed, unrounded numerical results.\n');
tex = [tex,sprintf('\\begin{table}[!htbp]\n\\centering\n\\scriptsize\n')];
tex = [tex,sprintf('\\setlength{\\tabcolsep}{2.2pt}\n')];
for j = 1:np
    if np == 2
        tex = [tex,sprintf('\\begin{minipage}[t]{0.48\\linewidth}\n\\centering\n')]; %#ok<AGROW>
    end
    if ~isempty(titles{j})
        tex = [tex,sprintf(['\\parbox[t][2\\baselineskip][t]{\\linewidth}{', ...
            '\\centering\\textbf{%s}}\\par\\smallskip\n'],titles{j})]; %#ok<AGROW>
    end
    tex = [tex,sprintf('\\begin{adjustbox}{max width=\\linewidth}\n'), ...
        panels{j},sprintf('\\end{adjustbox}\n')]; %#ok<AGROW>
    if np == 2
        tex = [tex,sprintf('\\end{minipage}%%\n')]; %#ok<AGROW>
        if j == 1, tex = [tex,sprintf('\\hfill%%\n')]; end %#ok<AGROW>
    end
end
tex = [tex,sprintf('\\caption{%s}\n\\label{%s}\n\\end{table}\n',caption,label)];
write_text_file(file,tex);
fprintf('Saved LaTeX: %s\n',file);
end
