function check_output_files(paths)
%CHECK_OUTPUT_FILES Confirm all numbered artifacts exist and basic formats are valid.
for k = 1:4
    file = fullfile(paths.tables,sprintf('table%d.tex',k));
    assert(exist(file,'file') == 2,'Missing LaTeX output: %s',file);
    text = fileread(file);
    assert(~isempty(strfind(text,'\begin{table}')) && ~isempty(strfind(text,'\end{table}'))); %#ok<STREMP>
    assert(isempty(strfind(text,'NaN')) && isempty(strfind(text,'Inf'))); %#ok<STREMP>
end
for stem = {'figure1','figure1a','figure1b','figure1c'}
    for extension = {'.pdf','.eps','.png'}
        file = fullfile(paths.figures,[stem{1},extension{1}]);
        info = dir(file);
        assert(~isempty(info) && info.bytes > 500,'Missing/empty figure: %s',file);
    end
    info = imfinfo(fullfile(paths.figures,[stem{1},'.png']));
    assert(info.Width > 2000 && info.Height > 1000,'PNG export is unexpectedly small.');
end
fprintf('All four numbered LaTeX tables and Figure 1 exports are present.\n');
end
