function write_text_file(file, text)
%WRITE_TEXT_FILE Checked text output; no implicit working-directory dependency.
fid = fopen(file,'w');
if fid < 0, error('MARBF:OutputFile','Cannot write %s.',file); end
closer = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end
