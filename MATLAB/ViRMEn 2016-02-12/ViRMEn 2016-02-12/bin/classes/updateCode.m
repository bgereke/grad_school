function exper = updateCode(exper)

fl = func2str(exper.experimentCode);
if exist([fl '.m'],'file')
    exper.codeText = {};
    fid = fopen([fl '.m']);
    while 1
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        exper.codeText{end+1} = tline;
    end
    fclose(fid);
else
    exper.codeText = {};
end