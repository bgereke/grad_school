function virmenMake

% ----------------------------------------------
% ViRMEn Makefile
% ----------------------------------------------
% Compile all .c files in the Virmen directory,
% only including GLFW in files that need it.
% ----------------------------------------------

mfile = mfilename('fullpath');
path = fileparts(mfile);
backupPath = pwd;
cd(path);
try
    compilationFunction;
catch err%#ok<CTCH>
    rethrow(err)
end
cd(backupPath);


function compilationFunction

% store some useful paths
binPath = pwd;
cd engine;
enginePath = pwd;
cd ..;

% get a cell array of all Virmen subdirectories
cd ..;
dirString = genpath(pwd);
dirArray = regexp(dirString, pathsep, 'split');
dirArray = dirArray(2:end-1);

% crawl all sub-directories for c files
for i = 1:length(dirArray)
   currentDir = dirArray{i};
   cd(currentDir);
   cFilesArray = struct2cell([dir('*.c'); dir('*.cpp')]);
   cFiles = cFilesArray(1,:);
   
   % compile all c files in this directory
   for j = 1:length(cFiles)
        f = cFiles{j};
        
        wasCopied = false;   % flag so we remember to delete copies
        
        % check if this file requires GLFW
        fid = fopen(f, 'r'); 
        isGLFW=false; 
        for ndx = 1:5 
            str = fgetl(fid); 
            if ~isempty(strfind(str,'GLFW')) 
                isGLFW=true;
            end
        end
        fclose(fid);
        
        % compile the file
        if ~isGLFW
            disp(['Compiling ', f, ' without GLFW...']);
            mex(f);
        else
            disp(['Compiling ', f, ' with GLFW...']);
            % copy file to the engine directory to compile with GLFW,
            % unless we are already in the engine directory
            if (~strcmp(pwd, enginePath))
                wasCopied = true;
                copyfile(f, enginePath);
                cd(enginePath);
            end
            % try to compile a couple times in case of spurious error
            for attempts = 1:3
                try
                    % compile depending on architecture
                    if strcmp(computer, 'PCWIN')
                        mex('-LGLFW_32','-lglfw3','-lopengl32','-outdir',currentDir,f);
                    elseif strcmp(computer, 'PCWIN64')
                        mex('-LGLFW','-lglfw3','-lopengl32','-outdir', currentDir, f);
                    elseif strcmp(computer, 'MACI64') % '-Duint16_t=uint16_T',
                        mex('-L./GLFW','-v','-lglfw.3','LDFLAGS="\$LDFLAGS','-framework','Cocoa','-framework','OpenGL','-framework','IOKit','-framework','CoreVideo"','-outdir',currentDir,f);
                    else
                        if wasCopied == true
                            delete(f);
                        end
                        error(strcat('Unsupported compilation architecture.', ...
                              'See ViRMEn manual for compilation instructions.'));
                    end
                    break; % if it worked, move on
                catch err
                    if (attempts > 2)
                        if (wasCopied == true)
                            delete(f);
                        end
                        rethrow(err);
                    else
                        continue;
                    end
                end
            end
            
            % if we copied into engine directory, delete the copy and go back
            % to the current directory to continue the crawl
            if (wasCopied == true)
                delete(f);
                cd(currentDir);
            end
            
        end
        
        disp(['Successfully compiled ', f]);
        disp('----------------------------');
   end

end

% navigate back to the bin directory where we started
cd(binPath);