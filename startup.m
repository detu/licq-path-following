% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


% disable warning message
warning off;

if(isunix)
   addpath(genpath('/home/detu/work/casadi355'));
end

