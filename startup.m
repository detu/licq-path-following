% Startup code should be run at the beginning. 
% This code is to include all the directories (set path). 

folder = fileparts(which(mfilename));
addpath(genpath(folder));

% windows path
if (ispc)
    addpath('C:\Data\work\casadi-matlabR2014b-v3.1.0-rc1');
end
% linux path
if(isunix)
  addpath(genpath('/home/detu/work/casadi355'));
end
clear p;
cd ..
