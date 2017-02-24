function root = ROOTDIR
%Retrieve full path of Toolbox installation directory.
%
% SYNOPSIS:
%   root = ROOTDIR
%
% PARAMETERS:
%   none.
%
% RETURNS:
%   root - Full path to code installation directory.

%{
#COPYRIGHT#
%}

nm = mfilename('fullpath');
ix = strfind(nm, filesep);
if ~isempty(ix),
   root = nm(1 : ix(end-1));
else
   root = nm;
end
