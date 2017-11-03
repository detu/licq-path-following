function [u_nlp_opt,plotState] = plotStatesGL(data, lb, ub, N)
%PLOTSTATESN Summary of this function goes here
% 
% [OUTPUTARGS] = PLOTSTATESN(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/27 22:42:29 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

%global nk nx nu d ns;
nk = 24;
nx = 8;
nz = 30;
nu = 2;
d  = 3;
ns = 0;

% optimized initial state
x0_opt       = data(1:nx+nz);  
data(1:nx+nz)= [];
data         = reshape(data, (nu + (nx+nz+ns)*d + (nx+ns)), nk);
u_nlp_opt    = data(1:nu,1:nk);
data(1:nu,:) = [];     % remove optimized controls



% prepare a matrix for plotting
nState    = (nx+nz+ns) + nk*(d+1)*(nx+nz+ns);
nPoint    = nState/(nx+nz+ns);
plotState = zeros(nx+nz+ns,nPoint);
plotState(1:nx+nz,1) = x0_opt;

% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:nk
    temp    = [data(:,i);data(85:114,i)];
    numCol  = size(temp,1);
    numRow  = numCol/(nx+nz+ns);
    temp    = reshape(temp,nx+nz+ns,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;
    
    sInd    = numRow + sInd;
end

end
