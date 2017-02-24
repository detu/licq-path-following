function x0 = computeInitGuess(xA, xB, numHor)
%COMPUTEINITGUESS Summary of this function goes here
% 
% [OUTPUTARGS] = COMPUTEINITGUESS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/08/25 15:25:13 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

nc = 3;                  % number of collocation point
nx = size(xA,1);
%x0 = zeros(nx,numHor);
x0 = zeros(nx,numHor*(nc+1));

for i=1:numHor
    t       = i/numHor;
    %x0(:,i) = (1-t)*xA + t*xB;
    iS  = 4*i - 3;
    iE  = 4*i;
    temp = (1-t)*xA + t*xB;
    x0(:,iS:iE) = repmat(temp,1,4);
end

end
