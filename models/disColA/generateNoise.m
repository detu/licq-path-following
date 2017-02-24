%GENERATENOISE Summary of this function goes here
% 
% [OUTPUTARGS] = GENERATENOISE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/10/23 21:09:58 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

%mpciterations = 150;
mpciterations = 500;
noiseLevel    = 0.01; % 1 percent
%noiseLevel    = 0.03; % 3 percent
%noiseLevel    = 0.05; % 5 percent
load CstrDistXinit.mat;
xf      = Xinit(1:84);
xholdup = xf(43:end);

noise = zeros(42,mpciterations);
for i=1:mpciterations
    noise(:,i) = noiseLevel * xholdup .*randn(42,1);
end

save noise1pct.mat noise;
%save noise3pct.mat noise;
%save noise5pct.mat noise;
