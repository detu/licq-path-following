%RUN_PROB4 Summary of this function goes here
% 
% caller for problem 4
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/09/16 22:22:59 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

format long;
% addpath('C:\Users\detu\Documents\casadi-matlabR2014b-v2.4.1')
% import casadi.*

xstart = [0.496811686457911; 0.608467572142176];
ystart = [1.2169;0];
%[primal, dual, info] = jpredictor_casadi(@(p)prob4_casadi(p) , [1;-4], [8;1], xstart, [1.2169;0], .1, [], [], 0);

% ystart = [0.735758882342878; 23.729329433526743];
% xstart = [1.000000000000000; 0.367879441171442];

%xstart = [0.496811686457911; 0.608467572142176];
%ystart = [0.735758900642395; 23.729329785913365];
% xstart = [-4.000000000000000; 54.598150033144236];
% ystart = 1.0e+03 * [0.109196300945695; -5.913916202343740];
% 
% [primal, dual, info] = jpredictor_casadi(@(p)prob4_casadi(p), [8;1], [1;-4], xstart, dual, .1, [], [], 0);

nlpRun = 1;
% %[primal, dual, info] = jpredictor_casadi_nlp(@(p)prob4_casadi, [8;1], [1;-4], xstart, ystart, .1, [], [], 0, nlpRun);
[primal, dual, info] = jpredictor_casadi_nlp(@(p)prob4_casadi, [1;-4], [8;1], xstart, ystart, .1, [], [], 0, nlpRun);


% ploting 
numInfo = numel(info);
t  = zeros(1,numInfo);
x1 = zeros(1,numInfo);
x2 = zeros(1,numInfo);
for i=1:numInfo
    t(1,i)  = info(i).t;
    x1(1,i) = info(i).x(1,:);
    x2(1,i) = info(i).x(2,:);
end

xsolution = [x1(end) x2(end)];

C(:,1) = {'LineWidth'; 2};
C(:,2) = {'MarkerEdgeColor'; 'k'};
C(:,3) = {'MarkerFaceColor'; 'g'};
%%

%figure(1)
subplot(2,1,1)
plot(t,x1,'--rs', C{:});
xlabel('t (parameter)'); 
ylabel('x1');
title('x1 vs. t');

%figure(2);
subplot(2,1,2)
plot(t,x2,'--rs', C{:});
xlabel('t (parameter)');
ylabel('x2');
title('x2 vs. t');

