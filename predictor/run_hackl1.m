%RUN_HACKL1 Summary of this function goes here
% 
% [OUTPUTARGS] = RUN_HACKL1(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/14 19:51:31 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

%[primal, dual, info] = jpredictor(@(p)hackl1(p) , 0, 1, [0;1], [0.5;0], 1);
%[primal, dual, info] = jpredictorn1(@(p)hackl1(p) , 0, 1, [0;1], [0.5;0], 1, 10, 1);
[primal, dual, info] = jpredictor_tr(@(p)hackl1(p), 0, 1, [0;1], [0.5;0], .1, 0);

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

% figure(3);
% plot(x1,x2,'--rs', C{:});
% xlabel('x1');
% ylabel('x2');
% title('x1 vs. x2');