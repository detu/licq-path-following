%RUN_HACKL2 Summary of this function goes here
% 
% [OUTPUTARGS] = RUN_HACKL2(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2015/10/14 22:44:19 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

lb = [0;0];
ub = [inf;inf];
%[primal, dual, info] = jpredictor(@(p)hackl2(p) , 0, 1, [0;0], [2;2;0;0], 1);
%[primal, dual, info] = jpredictorn1(@(p)hackl2(p) , 0, 1, [0;0], [2;2;0;0], 1, 10, 0);
%[primal, dual, info] = jpredictor_tr(@(p)hackl2(p), 0, 1, [0;0], [2;2;0;0], 1, 0);
[primal, dual, info] = jpredictor_tr(@(p)hackl2b(p), 0, 1, [0;0], [0;0], .2, lb, ub, 0);

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
