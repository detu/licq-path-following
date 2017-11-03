%RUN_KYPARISIS Summary of this function goes here
% 
% [OUTPUTARGS] = RUN_KYPARISIS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/05 14:36:17 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

% QUESTION: why work fine with small initial delta_t ??? OR close to 1 !

%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [-1;-1], [1;0], [0;0], [0;0;1]); % WRONG initial dual variables!!! 
%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [-1;-1], [-0.5;0], [0;0], [1;1;0], 0.8);
%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [-1;-1], [1;0], [0;0], [1;1;0], 0.3);
%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [.5;.5], [0.1;0.2], [.5;.5], [0;0;1], 0.5);
%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [-1;-1], [0;0], [0;0], [1;1;0], 0.001);

%[primal, dual, info] = jpredictor(@(p)kyparisis(p) , [.5;.5], [0.1234;0.5678], [.5;.5], [0;0;1], 1);
[primal, dual, info] = jpredictorn1(@(p)kyparisis(p) , [.5;.5], [0.1234;0.5678], [.5;.5], [0;0;1], 1, 5, 0);

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

C(:,1) = {'LineWidth'; 2};
C(:,2) = {'MarkerEdgeColor'; 'k'};
C(:,3) = {'MarkerFaceColor'; 'g'};

% figure(1)
% plot(t,x1,'--rs', C{:});
% xlabel('t (parameter)'); 
% ylabel('x1');
% title('x1 vs. t');
% 
% figure(2);
% plot(t,x2,'--rs', C{:});
% xlabel('t (parameter)');
% ylabel('x2');
% title('x2 vs. t');
% 
% figure(3);
% plot(x1,x2,'--rs', C{:});
% xlabel('x1');
% ylabel('x2');
% title('x1 vs. x2');

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

