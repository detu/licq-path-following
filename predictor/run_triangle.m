%RUN_TRIANGLE Summary of this function goes here
% 
% [OUTPUTARGS] = RUN_TRIANGLE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2015/10/19 15:53:31 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2015

xstart = [1.5; 0.0; 1.5; 1.0; 0.0 ];
ystart = 1e-7*[-0.198682149251302; 0; -0.074505805969238]; 

%[primal, dual, info] = jpredictor(@(p)triangle(p) , 0, 1, xstart, ystart, 0.75);
[primal, dual, info] = jpredictor_tr(@(p)triangle(p) , 0, 1, xstart, ystart, .1, [], [], 0);

% ploting 
numInfo = numel(info);
t  = zeros(1,numInfo);
x1 = zeros(1,numInfo);
x2 = zeros(1,numInfo);
x3 = zeros(1,numInfo);
x4 = zeros(1,numInfo);
x5 = zeros(1,numInfo);
for i=1:numInfo
    t(1,i)  = info(i).t;
    x1(1,i) = info(i).x(1,:);
    x2(1,i) = info(i).x(2,:);
    x3(1,i) = info(i).x(3,:);
    x4(1,i) = info(i).x(4,:);
    x5(1,i) = info(i).x(5,:);
end

C(:,1) = {'LineWidth'; 2};
C(:,2) = {'MarkerEdgeColor'; 'k'};
C(:,3) = {'MarkerFaceColor'; 'g'};

%figure(1)
subplot(5,1,1)
plot(t,x1,'--rs', C{:});
xlabel('t (parameter)'); 
ylabel('x1');
title('x1 vs. t');

%figure(2);
subplot(5,1,2)
plot(t,x2,'--rs', C{:});
xlabel('t (parameter)');
ylabel('x2');
title('x2 vs. t');

subplot(5,1,3)
plot(t,x3,'--rs', C{:});
xlabel('t (parameter)');
ylabel('x3');
title('x3 vs. t');

subplot(5,1,4)
plot(t,x4,'--rs', C{:});
xlabel('t (parameter)');
ylabel('x4');
title('x4 vs. t');

subplot(5,1,5)
plot(t,x5,'--rs', C{:});
xlabel('t (parameter)');
ylabel('x5');
title('x5 vs. t');
