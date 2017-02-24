function plotStates(data, lb, ub)
%PLOTSTATES Summary of this function goes here
%
% PLOT ALL STATES during MPC iterations.
% ( and also the control, too?)
%
% [OUTPUTARGS] = PLOTSTATES(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/13 15:32:40 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

%nk = 100;  % hardcode! should be from optDistillation output
%nx = 82;  % this also hardcode. number of states
%nu = 2;   % number of control
%d  = 3;   % number of collocation points

global nk nx nu d;

% optimized initial state
x0_opt       = data(1:nx);
data(1:nx)   = [];
data         = reshape(data, nu + nx*d + nx, nk);
u_nlp_opt    = data(1:nu,1:nk);
data(1:nu,:) = [];     % remove optimized controls

lb0          = lb(1:nx);
lb(1:nx)     = [];
lb           = reshape(lb, nu + nx*d + nx, nk);
lbU          = lb(1:nu,1:nk);
lb(1:nu,:)   = [];

ub0          = ub(1:nx);
ub(1:nx)     = [];
ub           = reshape(ub, nu + nx*d + nx, nk);
ubU          = ub(1:nu,1:nk);
ub(1:nu,:)   = [];

% plot optimized control
clf
figure(1);
plot(u_nlp_opt(1,:),'LineWidth',2.5);
hold on
plot(u_nlp_opt(2,:),'LineWidth',2.5);
hold on;
plot(lbU(1,:),'b--.');
hold on;
plot(lbU(2,:),'r--.');
hold on;
plot(ubU(1,:),'b--.');
hold on;
plot(ubU(2,:),'r--.');
legend('u[1]-LT','u[2]-VB');
ylabel('control inputs');
xlabel('time [minutes]');
title('optimized control inputs');

% prepare a matrix for plotting
nState    = nx + nk*(d+1)*nx;
nPoint    = nState/nx;
plotState = zeros(nx,nPoint);
plotState(:,1) = x0_opt;

plotLb         = zeros(nx,nPoint);
plotLb(:,1)    = lb0;

plotUb         = zeros(nx,nPoint);
plotUb(:,1)    = ub0;

% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:nk
    temp    = data(:,i);
    numCol  = size(temp,1);
    numRow  = numCol/nx;
    temp    = reshape(temp,nx,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;

    tempLb = lb(:,i);
    tempLb = reshape(tempLb,nx,numRow);
    plotLb(:,sInd:(numRow+sInd-1)) = tempLb;

    tempUb = ub(:,i);
    tempUb = reshape(tempUb,nx,numRow);
    plotUb(:,sInd:(numRow+sInd-1)) = tempUb;

    sInd    = numRow + sInd;
end

% load steady-state solutions
load u_opt_ss.mat;

% % plot steady-state and results from IPOPT. make steady-state solution as
% % horizontal bold lines..
x = linspace(0,nk,nPoint);
% disp('State no. 1')
% figure(2);
% for i=1:nx
%     clf;
%     plot(x, plotState(i,:),'-.b*');
%     hold on;
%     plot(x, xf(i)*ones(1,nPoint),'-k');
%     hold on;
%     plot(x, plotLb(i,:),'-r');
%     hold on;
%     plot(x, plotUb(i,:),'-r');
%     %axis([0 nk -0.2 1.2])
%     title(['light component at tray no. ' int2str(i)])
%     xlabel('time [minutes]');
%     ylabel('concentration or hold-up');
%     %keyboard;
%     %pause;
%     w = waitforbuttonpress;
%     if w == 0
%         %disp('Button click')
%         fprintf('State no. %d \n', i+1);
%     else
%         disp('Key press')
%     end
% end

% plot light component composition:
% trays no: 1(42), 5(46), 15(56), 21(62), 25(66), 35(76), and 41(82)
figure(2);
clf;
plot(x, plotState(1,:));
hold on;
plot(x, xf(1)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(5,:));
hold on;
plot(x, xf(5)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(15,:));
hold on;
plot(x, xf(15)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(21,:));
hold on;
plot(x, xf(21)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(25,:));
hold on;
plot(x, xf(25)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(35,:));
hold on;
plot(x, xf(35)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(41,:));
hold on;
plot(x, xf(41)*ones(1,nPoint),'-k');

legend('comp. tray #1','comp. tray #5', 'comp. tray #15', 'comp. tray #21', 'comp. tray #25', 'comp. tray #35', 'comp. tray #41');
ylabel('light component composition');
xlabel('time [minutes]');
title('optimized light component compositions');

% plot hold on
figure(3);
clf;
plot(x, plotState(42,:));
hold on;
plot(x, xf(42)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(46,:));
hold on;
plot(x, xf(46)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(56,:));
hold on;
plot(x, xf(56)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(62,:));
hold on;
plot(x, xf(62)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(66,:));
hold on;
plot(x, xf(66)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(76,:));
hold on;
plot(x, xf(76)*ones(1,nPoint),'-k');
hold on;
plot(x, plotState(82,:));
hold on;
plot(x, xf(82)*ones(1,nPoint),'-k');

legend('hold-up #1','hold-up #5', 'hold-up #15', 'hold-up #21', 'hold-up #25', 'hold-up #35', 'hold-up #41');
ylabel('hold-up');
xlabel('time [minutes]');
title('optimized hold-ups');

% Specifically plot xB and xD
figure(4);
clf;
plot(x, plotState(1,:),'-.b*');
hold on;
plot(x, xf(1)*ones(1,nPoint),'-k');
hold on;
% plot(x, plotLb(1,:),'-r');
% hold on;
% plot(x, plotUb(1,:),'-r');
legend('comp. tray #1 [1-xB]');
ylabel('light component composition');
xlabel('time [minutes]');
title('optimized light component compositions');

figure(5);
clf;
plot(x, plotState(41,:),'-.b*');
hold on;
plot(x, xf(41)*ones(1,nPoint),'-k');
% hold on;
% plot(x, plotLb(41,:),'-r');
% hold on;
% plot(x, plotUb(41,:),'-r');
legend('comp. tray #41 [xD]');
ylabel('light component composition');
xlabel('time [minutes]');
title('optimized light component compositions');

end