%PLOTCOLLOCATIONDISCOLA Summary of this function goes here
%
% Load collocation results from optDistillation.m and plot the states
% values compare them with steady-state results.
% 
% [OUTPUTARGS] = PLOTCOLLOCATIONDISCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2016/04/12 16:20:08 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

load nlp0.mat;

nk = 50;  % hardcode! should be from optDistillation output
nx = 82;  % this also hardcode. number of states
nu = 2;   % number of control
d  = 3;   % number of collocation points

% optimized initial state
x0_opt       = x_opt(1:nx);  % load optimized decision variables from IPOPT (orthogonal collocation)
x_opt(1:nx)  = [];
x_opt        = reshape(x_opt, nu + nx*d + nx, nk);
u_nlp_opt    = x_opt(1:2,1:nk);
x_opt(1:2,:) = [];     % remove optimized controls

% prepare a matrix for plotting
nState    = nx + nk*(d+1)*nx;    
nPoint    = nState/nx;
plotState = zeros(nx,nPoint);
plotState(:,1) = x0_opt;

% extract states from each collocation point and each time horizon
sInd    = 2; % initial index row
for i=1:nk
    temp    = x_opt(:,i);
    numCol  = size(temp,1);
    numRow  = numCol/nx;
    temp    = reshape(temp,nx,numRow);
    plotState(:,sInd:(numRow+sInd-1)) = temp;
    sInd    = numRow + sInd;
end

% load steady-state solutions
load u_opt_ss.mat;

% plot steady-state and results from IPOPT. make steady-state solution as
% horizontal bold lines.. 
x = linspace(0,nk,nPoint);
disp('State no. 1')
for i=1:nx
    clf;
    plot(x, plotState(i,:),'-.b*'); 
    hold on;
    plot(x, xf(i)*ones(1,nPoint),'-k');
    axis([0 nk 0.0 1.0])
    title(['light component composition at tray. ' int2str(i)])
    xlabel('time [minutes]');
    ylabel('concentration or hold-up');
    %keyboard;
    %pause;
    w = waitforbuttonpress;
%     if w == 0
%         %disp('Button click')
%         fprintf('State no. %d \n', i+1);
%     else
%         disp('Key press')
%     end
end