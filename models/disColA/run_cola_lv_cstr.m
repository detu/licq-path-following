clear all;
% load cstrDistInit.mat;
% xmeasure      = Xinit;
[t,x]   = ode15s('cola_lv_cstr',[0 200000],0.5*ones(84,1));
%[t,x]   = ode15s('cola_lv_cstr',[0 200000], xmeasure);
lengthx = size(x); 
Xinit   = x(lengthx(1),:)';

save cstrDistInit11.mat Xinit;