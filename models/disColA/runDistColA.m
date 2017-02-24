%RUNDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = RUNDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: detu $	$Date: 2016/01/24 22:35:30 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2016

% TODO:
% - compare with NLP solution !
% - plot results to compare PF and IPOPT 

%===============================================================================
% Initiate guess both for primal and dual variables
%===============================================================================
load nlp0.mat;
xstart = x_opt;
ystart = y_opt;


p_init  = 0.50;   %concentration in reactor
p_final = 0.55;   % WARNING: other values may give QP error: -6 (non-convex problem)
%p_final = 0.60;
%p_final = 0.501;
%p_final = 0.504;
%delta_t = 0.2;
delta_t = 0.1;
%delta_t = 0.05;

nlpRun = 1;
%nlpRun = 0;

lb_init = lb;
ub_init = ub;

if (nlpRun == 1)
    % return primal and dual variables from NLP solver and PF algorithm
    [primal, dual, info, nlp] = jpredictor_casadi_nlp(@(p)distColA_casadi(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, nlpRun);
else
    % just return primal and dual variables from PF algorithm
    [primal, dual, info] = jpredictor_casadi_nlp(@(p)distColA_casadi(p), p_init, p_final, xstart, ystart, delta_t, lb_init, ub_init, 0, nlpRun);
end


%===============================================================================
% PLOTTING
%===============================================================================
nx = 82;
nu = 2;
nk = 10; 
tf = 10.0;
[~,C,D,d] = collocationSetup();
tgrid   = linspace(0,tf,nk+1);
tgrid_u = linspace(0,tf,nk);

u1_opt_pf  = primal(83:(nx+nu)+nx*d:end);
u2_opt_pf  = primal(84:(nx+nu)+nx*d:end);

x1_opt_pf  = primal(1:(nx+nu)+nx*d:end);
x2_opt_pf  = primal(2:(nx+nu)+nx*d:end);
x21_opt_pf = primal(21:(nx+nu)+nx*d:end);
x40_opt_pf = primal(40:(nx+nu)+nx*d:end);
x41_opt_pf = primal(41:(nx+nu)+nx*d:end);
if (nlpRun == 1)
    x1_opt_nlp  = nlp.primal(1:(nx+nu)+nx*d:end);
    x2_opt_nlp  = nlp.primal(2:(nx+nu)+nx*d:end);
    x21_opt_nlp = nlp.primal(21:(nx+nu)+nx*d:end);
    x40_opt_nlp = nlp.primal(40:(nx+nu)+nx*d:end);
    x41_opt_nlp = nlp.primal(41:(nx+nu)+nx*d:end);
    
    u1_opt_nlp  = nlp.primal(83:(nx+nu)+nx*d:end);
    u2_opt_nlp  = nlp.primal(84:(nx+nu)+nx*d:end);
end


% x1_opt  = primal(1:(nx+nu)+nx*d:end);
% x2_opt  = primal(2:(nx+nu)+nx*d:end);
% x3_opt  = primal(3:(nx+nu)+nx*d:end);
% x4_opt  = primal(4:(nx+nu)+nx*d:end);
% x5_opt  = primal(5:(nx+nu)+nx*d:end);
% x6_opt  = primal(6:(nx+nu)+nx*d:end);
% x7_opt  = primal(7:(nx+nu)+nx*d:end);
% x8_opt  = primal(8:(nx+nu)+nx*d:end);
% x9_opt  = primal(9:(nx+nu)+nx*d:end);
% x10_opt = primal(10:(nx+nu)+nx*d:end);
% x11_opt = primal(11:(nx+nu)+nx*d:end);
% x12_opt = primal(12:(nx+nu)+nx*d:end);
% x13_opt = primal(13:(nx+nu)+nx*d:end);
% x14_opt = primal(14:(nx+nu)+nx*d:end);
% x15_opt = primal(15:(nx+nu)+nx*d:end);
% x16_opt = primal(16:(nx+nu)+nx*d:end);
% x17_opt = primal(17:(nx+nu)+nx*d:end);
% x18_opt = primal(18:(nx+nu)+nx*d:end);
% x19_opt = primal(19:(nx+nu)+nx*d:end);
% x20_opt = primal(20:(nx+nu)+nx*d:end);
% x21_opt = primal(21:(nx+nu)+nx*d:end);
% x22_opt = primal(22:(nx+nu)+nx*d:end);
% x23_opt = primal(23:(nx+nu)+nx*d:end);
% x24_opt = primal(24:(nx+nu)+nx*d:end);
% x25_opt = primal(25:(nx+nu)+nx*d:end);
% x26_opt = primal(26:(nx+nu)+nx*d:end);
% x27_opt = primal(27:(nx+nu)+nx*d:end);
% x28_opt = primal(28:(nx+nu)+nx*d:end);
% x29_opt = primal(29:(nx+nu)+nx*d:end);
% x30_opt = primal(30:(nx+nu)+nx*d:end);
% x31_opt = primal(31:(nx+nu)+nx*d:end);
% x32_opt = primal(32:(nx+nu)+nx*d:end);
% x33_opt = primal(33:(nx+nu)+nx*d:end);
% x34_opt = primal(34:(nx+nu)+nx*d:end);
% x35_opt = primal(35:(nx+nu)+nx*d:end);
% x36_opt = primal(36:(nx+nu)+nx*d:end);
% x37_opt = primal(37:(nx+nu)+nx*d:end);
% x38_opt = primal(38:(nx+nu)+nx*d:end);
% x39_opt = primal(39:(nx+nu)+nx*d:end);
% x40_opt = primal(40:(nx+nu)+nx*d:end);
% x41_opt = primal(41:(nx+nu)+nx*d:end);
% x42_opt = primal(42:(nx+nu)+nx*d:end);
% x43_opt = primal(43:(nx+nu)+nx*d:end);
% x44_opt = primal(44:(nx+nu)+nx*d:end);
% x45_opt = primal(45:(nx+nu)+nx*d:end);
% x46_opt = primal(46:(nx+nu)+nx*d:end);
% x47_opt = primal(47:(nx+nu)+nx*d:end);
% x48_opt = primal(48:(nx+nu)+nx*d:end);
% x49_opt = primal(49:(nx+nu)+nx*d:end);
% x50_opt = primal(50:(nx+nu)+nx*d:end);
% x51_opt = primal(51:(nx+nu)+nx*d:end);
% x52_opt = primal(52:(nx+nu)+nx*d:end);
% x53_opt = primal(53:(nx+nu)+nx*d:end);
% x54_opt = primal(54:(nx+nu)+nx*d:end);
% x55_opt = primal(55:(nx+nu)+nx*d:end);
% x56_opt = primal(56:(nx+nu)+nx*d:end);
% x57_opt = primal(57:(nx+nu)+nx*d:end);
% x58_opt = primal(58:(nx+nu)+nx*d:end);
% x59_opt = primal(59:(nx+nu)+nx*d:end);
% x60_opt = primal(60:(nx+nu)+nx*d:end);
% x61_opt = primal(61:(nx+nu)+nx*d:end);
% x62_opt = primal(62:(nx+nu)+nx*d:end);
% x63_opt = primal(63:(nx+nu)+nx*d:end);
% x64_opt = primal(64:(nx+nu)+nx*d:end);
% x65_opt = primal(65:(nx+nu)+nx*d:end);
% x66_opt = primal(66:(nx+nu)+nx*d:end);
% x67_opt = primal(67:(nx+nu)+nx*d:end);
% x68_opt = primal(68:(nx+nu)+nx*d:end);
% x69_opt = primal(69:(nx+nu)+nx*d:end);
% x70_opt = primal(70:(nx+nu)+nx*d:end);
% x71_opt = primal(71:(nx+nu)+nx*d:end);
% x72_opt = primal(72:(nx+nu)+nx*d:end);
% x73_opt = primal(73:(nx+nu)+nx*d:end);
% x74_opt = primal(74:(nx+nu)+nx*d:end);
% x75_opt = primal(75:(nx+nu)+nx*d:end);
% x76_opt = primal(76:(nx+nu)+nx*d:end);
% x77_opt = primal(77:(nx+nu)+nx*d:end);
% x78_opt = primal(78:(nx+nu)+nx*d:end);
% x79_opt = primal(79:(nx+nu)+nx*d:end);
% x80_opt = primal(80:(nx+nu)+nx*d:end);
% x81_opt = primal(81:(nx+nu)+nx*d:end);
% x82_opt = primal(82:(nx+nu)+nx*d:end);
% 
% x1_opt_nlp  = nlp.primal(1:(nx+nu)+nx*d:end);
% x2_opt_nlp  = nlp.primal(2:(nx+nu)+nx*d:end);
% x3_opt_nlp  = nlp.primal(3:(nx+nu)+nx*d:end);
% x4_opt_nlp  = nlp.primal(4:(nx+nu)+nx*d:end);
% x5_opt_nlp  = nlp.primal(5:(nx+nu)+nx*d:end);
% x6_opt_nlp  = nlp.primal(6:(nx+nu)+nx*d:end);
% x7_opt_nlp  = nlp.primal(7:(nx+nu)+nx*d:end);
% x8_opt_nlp  = nlp.primal(8:(nx+nu)+nx*d:end);
% x9_opt_nlp  = nlp.primal(9:(nx+nu)+nx*d:end);
% x10_opt_nlp = nlp.primal(10:(nx+nu)+nx*d:end);
% x11_opt_nlp = nlp.primal(11:(nx+nu)+nx*d:end);
% x12_opt_nlp = nlp.primal(12:(nx+nu)+nx*d:end);
% x13_opt_nlp = nlp.primal(13:(nx+nu)+nx*d:end);
% x14_opt_nlp = nlp.primal(14:(nx+nu)+nx*d:end);
% x15_opt_nlp = nlp.primal(15:(nx+nu)+nx*d:end);
% x16_opt_nlp = nlp.primal(16:(nx+nu)+nx*d:end);
% x17_opt_nlp = nlp.primal(17:(nx+nu)+nx*d:end);
% x18_opt_nlp = nlp.primal(18:(nx+nu)+nx*d:end);
% x19_opt_nlp = nlp.primal(19:(nx+nu)+nx*d:end);
% x20_opt_nlp = nlp.primal(20:(nx+nu)+nx*d:end);
% x21_opt_nlp = nlp.primal(21:(nx+nu)+nx*d:end);
% x22_opt_nlp = nlp.primal(22:(nx+nu)+nx*d:end);
% x23_opt_nlp = nlp.primal(23:(nx+nu)+nx*d:end);
% x24_opt_nlp = nlp.primal(24:(nx+nu)+nx*d:end);
% x25_opt_nlp = nlp.primal(25:(nx+nu)+nx*d:end);
% x26_opt_nlp = nlp.primal(26:(nx+nu)+nx*d:end);
% x27_opt_nlp = nlp.primal(27:(nx+nu)+nx*d:end);
% x28_opt_nlp = nlp.primal(28:(nx+nu)+nx*d:end);
% x29_opt_nlp = nlp.primal(29:(nx+nu)+nx*d:end);
% x30_opt_nlp = nlp.primal(30:(nx+nu)+nx*d:end);
% x31_opt_nlp = nlp.primal(31:(nx+nu)+nx*d:end);
% x32_opt_nlp = nlp.primal(32:(nx+nu)+nx*d:end);
% x33_opt_nlp = nlp.primal(33:(nx+nu)+nx*d:end);
% x34_opt_nlp = nlp.primal(34:(nx+nu)+nx*d:end);
% x35_opt_nlp = nlp.primal(35:(nx+nu)+nx*d:end);
% x36_opt_nlp = nlp.primal(36:(nx+nu)+nx*d:end);
% x37_opt_nlp = nlp.primal(37:(nx+nu)+nx*d:end);
% x38_opt_nlp = nlp.primal(38:(nx+nu)+nx*d:end);
% x39_opt_nlp = nlp.primal(39:(nx+nu)+nx*d:end);
% x40_opt_nlp = nlp.primal(40:(nx+nu)+nx*d:end);
% x41_opt_nlp = nlp.primal(41:(nx+nu)+nx*d:end);
% x42_opt_nlp = nlp.primal(42:(nx+nu)+nx*d:end);
% x43_opt_nlp = nlp.primal(43:(nx+nu)+nx*d:end);
% x44_opt_nlp = nlp.primal(44:(nx+nu)+nx*d:end);
% x45_opt_nlp = nlp.primal(45:(nx+nu)+nx*d:end);
% x46_opt_nlp = nlp.primal(46:(nx+nu)+nx*d:end);
% x47_opt_nlp = nlp.primal(47:(nx+nu)+nx*d:end);
% x48_opt_nlp = nlp.primal(48:(nx+nu)+nx*d:end);
% x49_opt_nlp = nlp.primal(49:(nx+nu)+nx*d:end);
% x50_opt_nlp = nlp.primal(50:(nx+nu)+nx*d:end);
% x51_opt_nlp = nlp.primal(51:(nx+nu)+nx*d:end);
% x52_opt_nlp = nlp.primal(52:(nx+nu)+nx*d:end);
% x53_opt_nlp = nlp.primal(53:(nx+nu)+nx*d:end);
% x54_opt_nlp = nlp.primal(54:(nx+nu)+nx*d:end);
% x55_opt_nlp = nlp.primal(55:(nx+nu)+nx*d:end);
% x56_opt_nlp = nlp.primal(56:(nx+nu)+nx*d:end);
% x57_opt_nlp = nlp.primal(57:(nx+nu)+nx*d:end);
% x58_opt_nlp = nlp.primal(58:(nx+nu)+nx*d:end);
% x59_opt_nlp = nlp.primal(59:(nx+nu)+nx*d:end);
% x60_opt_nlp = nlp.primal(60:(nx+nu)+nx*d:end);
% x61_opt_nlp = nlp.primal(61:(nx+nu)+nx*d:end);
% x62_opt_nlp = nlp.primal(62:(nx+nu)+nx*d:end);
% x63_opt_nlp = nlp.primal(63:(nx+nu)+nx*d:end);
% x64_opt_nlp = nlp.primal(64:(nx+nu)+nx*d:end);
% x65_opt_nlp = nlp.primal(65:(nx+nu)+nx*d:end);
% x66_opt_nlp = nlp.primal(66:(nx+nu)+nx*d:end);
% x67_opt_nlp = nlp.primal(67:(nx+nu)+nx*d:end);
% x68_opt_nlp = nlp.primal(68:(nx+nu)+nx*d:end);
% x69_opt_nlp = nlp.primal(69:(nx+nu)+nx*d:end);
% x70_opt_nlp = nlp.primal(70:(nx+nu)+nx*d:end);
% x71_opt_nlp = nlp.primal(71:(nx+nu)+nx*d:end);
% x72_opt_nlp = nlp.primal(72:(nx+nu)+nx*d:end);
% x73_opt_nlp = nlp.primal(73:(nx+nu)+nx*d:end);
% x74_opt_nlp = nlp.primal(74:(nx+nu)+nx*d:end);
% x75_opt_nlp = nlp.primal(75:(nx+nu)+nx*d:end);
% x76_opt_nlp = nlp.primal(76:(nx+nu)+nx*d:end);
% x77_opt_nlp = nlp.primal(77:(nx+nu)+nx*d:end);
% x78_opt_nlp = nlp.primal(78:(nx+nu)+nx*d:end);
% x79_opt_nlp = nlp.primal(79:(nx+nu)+nx*d:end);
% x80_opt_nlp = nlp.primal(80:(nx+nu)+nx*d:end);
% x81_opt_nlp = nlp.primal(81:(nx+nu)+nx*d:end);
% x82_opt_nlp = nlp.primal(82:(nx+nu)+nx*d:end);
% 
% plot(tgrid,x1_opt, tgrid,x1_opt_nlp); hold on;
% plot(tgrid,x2_opt, tgrid,x2_opt_nlp); hold on;
% plot(tgrid,x3_opt, tgrid,x3_opt_nlp); hold on;
% plot(tgrid,x4_opt, tgrid,x4_opt_nlp); hold on;
% plot(tgrid,x5_opt, tgrid,x5_opt_nlp); hold on;
% plot(tgrid,x6_opt, tgrid,x6_opt_nlp); hold on;
% plot(tgrid,x7_opt, tgrid,x7_opt_nlp); hold on;
% plot(tgrid,x8_opt, tgrid,x8_opt_nlp); hold on;
% plot(tgrid,x9_opt, tgrid,x9_opt_nlp); hold on;
% plot(tgrid,x10_opt, tgrid,x10_opt_nlp); hold on;
% plot(tgrid,x11_opt, tgrid,x11_opt_nlp); hold on;
% plot(tgrid,x12_opt, tgrid,x12_opt_nlp); hold on;
% plot(tgrid,x13_opt, tgrid,x13_opt_nlp); hold on;
% plot(tgrid,x14_opt, tgrid,x14_opt_nlp); hold on;
% plot(tgrid,x15_opt, tgrid,x15_opt_nlp); hold on;
% plot(tgrid,x16_opt, tgrid,x16_opt_nlp); hold on;
% plot(tgrid,x17_opt, tgrid,x17_opt_nlp); hold on;
% plot(tgrid,x18_opt, tgrid,x18_opt_nlp); hold on;
% plot(tgrid,x19_opt, tgrid,x19_opt_nlp'); hold on;
% plot(tgrid,x20_opt, tgrid,x20_opt_nlp'); hold on;
% plot(tgrid,x21_opt, tgrid,x21_opt_nlp'); 
% legend('x[1]','x[2]','x[3]','x[4]', 'x[5]', 'x[6]','x[7]','x[8]','x[9]', 'x[10]', 'x[11]','x[12]','x[13]','x[14]', 'x[15]', 'x[16]','x[17]','x[18]','x[19]', 'x[20]', 'x[21]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[1] - x[21]');
% 
% figure(4)
% plot(tgrid,x22_opt, tgrid,x22_opt_nlp','LineWidth',2.5); hold on;
% plot(tgrid,x23_opt, tgrid,x23_opt_nlp','LineWidth',2.5); hold on;
% plot(tgrid,x24_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x25_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x26_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x27_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x28_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x29_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x30_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x31_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x32_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x33_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x34_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x35_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x36_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x37_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x38_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x39_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x40_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x41_opt, 'LineWidth',2.5); 
% legend('x[22]','x[23]','x[24]', 'x[25]', 'x[26]','x[27]','x[28]','x[29]', 'x[30]', 'x[31]','x[32]','x[33]','x[34]', 'x[35]', 'x[36]','x[37]','x[38]','x[39]', 'x[40]', 'x[41]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[22] - x[41]');
% 
% figure(5)
% plot(tgrid,x42_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x43_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x44_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x45_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x46_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x47_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x48_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x49_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x50_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x51_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x52_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x53_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x54_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x55_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x56_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x57_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x58_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x59_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x60_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x61_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x62_opt, 'LineWidth',2.5); 
% legend('x[42]','x[43]','x[44]', 'x[45]', 'x[46]','x[47]','x[48]','x[49]', 'x[50]', 'x[51]','x[52]','x[53]','x[54]', 'x[55]', 'x[56]','x[57]','x[58]','x[59]', 'x[60]', 'x[61]', 'x[62]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[42] - x[62]');
% 
% figure(6)
% plot(tgrid,x63_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x64_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x65_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x66_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x67_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x68_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x69_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x70_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x71_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x72_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x73_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x74_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x75_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x76_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x77_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x78_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x79_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x80_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x81_opt, 'LineWidth',2.5); hold on;
% plot(tgrid,x82_opt, 'LineWidth',2.5);
% legend('x[63]','x[64]', 'x[65]', 'x[66]','x[67]','x[68]','x[69]', 'x[70]', 'x[71]','x[72]','x[73]','x[74]', 'x[75]', 'x[76]','x[77]','x[78]','x[79]', 'x[80]', 'x[81]', 'x[82]');
% xlabel('time [minutes]');
% ylabel('state variables');
% title('Distillation column A states: x[63] - x[82]');


%clf()
%figure(1)
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
hold on;
if (nlpRun == 1)
%     plot(tgrid, x1_opt_pf, 'b-', tgrid, x1_opt_nlp, 'r--', 'LineWidth',2.5); hold on;
%     plot(tgrid, x2_opt_pf, 'b-', tgrid, x2_opt_nlp, 'r--', 'LineWidth',2.5); hold on;
%     plot(tgrid,x21_opt_pf, 'b-', tgrid, x21_opt_nlp,'r--', 'LineWidth',2.5); hold on;
%     plot(tgrid,x40_opt_pf, 'b-', tgrid, x40_opt_nlp,'r--', 'LineWidth',2.5); hold on;
%     plot(tgrid,x41_opt_pf, 'b-', tgrid, x41_opt_nlp,'r--', 'LineWidth',2.5);
%     legend('x[1]-reboiler_pf','x[1]-reboiler_nlp','x[2]-pf','x[2]-nlp','x[21]-feeder-pf','x[21]-feeder-nlp', 'x[40]-pf', 'x[40]-nlp', 'x[41]-condenser-pf', 'x[41]-condenser-nlp');
    hX1_pf   = plot(tgrid, x1_opt_pf); 
    set(hX1_pf,'Color','m','LineStyle','-','LineWidth',1.5,'Marker','o');
    hX1_nlp  = plot(tgrid, x1_opt_nlp); 
    set(hX1_nlp,'Color','m','LineStyle',':','LineWidth',1.5,'Marker','d');
    hX2_pf   = plot(tgrid, x2_opt_pf); 
    set(hX2_pf,'Color','b','LineStyle','-','LineWidth',1.5,'Marker','o');
    hX2_nlp  = plot(tgrid, x2_opt_nlp);
    set(hX2_nlp,'Color','b','LineStyle',':','LineWidth',1.5,'Marker','d');
    hX21_pf  = plot(tgrid, x21_opt_pf); 
    set(hX21_pf,'Color','g','LineStyle','-','LineWidth',1.5,'Marker','o');
    hX21_nlp = plot(tgrid, x21_opt_nlp); 
    set(hX21_nlp,'Color','g','LineStyle',':','LineWidth',1.5,'Marker','d');
    hX40_pf  = plot(tgrid, x40_opt_pf); 
    set(hX40_pf,'Color','r','LineStyle','-','LineWidth',1.5,'Marker','o');
    hX40_nlp = plot(tgrid, x40_opt_nlp);
    set(hX40_nlp,'Color','r','LineStyle',':','LineWidth',1.5,'Marker','d');
    hX41_pf  = plot(tgrid, x41_opt_pf); 
    set(hX41_pf,'Color','c','LineStyle','-','LineWidth',1.5,'Marker','o');
    hX41_nlp = plot(tgrid, x41_opt_nlp);
    set(hX41_nlp,'Color','c','LineStyle',':','LineWidth',1.5,'Marker','d');
    
else
    plot(tgrid,x1_opt_pf, 'LineWidth',2.5); hold on;
    plot(tgrid,x2_opt_pf, 'LineWidth',2.5); hold on;
    plot(tgrid,x21_opt_pf, 'LineWidth',2.5); hold on;
    plot(tgrid,x40_opt_pf, 'LineWidth',2.5); hold on;
    plot(tgrid,x41_opt_pf, 'LineWidth',2.5);
    legend('x[1]-reboiler','x[2]','x[21]-feeder','x[40]', 'x[41]-condenser');
end

%xlabel('time [minutes]');
%ylabel('state variables');
hTitle  = title ('State variables comparison PF - NLP');
hXLabel = xlabel('Time (second)'                      );
hYLabel = ylabel('Concentration (-)'                  );
hLegend = legend( ...
  [hX1_pf, hX1_nlp, hX2_pf, hX2_nlp, hX21_pf, hX21_nlp, hX40_pf, hX40_nlp, hX41_pf, hX41_nlp], ...
  'PF:  x[1]-reboiler' , ...
  'NLP: x[1]-reboiler' , ...
  'PF:  x[2]'          , ...
  'NLP: x[2]'          , ...
  'PF:  x[21]-feeder'  , ...
  'NLP: x[21]-feeder'  , ...
  'PF:  x[40]'         , ...
  'NLP: x[40]'         , ...
  'PF:  x[41]-condenser', ...
  'NLP: x[41]-condenser', ...
  'location', 'NorthWest' );

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:1.0, ...
  'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 5-states.eps
%close;

figure(2)
hold on;
if (nlpRun == 1)
    %plot(tgrid_u, u1_opt_pf, 'b-', tgrid_u,  u1_opt_nlp, 'r--', 'LineWidth',2.5); hold on;
    %plot(tgrid_u, u2_opt_pf, 'b-', tgrid_u,  u2_opt_nlp, 'r--', 'LineWidth',2.5);
    %legend('PF: u[1]-LT','NLP: u[1]-LT','PF: u[2]-VB', 'NLP: u[2]-VB');
    hU1_pf   = plot(tgrid_u, u1_opt_pf); 
    set(hU1_pf,'Color','m','LineStyle','-','LineWidth',1.5,'Marker','o');
    hU1_nlp  = plot(tgrid_u, u1_opt_nlp); 
    set(hU1_nlp,'Color','m','LineStyle',':','LineWidth',1.5,'Marker','d');
    hU2_pf   = plot(tgrid_u, u2_opt_pf); 
    set(hU2_pf,'Color','b','LineStyle','-','LineWidth',1.5,'Marker','o');
    hU2_nlp  = plot(tgrid_u, u2_opt_nlp);
    set(hU2_nlp,'Color','b','LineStyle',':','LineWidth',1.5,'Marker','d');
else
    plot(tgrid_u, u1_opt_pf, 'LineWidth',2.5); hold on;
    plot(tgrid_u, u2_opt_pf, 'LineWidth',2.5);
    legend('u[1]-LT','u[2]-VB');
end

% ylabel('control inputs');
% xlabel('time [minutes]');
hTitle1  = title ('Control inputs comparison PF - NLP');
hXLabel1 = xlabel('Time (second)'                      );
hYLabel1 = ylabel('Control inputs (-)'                  );
hLegend1 = legend( ...
  [hU1_pf, hU1_nlp, hU2_pf, hU2_nlp], ...
  'PF:  u[1]-LT' , ...
  'NLP: u[1]-LT' , ...
  'PF:  u[2]-VB'          , ...
  'NLP: u[2]-VB'          , ...
  'location', 'NorthWest' );

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.2:5.0, ...
  'LineWidth'   , 1         );

