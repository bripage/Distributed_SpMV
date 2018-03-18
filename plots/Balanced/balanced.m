clear; clc;
Nfiles = 10;
PHI = zeros(Nfiles,1);
B = zeros(Nfiles,1);
Ms = [1270432 525825 643994 1219574 259789 147900 1062400 221119 589446 415863 227632 245874 79171 71505 83334 155331 116158 1102824 381689 97578 63838 72000 38120 45101 22283];

for i = 1:Nfiles
% load experimental data
fileID = fopen('Mval.txt', 'w');
fprintf(fileID, '%f', Ms(i));
fclose all;

str = ['m',num2str(i),'l.txt'];
EXP = load(str);
P_EXP_L = EXP(:,1);
BDR_EXP_L = EXP(:,2);

str = ['m',num2str(i),'r.txt'];
EXP = load(str);
P_EXP_R = EXP(:,1);
BDR_EXP_R = EXP(:,2);

ftL = fittype('BDRL(x,B,phi)'); 
ftR = fittype('BDRR(x,B,phi)'); 

Bg = 60129542144;
phig = 50000000;
fL = fit( P_EXP_L, BDR_EXP_L, ftL, 'StartPoint', [Bg phig], 'TolFun', 1e-20, 'TolX', 1e-20, 'MaxFunEvals', 1000, 'Exclude', P_EXP_L > 64);
fR = fit( P_EXP_R, BDR_EXP_R, ftR, 'StartPoint', [Bg phig], 'TolFun', 1e-20, 'TolX', 1e-20, 'MaxFunEvals', 1000, 'Exclude', P_EXP_R < 81);
PHI_L = fL.phi;
PHI_R = fR.phi;
B_L = fL.B;
B_R = fR.B;

%p1 = plot( f, 'k', P_EXP, BDR_EXP , 'ko');
%set(p1, 'LineWidth', 2);
%set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
%xlabel('\# of MPI Processes', 'interpreter', 'latex', 'FontSize', 18);
%ylabel('Time', 'interpreter', 'latex', 'FontSize', 18);
%title('MPI Reduction Time', 'FontSize', 24);
%matrixFig = ['m',num2str(i)];
%print(matrixFig, '-dpdf');


d = 64; %size of double in bits
t_barL = (B_L*BDR_EXP_L./(((Ms(i)./P_EXP_L)*d)+PHI_L));
t_barR = (B_R*BDR_EXP_R./(((Ms(i)./P_EXP_R)*d)+PHI_R));

pl = plot(P_EXP_L,t_barL, 'o');
hold on;
pl = plot(P_EXP_R,t_barR, 'o');
hold on;

set(pl, 'LineWidth', 2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
xlabel('\# of MPI Processes', 'interpreter', 'latex', 'FontSize', 18);
ylabel('Time', 'interpreter', 'latex', 'FontSize', 18);
title('Normalized MPI Reduction Time', 'FontSize', 24);

end

x1 = linspace(1,100,100)';
y1 = log2(x1);
x2 = linspace(81,121,100)';
y2 = 1.5*log2(x2-81);
xlim([0 130]);
ylim([1 12]);
plot(x1, y1, 'k');
plot(x2, y2, 'k');
%print('timeplot_naive', '-dpdf');
