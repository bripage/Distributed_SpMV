clear; clc;
Nfiles = 25;
PHI = zeros(Nfiles,1);
Ms = [1270432 525825 643994 1219574 259789 147900 1062400 221119 589446 415863 227632 245874 79171 71505 83334 155331 116158 1102824 381689 97578 63838 72000 38120 45101 22283];
Ps = [4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64 4 9 16 25 36 49 64];
Ts = zeros(length(Ps), 1);

for i = 1:Nfiles
% load experimental data
fileID = fopen('Mval.txt', 'w');
fprintf(fileID, '%f', Ms(i));
fclose all;

str = ['m',num2str(i),'.txt'];
EXP = load(str);
P_EXP = EXP(:,1);
SDR_EXP = EXP(:,2);

for j = 1:7
    Ts(((i-1)*7)+j) = SDR_EXP(j);
end

ft = fittype('SDR(x,phi)'); 

phig = 1;
f = fit( P_EXP, SDR_EXP, ft, 'StartPoint', [phig], 'TolFun', 1e-20, 'TolX', 1e-20, 'MaxFunEvals', 10000);
PHI(i) = f.phi;

end

d = 64; %size of double in bits
B = 60129542144; %network bandwidth in bits
PHI_mean = mean(PHI);
t_bar = zeros(length(Ts),1);
for i = 1:length(Ms)
    for j = 1:7
        t_bar(((i-1)*7)+j) = (B/(((Ms(i)/sqrt(Ps(i)))*d)+PHI_mean));
    end
end

for i = 1:length(Ts)
        Ts(i) = Ts(i) * t_bar(i);
end

x = [1 4 9 16 25 36 49 64];
y = [2 4 9.51 16 23.22 31.02 39.302 48];
pl = plot(x, y, 'k', Ps, Ts , 'ko');

set(pl, 'LineWidth', 2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
xlabel('\# of MPI Processes', 'interpreter', 'latex', 'FontSize', 18);
ylabel('Time', 'interpreter', 'latex', 'FontSize', 18);
title('Normalized MPI Reduction Time', 'FontSize', 24);
print('timeplot_naive', '-dpdf');

