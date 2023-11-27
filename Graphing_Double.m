clear
tic

cd '/Users/aishiguha/Documents/MATLAB/Code'

r = 0.9;
gamma2 = 2.3;
gamma1 = (0:0.1:3)'; 

F_d = [];
B_d = [];


cd '/Users/aishiguha/Documents/MATLAB/Data'

trials = 10;

for i = 1:trials
    filename = 'Double_r-'+string(r) +'_g2-' +string(gamma2) + '_g1-'+string(gamma1(length(gamma1)))+','+string(gamma1(2)-gamma1(1))+','+string(gamma1(1)) + '_trial-' + string(i) +'.mat';
    load(filename, 'Fidelity','BellStat')
    F_d = [F_d ; Fidelity];
    B_d = [B_d ; BellStat];
end


%Fidelity
% Mean_Fd = mean(F_d);
% sd_Fd = std(F_d);
% Upper_Fd = Mean_Fd + (sd_Fd)/(sqrt(trials));
% Lower_Fd = Mean_Fd - (sd_Fd)/(sqrt(trials));
% 
% [Max_Fd, index] = max(Mean_Fd);
% g1_max = gamma1(index);

%Bell Stat
Mean_Bd = mean(B_d);
sd_Bd = std(B_d);
Upper_Bd = Mean_Bd + (sd_Bd)/(sqrt(trials));
Lower_Bd = Mean_Bd - (sd_Bd)/(sqrt(trials));

[Max_Bd, index] = max(Mean_Bd);
g1_max = gamma1(index);

% hold on
% plot(gamma1, Mean_Fd, 'ro-')
% plot(gamma1, Upper_Fd, 'b-')
% plot(gamma1, Lower_Fd, 'b-')
% 
% plot(g1_max,Max_Fd,'r.', 'MarkerSize',30)
% xlabel('$\gamma_{\rm A}$','Interpreter','LaTeX','FontSize',20)
% ylabel('Average Fidelity')
% hold off
% 
% disp(Max_Fd)
% disp(Upper_Fd(index) - Max_Fd)
% 
% figname = 'FigDouble-Fidvg1_r-' + string(r) + '_g2-' + string(gamma2) + '.fig';


% 

hold on
plot(gamma1, Mean_Bd, 'bo-')
plot(gamma1, Upper_Bd, 'r-')
plot(gamma1, Lower_Bd, 'r-')

%plot(g1_max,Max_Bd,'r.', 'MarkerSize',30)
xlabel('$\gamma_{\rm A}$','Interpreter','LaTeX','FontSize',20)
ylabel('Average Bell Statistic')

%horizontal line
yline(2,'--');
yline(2*sqrt(2),'--')


hold off

disp(Max_Bd)
disp(Upper_Bd(index) - Max_Bd)
disp(g1_max)

figname = 'FigDouble-BellStatvg1_r-' + string(r) + '_g2-' + string(gamma2) + '.fig';


% 


cd '/Users/aishiguha/Documents/MATLAB/Figures'


saveas(gcf, figname)


cd '/Users/aishiguha/Documents/MATLAB/Code'

toc