clc
clear all
close all
lw=1;
f_size=10;
load('results\Results_Problem3_2_1.mat')

figure
plot(10*log10(MSE_FCLMS1),'.-','linewidth',lw)
hold on
plot(10*log10(MSE_FCLMS2),'--','linewidth',lw)
plot(10*log10(MSE_FCLMS3),'o-','linewidth',lw,'MarkerIndices',1:10:length(MSE_FCLMS3))
plot(10*log10(MSE_FCLMS4),'s-','linewidth',lw,'MarkerIndices',1:10:length(MSE_FCLMS4))
plot(10*log10(MSE_FCLMS5),'+-','linewidth',lw)
plot(10*log10(MSE_FCLMS6),'d-','linewidth',lw)
plot(10*log10(MSE_CLMS),':.','linewidth',lw)
ylim([-15 40])
xlim([0 450])
xlabel('No. of iterations')
ylabel('\DeltaW (dB)')
grid minor
set(gca,'FontSize',f_size)
h_leg=legend('FCLMS (f=0.9)','FCLMS (f=0.8)','FCLMS (f=0.7)','FCLMS (f=0.6)','FCLMS (f=0.5)','FCLMS (f=0.4)','CLMS / FCLMS (f=1)','location','northeast');
% title('Problem # 3.2.1')
grid minor
set(gca,'FontSize',f_size)
set(h_leg,'FontSize',f_size)

saveas(gcf,strcat('figures\Problem_3_2_1.png'),'png')
saveas(gcf,strcat('figures\Problem_3_2_1.eps'),'psc2')

