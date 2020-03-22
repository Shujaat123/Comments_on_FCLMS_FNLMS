clc
clear all
close all
lw=1;
f_size=10;
load('results\Results_Problem3_1_2.mat')

figure
plot(10*log10(MSE_FNLMS1),'.-','linewidth',lw)
hold on
plot(10*log10(MSE_FNLMS2),'--','linewidth',lw)
plot(10*log10(MSE_FNLMS3),'o-','linewidth',lw)
plot(10*log10(MSE_FNLMS4),'s-','linewidth',lw)
plot(10*log10(MSE_FNLMS5),'+-','linewidth',lw)
plot(10*log10(MSE_FNLMS6),'d-','linewidth',lw)
plot(10*log10(MSE_NLMS2),':.','linewidth',lw)
ylim([-10 80])
xlim([0 90])
xlabel('No. of iterations')
ylabel('\DeltaW (dB)')
grid minor
set(gca,'FontSize',f_size)
h_leg=legend('FNLMS (f=0.9)','FNLMS (f=0.8)','FNLMS (f=0.7)','FNLMS (f=0.6)','FNLMS (f=0.5)','FNLMS (f=0.4)','NLMS','location','northeast');
% title('Problem # 3.1.2')
grid minor
set(gca,'FontSize',f_size)
set(h_leg,'FontSize',f_size)

saveas(gcf,strcat('figures\Problem_3_1_2a.png'),'png')
saveas(gcf,strcat('figures\Problem_3_1_2a.eps'),'psc2')


figure
plot(10*log10(MSE_FNLMS1),'.-','linewidth',lw)
hold on
plot(10*log10(MSE_FNLMS2),'--','linewidth',lw)
plot(10*log10(MSE_FNLMS3),'o-','linewidth',lw)
plot(10*log10(MSE_FNLMS4),'s-','linewidth',lw)
plot(10*log10(MSE_FNLMS5),'+-','linewidth',lw)
plot(10*log10(MSE_FNLMS6),'d-','linewidth',lw)
plot(10*log10(MSE_NLMS2),':.','linewidth',lw)
ylim([-10 80])
xlim([0 90])
xlabel('No. of iterations')
ylabel('\DeltaW (dB)')
grid minor
set(gca,'FontSize',f_size)
h_leg=legend('FNLMS (f=0.9)','FNLMS (f=0.8)','FNLMS (f=0.7)','FNLMS (f=0.6)','FNLMS (f=0.5)','FNLMS (f=0.4)','NLMS / FNLMS (f=1)','location','northeast');
% title('Problem # 3.1.2')
grid minor
set(gca,'FontSize',f_size)
set(h_leg,'FontSize',f_size)

saveas(gcf,strcat('figures\Problem_3_1_2b.png'),'png')
saveas(gcf,strcat('figures\Problem_3_1_2b.eps'),'psc2')

