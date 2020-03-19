clc;
clear all;
close all;

lw = 1;     % linewidth of plot
N = 1e3;    % number of samples
runs = 1000; % number of independent simulations
f_size = 10; % font size of graphs

eta_NLMS = 10e-1; % step-size of NLMS
eta_FNLMS = 5e-1; % step-size of FNLMS
f_FNLMS = [0.9 0.8 0.7 0.6 0.5 0.4]; % fractional powers

h = [1 2 2 2 1 1 2 2 3 1 1 2 2 2 1 2 1 2 2 2 1 1 2 2 2 1 1 3 2 2]'; % system impulse response

for k = 1 : runs

x = randn(1,N);
d = filter(h,1,x);

% W_FNLMS1 = randn(size(h));
W_FNLMS1 = zeros(size(h));
W_FNLMS2 = W_FNLMS1;
W_FNLMS3 = W_FNLMS1;
W_FNLMS4 = W_FNLMS1;
W_FNLMS5 = W_FNLMS1;
W_FNLMS6 = W_FNLMS1;
W_FNLMS7 = W_FNLMS1;
U = zeros(size(h));

W_NLMS = W_FNLMS1;

    for n = 1 : N
        U(2:end,1) = U(1:end-1,1);
        U(1,1) = x(n);

%% f1
        y_FNLMS1 = (W_FNLMS1')*U;
        e_FNLMS1 = d(n) - y_FNLMS1;      

        W_FNLMS1 = W_FNLMS1 +  eta_FNLMS*(e_FNLMS1)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS1 = W_FNLMS1 +  (eta_FNLMS/gamma(2-f_FNLMS(1)))*(e_FNLMS1)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS1.^(1-f_FNLMS(1)));
        
%% f2
        y_FNLMS2 = (W_FNLMS2')*U;
        e_FNLMS2 = d(n) - y_FNLMS2;      

        W_FNLMS2 = W_FNLMS2 +  eta_FNLMS*(e_FNLMS2)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS2 = W_FNLMS2 +  (eta_FNLMS/gamma(2-f_FNLMS(2)))*(e_FNLMS2)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS2.^(1-f_FNLMS(2)));
        
%% f3
        y_FNLMS3 = (W_FNLMS3')*U;
        e_FNLMS3 = d(n) - y_FNLMS3;      

        W_FNLMS3 = W_FNLMS3 +  eta_FNLMS*(e_FNLMS3)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS3 = W_FNLMS3 +  (eta_FNLMS/gamma(2-f_FNLMS(3)))*(e_FNLMS3)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS3.^(1-f_FNLMS(3)));
        
   %% f4
        y_FNLMS4 = (W_FNLMS4')*U;
        e_FNLMS4 = d(n) - y_FNLMS4;      

        W_FNLMS4 = W_FNLMS4 +  eta_FNLMS*(e_FNLMS4)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS4 = W_FNLMS4 +  (eta_FNLMS/gamma(2-f_FNLMS(4)))*(e_FNLMS4)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS4.^(1-f_FNLMS(4)));
        
   %% f5
        y_FNLMS5 = (W_FNLMS5')*U;
        e_FNLMS5 = d(n) - y_FNLMS5;      

        W_FNLMS5 = W_FNLMS5 +  eta_FNLMS*(e_FNLMS5)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS5 = W_FNLMS5 +  (eta_FNLMS/gamma(2-f_FNLMS(5)))*(e_FNLMS5)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS5.^(1-f_FNLMS(5)));

   %% f6
        y_FNLMS6 = (W_FNLMS6')*U;
        e_FNLMS6 = d(n) - y_FNLMS6;      

        W_FNLMS6 = W_FNLMS6 +  eta_FNLMS*(e_FNLMS6)*U.*(1/(norm(U)^2+1e-10));
        W_FNLMS6 = W_FNLMS6 +  (eta_FNLMS/gamma(2-f_FNLMS(6)))*(e_FNLMS6)*U.*(1/(norm(U)^2+1e-10)).*(W_FNLMS6.^(1-f_FNLMS(6)));
        
        J_FNLMS1(k,n) = mean(abs(h-W_FNLMS1));
        J_FNLMS2(k,n) = mean(abs(h-W_FNLMS2));
        J_FNLMS3(k,n) = mean(abs(h-W_FNLMS3));
        J_FNLMS4(k,n) = mean(abs(h-W_FNLMS4));
        J_FNLMS5(k,n) = mean(abs(h-W_FNLMS5));
        J_FNLMS6(k,n) = mean(abs(h-W_FNLMS6));
    
   %% NLMS
        y_NLMS = (W_NLMS')*U;
        e_NLMS = d(n) - y_NLMS;      
        
        W_NLMS = W_NLMS +  eta_NLMS*(e_NLMS)*U.*(1/(norm(U)^2+1e-10));
        
        J_NLMS(k,n) = mean(abs(h-W_NLMS));
    end
end

MSE_FNLMS1 = mean(J_FNLMS1,1);
10*log10(MSE_FNLMS1(end))

MSE_FNLMS2 = mean(J_FNLMS2,1);
10*log10(MSE_FNLMS2(end))

MSE_FNLMS3 = mean(J_FNLMS3,1);
10*log10(MSE_FNLMS3(end))

MSE_FNLMS4 = mean(J_FNLMS4,1);
10*log10(MSE_FNLMS4(end))

MSE_FNLMS5 = mean(J_FNLMS5,1);
10*log10(MSE_FNLMS5(end))

MSE_FNLMS6 = mean(J_FNLMS6,1);
10*log10(MSE_FNLMS6(end))

MSE_NLMS = mean(J_NLMS,1);
10*log10(MSE_NLMS(end))

% Plot Results
figure
plot(10*log10(MSE_FNLMS1),'linewidth',lw)
hold on
plot(10*log10(MSE_FNLMS2),'linewidth',lw)
plot(10*log10(MSE_FNLMS3),'linewidth',lw)
plot(10*log10(MSE_FNLMS4),'linewidth',lw)
plot(10*log10(MSE_FNLMS5),'linewidth',lw)
plot(10*log10(MSE_FNLMS6),'linewidth',lw)
plot(10*log10(MSE_NLMS),'b','linewidth',lw)
ylim([-10 80])
xlim([0 120])
xlabel('No. of iterations')
ylabel('\DeltaW (dB)')
grid minor
set(gca,'FontSize',f_size)
legend('Fractional-order variant of NLMS (f=0.9)','Fractional-order variant of NLMS (f=0.8)','Fractional-order variant of NLMS (f=0.7)','Fractional-order variant of NLMS (f=0.6)','Fractional-order variant of NLMS (f=0.5)','Fractional-order variant of NLMS (f=0.4)','NLMS / Fractional-order variant of NLMS (f=1)');
title('Problem # 3.1.2')


