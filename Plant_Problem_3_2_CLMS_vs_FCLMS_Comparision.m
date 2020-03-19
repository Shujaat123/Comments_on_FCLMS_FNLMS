clc;
clear all;
close all;

lw =1;          % line width in plot
f_size=10;      % font size of plot text
runs = 1000;     % number of independent simulations

N = 1e3;        % number of samples

% simulation parameters for Fractional variance of CLMS (FCLMS) [7], and CLMS
eta_FCLMS = 2e-2;   % step-size of FCLMS
f_FCLMS = [0.9  0.8  0.7  0.6  0.5  0.4]; % fractional powers

eta_CLMS = 4e-2; % step-size of CLMS

SNR = 10; % noise level

h = [-10:1:10]'; % system impulse response

% initialize MSE vectors
MSE_FCLMS1 = 0;
MSE_FCLMS2 = 0;
MSE_FCLMS3 = 0;
MSE_FCLMS4 = 0;
MSE_FCLMS5 = 0;
MSE_FCLMS6 = 0;
MSE_CLMS = 0;

for k = 1 : runs

x = complex(0.707*randn(1,N),0.707*randn(1,N));
r = filter(h,1,x);
d = awgn(r,SNR);

% initialize random weights
W_FCLMS1 = randn(size(h));
W_FCLMS2 = W_FCLMS1;
W_FCLMS3 = W_FCLMS1;
W_FCLMS4 = W_FCLMS1;
W_FCLMS5 = W_FCLMS1;
W_FCLMS6 = W_FCLMS1;
W_FCLMS7 = W_FCLMS1;

W_CLMS = W_FCLMS1;

% moving window for filteration
U = zeros(size(h));

    for n = 1 : N
        U(2:end,1) = U(1:end-1,1);
        U(1,1) = x(n);

        %% FCLMS
%% f1
        y_FCLMS1 = (W_FCLMS1')*U;
        e_FCLMS1 = d(n) - y_FCLMS1;      

        W_FCLMS1 = W_FCLMS1 +  eta_FCLMS*conj(e_FCLMS1)*U;
        W_FCLMS1 = W_FCLMS1 +  (eta_FCLMS/gamma(2-f_FCLMS(1)))*conj(e_FCLMS1)*U.*(W_FCLMS1.^(1-f_FCLMS(1)));
        
%% f2
        y_FCLMS2 = (W_FCLMS2')*U;
        e_FCLMS2 = d(n) - y_FCLMS2;      

        W_FCLMS2 = W_FCLMS2 +  eta_FCLMS*conj(e_FCLMS2)*U;
        W_FCLMS2 = W_FCLMS2 +  (eta_FCLMS/gamma(2-f_FCLMS(2)))*conj(e_FCLMS2)*U.*(W_FCLMS2.^(1-f_FCLMS(2)));
        
%% f3
        y_FCLMS3 = (W_FCLMS3')*U;
        e_FCLMS3 = d(n) - y_FCLMS3;      

        W_FCLMS3 = W_FCLMS3 +  eta_FCLMS*conj(e_FCLMS3)*U;
        W_FCLMS3 = W_FCLMS3 +  (eta_FCLMS/gamma(2-f_FCLMS(3)))*conj(e_FCLMS3)*U.*(W_FCLMS3.^(1-f_FCLMS(3)));
        
   %% f4
        y_FCLMS4 = (W_FCLMS4')*U;
        e_FCLMS4 = d(n) - y_FCLMS4;      

        W_FCLMS4 = W_FCLMS4 +  eta_FCLMS*conj(e_FCLMS4)*U;
        W_FCLMS4 = W_FCLMS4 +  (eta_FCLMS/gamma(2-f_FCLMS(4)))*conj(e_FCLMS4)*U.*(W_FCLMS4.^(1-f_FCLMS(4)));
        
   %% f5
        y_FCLMS5 = (W_FCLMS5')*U;
        e_FCLMS5 = d(n) - y_FCLMS5;      

        W_FCLMS5 = W_FCLMS5 +  eta_FCLMS*conj(e_FCLMS5)*U;
        W_FCLMS5 = W_FCLMS5 +  (eta_FCLMS/gamma(2-f_FCLMS(5)))*conj(e_FCLMS5)*U.*(W_FCLMS5.^(1-f_FCLMS(5)));

   %% f6
        y_FCLMS6 = (W_FCLMS6')*U;
        e_FCLMS6 = d(n) - y_FCLMS6;      

        W_FCLMS6 = W_FCLMS6 +  eta_FCLMS*conj(e_FCLMS6)*U;
        W_FCLMS6 = W_FCLMS6 +  (eta_FCLMS/gamma(2-f_FCLMS(6)))*conj(e_FCLMS6)*U.*(W_FCLMS6.^(1-f_FCLMS(6)));
        
        J_FCLMS1(n) = mean(abs(h-W_FCLMS1));
        J_FCLMS2(n) = mean(abs(h-W_FCLMS2));
        J_FCLMS3(n) = mean(abs(h-W_FCLMS3));
        J_FCLMS4(n) = mean(abs(h-W_FCLMS4));
        J_FCLMS5(n) = mean(abs(h-W_FCLMS5));
        J_FCLMS6(n) = mean(abs(h-W_FCLMS6));
        
        %% CLMS        
        y_CLMS = (W_CLMS')*U;
        e_CLMS = d(n) - y_CLMS;      
       
        W_CLMS = W_CLMS +  eta_CLMS*conj(e_CLMS)*U;
        
        J_CLMS(n) = mean(abs(h-W_CLMS));
    end
MSE_FCLMS1 = MSE_FCLMS1 + J_FCLMS1;    
MSE_FCLMS2 = MSE_FCLMS2 + J_FCLMS2;
MSE_FCLMS3 = MSE_FCLMS3 + J_FCLMS3;
MSE_FCLMS4 = MSE_FCLMS4 + J_FCLMS4;
MSE_FCLMS5 = MSE_FCLMS5 + J_FCLMS5;
MSE_FCLMS6 = MSE_FCLMS6 + J_FCLMS6;
MSE_CLMS = MSE_CLMS + J_CLMS;
end

MSE_FCLMS1 = MSE_FCLMS1./runs;
10*log10(MSE_FCLMS1(end))

MSE_FCLMS2 = MSE_FCLMS2./runs;
10*log10(MSE_FCLMS2(end))

MSE_FCLMS3 = MSE_FCLMS3./runs;
10*log10(MSE_FCLMS3(end))

MSE_FCLMS4 = MSE_FCLMS4./runs;
10*log10(MSE_FCLMS4(end))

MSE_FCLMS5 = MSE_FCLMS5./runs;
10*log10(MSE_FCLMS5(end))

MSE_FCLMS6 = MSE_FCLMS6./runs;
10*log10(MSE_FCLMS6(end))


MSE_CLMS = MSE_CLMS./runs;
10*log10(MSE_CLMS(end))

% Plot Results
figure
plot(10*log10(MSE_FCLMS1),'linewidth',lw)
hold on
plot(10*log10(MSE_FCLMS2),'linewidth',lw)
plot(10*log10(MSE_FCLMS3),'linewidth',lw)
plot(10*log10(MSE_FCLMS4),'linewidth',lw)
plot(10*log10(MSE_FCLMS5),'linewidth',lw)
plot(10*log10(MSE_FCLMS6),'linewidth',lw)
plot(10*log10(MSE_CLMS),'b','linewidth',lw)
ylim([-(SNR+10) 20])
xlabel('No. of iterations')
ylabel('\DeltaW (dB)')
grid minor
set(gca,'FontSize',f_size)
legend('Fractional-order variant of Complex LMS (f=0.9)','Fractional-order variant of Complex LMS (f=0.8)','Fractional-order variant of Complex LMS (f=0.7)','Fractional-order variant of Complex LMS (f=0.6)','Fractional-order variant of Complex LMS (f=0.5)','Fractional-order variant of Complex LMS (f=0.4)','Complex LMS / Fractional-order variant of Complex LMS (f=1)');
title('Problem # 3.2')
