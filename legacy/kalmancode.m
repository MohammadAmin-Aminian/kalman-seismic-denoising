clc;close;clear;
load('Sigs.mat')

dt = 1/40;
t = 0:dt:length(sig1)*dt-dt;  
xt = (1:length(sig1))*0.025;
xf = (-floor(length(sig1)/2):ceil(length(sig1)/2)-1)'/0.025/length(sig1);

kal=kalman(sig1,5e-5);

subplot(221)
plot(t,sig1);
title('Before kalman filter')
xlabel('Time(sec)')
ylabel('Amplitude')
ylim([-3e4 4e4])
xlim([0 660])
subplot(223)
plot(t,kal);
ylim([-3e4 4e4])
xlim([0 660])
title('After kalman filter')
xlabel('Time(sec)')
ylabel('Amplitude')
%
fmax=1./dt;
f=-fmax/2:(1./max(t)):fmax/2;
fsig=fftshift(fft(sig1));
fkal=fftshift(fft(kal));
subplot(222)
plot(f,fsig,'r');
xlim([-5 5])
title('Before kalman filter')
xlabel('Frequncy(hz)')
ylabel('Amplitude')
subplot(224)
plot(f,fkal,'r');
title('After kalman filter')
xlabel('Frequncy(hz)')
ylabel('Amplitude')
xlim([-5 5])

%%
close
subplot(121)
plot(t,sig1-kal');
title('Difference in time')
xlabel('Time(sec)')
ylabel('Amplitude')
ylim([-3e4 4e4])
xlim([0 660])
subplot(122)
plot(f,fsig-fkal','r');
xlim([-5 5])
title('Difference in Frequency')
xlabel('Frequncy(hz)')
ylabel('Amplitude')