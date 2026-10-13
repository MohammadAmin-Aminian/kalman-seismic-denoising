clc;
clear;
close all;
%%
% Synthetic data
dt = 2/1000; f0 = 20; tmax = 3; 
h = 20:20:1000;
tau = [0.5,0.9,1.5,1.7,2]; 
v = [800,1200,1600,2000,3000]; 
amp = [3,3,3,3,3]; 
% snr = 1000000000000 ; L = 10000;
snr = 1000; L = 1000;

d = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
A1 =[fliplr(d),d];

y = 0.4*randn(size(A1)); % random noise
A2=A1+y;% noisy data


% dt = 0.001;
% f0 = 25;
% tmax = 5;
% x = 0:5:500;
% tau = [.5,1.5,2]; 
% v = [400,420,440]; 
% amp = [5,-5,5]; 
% snr = 1000;                    
% L = 1000;
% A1 = 5*hyperbolic_events(dt,f0,tmax,x,tau,v,amp,snr,L); %synthetic data

%%
A3=sgolayfilt(A2,4,31); % savitzky filter

figure(1)
subplot(1,3,1)
wiggle(A1)
title('synthetic');
xlabel('number of traces');
ylabel('time');

subplot(1,3,2)
wiggle(A2)
title('noisy snr=4');
xlabel('number of traces');
ylabel('time');

subplot(1,3,3)
wiggle(A3);
title('filtered S-G snr=13');
xlabel('number of traces');
ylabel('time')
%%

A4=fx_decon(A2,.001,4,1,1,500); % FX filter
figure(4)
wiggle(A4,dt)
title('filtered FX snr=10')
xlabel('number of traces');
ylabel('time');

snr=10*log10(sum(sum(A1.^2))/sum(sum(y.^2))) %snr of noisy data
s1=A3-A1;
snrsg=10*log10(sum(sum(A3.^2))/sum(sum(s1.^2))) %snr of savitzky filter
s2=A4-A1;
snrfx=10*log10(sum(sum(A4.^2))/sum(sum(s2.^2))) % snr of FX filter

%%

[u,s,v]=svd(A2);  %svd method

R=rank(A2);

v1=v';
Image=zeros(length(A1),R,50); %eigenimage
for i=1:50;
    Image(:,:,i)=s(i,i)*u(:,i)*v1(i,:);
end
DD=sum(Image,3); %sum of eigenimage and denoising with svd

figure
wiggle(DD)
title(' filtered svd snr=6');

snrsvd=10*log10((sum(sum(A1.^2))/sum(sum((DD-A1).^2)))) % snr of svd

%%
% Kalman

q = 0.005;
xf=zeros(length(A1),width(A1));
kk=zeros(length(A1),width(A1));
pp=zeros(length(A1),width(A1));
for j=1:width(A1)
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(A2(:,j),q); % kalman filter
end

figure
wiggle(xf);
title('filtered Kalman snr=7')
xlabel('number of traces');
ylabel('time');
res = xf-A1;

snrKalman = 10*log10(sum(sum(xf.^2))/sum(sum(res.^2))); % snr of kalman
%%

[nt nx] = size(A2);

dt1 = 0.001;
fs=1/dt1;
% N=l000;
f = [-nt/2:nt/2-1]*fs/1000;

dd2 = fftshift(abs(fft(A1(:,:))));

for i = 1 : 1 : nt
f_d2(i,:) = max(dd2(i,:));
end

subplot(231);
plot(f,normalize(f_d2,'range') ,"LineWidth",1.5,'color','k')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
xlim([0 250])
% ylim([0 220])

title("Clean Data")

%

dd = fftshift(abs(fft(A2(:,:))));

for i = 1 : 1 : nt
f_dd(i,:) = max(dd(i,:));
end

subplot(232);
plot(f,normalize(f_dd,'range') ,"LineWidth",1.5,'color','b')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
xlim([0 250])
% ylim([0 220])

title("Noisy Data")


F = fftshift(abs(fft(DD(:,:))));

for i = 1 : 1 : nt
fk(i,:) = max(F(i,:));
end

subplot(233);
plot(f,normalize(fk,'range') ,"LineWidth",1.5,'color','k')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
xlim([0 250])
% ylim([0 220])

title("SVD method")

mm2 = fftshift(abs(fft(A3(:,:))));

for i = 1 : 1 : nt
myfil(i,:) = max(mm2(i,:));
end
subplot(234);
plot(f,normalize(myfil,'range') ,"LineWidth",1.5,'color','r')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
xlim([0 250])
% ylim([0 220])
title("savitzky filter")

mm = fftshift(abs(fft(xf(:,:))));

for i = 1 : 1 : nt
myfil(i,:) = max(mm(i,:));
end
subplot(235);
plot(f,normalize(myfil,'range') ,"LineWidth",1.5,'color','r')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
xlim([0 250])
% ylim([0 220])
title("Kalman")
