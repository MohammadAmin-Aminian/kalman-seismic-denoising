clc;clear;close all;
dt = 0.001;
f0 = 25;
tmax = 5;

x = 0:5:500;

tau = [.5,1.5,2]; 
v = [400,420,440]; 
amp = [5,-5,5]; 
snr = 1000;                    
L = 1000;


A1 = 5*hyperbolic_events(dt,f0,tmax,x,tau,v,amp,snr,L); %synthetic data
y=.08*randn(size(A1)); % random noise
A2=A1+y;% noisy data
A3=sgolayfilt(A2,4,31); % savitzky filter

figure(1)
plot_traces_normal(A1,dt)
title('synthetic');
xlabel('number of traces');
ylabel('time');
figure(2)
plot_traces_normal(A2,dt)
title('noisy snr=4');
xlabel('number of traces');
ylabel('time');


figure(3)
plot_traces_normal(A3,dt);
title('filtered S-G snr=13');
xlabel('number of traces');
ylabel('time')

A4=fx_decon(A2,.001,4,1,1,500); % FX filter
figure(4)
plot_traces_normal(A4,dt)
title('filtered FX snr=10')
xlabel('number of traces');
ylabel('time');

snr=10*log10(sum(sum(A1.^2))/sum(sum(y.^2))) %snr of noisy data
s1=A3-A1;
snrsg=10*log10(sum(sum(A3.^2))/sum(sum(s1.^2))) %snr of savitzky filter
s2=A4-A1;
snrfx=10*log10(sum(sum(A4.^2))/sum(sum(s2.^2))) % snr of FX filter



[u,s,v]=svd(A2);  %svd method

R=rank(A2);

v1=v';
Image=zeros(5001,R,50); %eigenimage
for i=1:50;
    Image(:,:,i)=s(i,i)*u(:,i)*v1(i,:);
end
DD=sum(Image,3); %sum of eigenimage and denoising with svd

figure
plot_traces_normal(DD,dt)
title(' filtered svd snr=6');

snrsvd=10*log10((sum(sum(A1.^2))/sum(sum((DD-A1).^2)))) % snr of svd


q=.005;

xf=zeros(5001,101);
kk=zeros(5001,101);
pp=zeros(5001,101);
for j=1:101
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(A2(:,j),q); % kalman filter
end

figure
plot_traces_normal(xf,dt);
title('filtered Kalman snr=7')
xlabel('number of traces');
ylabel('time');
res=xf-A1;
snrKalman=10*log10(sum(sum(xf.^2))/sum(sum(res.^2))) % snr of kalman

