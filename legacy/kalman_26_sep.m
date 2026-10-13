clc;
clear all;
close all;


%synthetic data
dt = 2/1000; f0 = 20; tmax = 3; 
h = 20:20:1000;
tau = [0.5,0.9,1.5,1.7,2]; 
v = [800,1200,1600,2000,3000]; 
amp = [3,3,3,3,3]; 
snr = 1000; L = 1000;

d = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
data =[fliplr(d),d]; % synthetic data

%%
y =0.04*randn(size(data)); % random noise
data_Noisy=data+y;% noisy data

Filter_sav=sgolayfilt(data_Noisy,4,27); % savitzky filter

filter_fx=fx_decon(data_Noisy,.001,8,0,0,200); % FX filter

[u,s,v]=svd(data_Noisy);  %svd method

R=rank(data_Noisy);

v1=v';
Image=zeros(1501,R,18); %eigenimage
for i=1:18;
    Image(:,:,i)=s(i,i)*u(:,i)*v1(i,:);
end
DD=sum(Image,3); %sum of eigenimage and denoising with svd
%%
% Kalman
q = 0.004;

xf=zeros(1501,100);
kk=zeros(1501,100);
pp=zeros(1501,100);
for j=1:100
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(data_Noisy(:,j),q); % kalman filter
end

snr=10*log10(sum(sum(data.^2))/sum(sum(y.^2))) %snr of noisy data
s1=Filter_sav-data;
snrsg=10*log10(sum(sum(Filter_sav.^2))/sum(sum(s1.^2))) %snr of savitzky filter
snrfx=10*log10(sum(sum(filter_fx.^2))/sum(sum((filter_fx-data).^2))) % snr of FX filter
snrsvd=10*log10((sum(sum(DD.^2))/sum(sum((DD-data).^2)))) % snr of svd

res=xf-data;
snrKalman=10*log10(sum(sum(xf.^2))/sum(sum(res.^2))) % snr of kalman

%%
% VRE
xfff=VRE(xf,2,10); 
resss=xfff-data;
snrVRE=10*log10(sum(sum(xfff.^2))/sum(sum(resss.^2)))

%%
figure
subplot(121)
wiggle(data)
title('Data (synthetic)');
xlabel('Number of traces');
ylabel('Time (ms)');

subplot(122)
wiggle(data_Noisy)
title('Noisy ');
xlabel('Number of traces');
ylabel('Time (ms)');

%

figure
subplot(221)
wiggle(Filter_sav);
title('filtered S-G ');
xlabel('number of traces');
ylabel('time')


subplot(222)
wiggle(filter_fx)
title(' FX')
xlabel('Number of traces');
ylabel('Time (ms)');




subplot(223)
wiggle(DD)
title('svd');

subplot(224)
wiggle(xf);
title(' Kalman')
xlabel('Number of traces');
ylabel('Time (ms)');



% figure
% subplot(121)
% wiggle(xf);
% title('filtered Kalman')
% xlabel('number of traces');
% ylabel('time');
% 
% 
% 
% subplot(122)
% wiggle(xfff)
% title('filtered Kalman with VRE')
% xlabel('number of traces');
% ylabel('time');

figure;
imagesc(data_Noisy)
colorbar
title('Data ');
xlabel('Offset');
ylabel('Time (ms)');

figure;
imagesc(xf)
title('Kalman ');
xlabel('Number of traces');
ylabel('Time (ms)');
colorbar



%%



[nt nx] = size(data_Noisy);

dt1 = 0.001;
fs=1/dt1;
% N=l000;
f = [-nt/2:nt/2-1]*fs/1000;

dd2 = fftshift(abs(fft(data(:,:))));

for i = 1 : 1 : nt
f_d2(i,:) = max(dd2(i,:));
end

figure
subplot(331);
plot(f,(f_d2) ,'LineWidth',1.5,'color','k')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])

title('Clean Data')

%

dd = fftshift(abs(fft(data_Noisy(:,:))));

for i = 1 : 1 : nt
f_dd(i,:) = max(dd(i,:));
end

subplot(332);
plot(f,(f_dd) ,'LineWidth',1.5,'color','b')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])

title('Noisy Data')


mm2 = fftshift(abs(fft(Filter_sav(:,:))));

for i = 1 : 1 : nt
myfil(i,:) = max(mm2(i,:));
end
subplot(333);
plot(f,(myfil) ,'LineWidth',1.5,'color','b')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])
title('savitzky filter')



jj=fftshift(abs(fft(filter_fx(:,:))));
for i = 1 : 1 : nt
fk(i,:) = max(jj(i,:));
end
subplot(334);
plot(f,(fk) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])

title('FX filter')




F = fftshift(abs(fft(DD(:,:))));

for i = 1 : 1 : nt
fk(i,:) = max(F(i,:));
end

subplot(335);
plot(f,(fk) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])

title('SVD method')


mm = fftshift(abs(fft(xf(:,:))));

for i = 1 : 1 : nt
myfil(i,:) = max(mm(i,:));
end
subplot(336);
plot(f,(myfil) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])
title('Kalman')

%qq=fftshift(abs(fft(xfff(:,:))));
%for i = 1 : 1 : nt
%myfil(i,:) = max(qq(i,:));
%end
%subplot(337);
%plot(f,normalize(myfil) ,'LineWidth',1.5,'color','r')
%xlabel('Frequency (Hz)')
%ylabel('Amplitude')
%xlim([0 250])
% ylim([0 220])
%title('Kalman with VRE')

% 
% sss=fftshift(abs(fft(y(:,:))));
% for i = 1 : 1 : nt
% myfil(i,:) = max(sss(i,:));
% end
subplot(338);
plot(f,(myfil) ,'LineWidth',1.5,'color','g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
ylim([0 150])
title('random noise')
%%
figure
plot(f,(f_d2) ,'LineWidth',1.5,'color','k')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

hold on


plot(f,(f_dd) ,'LineWidth',1.5,'color','b')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

hold on

plot(f,(fk) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

hold on

plot(f,(myfil) ,'LineWidth',1.5,'color','g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

legend('Clean Data','Noisy Data','Kalman')
%%
  
[S,k,f] = fk_spectra(data_Noisy,0.001,12.5,6);

figure
imagesc(k,fliplr(f),flipud(S));
title("FK Data")  
ylim([0 150])
xlabel("Wavenumber (1/km)")
ylabel("Frequency (Hz)")

[S,k,f] = fk_spectra(xf,0.001,12.5,6);
  
figure
imagesc(k,fliplr(f),flipud(S));
title("FK Kalman")
ylim([0 150])
xlabel("Wavenumber (1/km)")
ylabel("Frequency (Hz)")
%%

n = 50;
subplot(311)
plot(data_Noisy(:,n))
ylim([-3.5 3.5])
xlabel("Time (ms)")
ylabel("Amplitude (m)")
title("Data")

subplot(312)
plot(filter_fx(:,n))
ylim([-3.5 3.5])
xlabel("Time (ms)")
ylabel("Amplitude (m)")
title("FX")


subplot(313)
plot(xf(:,n))
ylim([-3.5 3.5])
xlabel("Time (ms)")
ylabel("Amplitude (m)")
title("Kalman")












