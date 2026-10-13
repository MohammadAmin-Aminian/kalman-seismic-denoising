[nt nx] = size(A1);

% Chirp + Gaussian noise 40%
t = 0 : 0.01 : 15.00;
t1 = 0.05;
f1 = 30;
y = chirp(t,f0,t1,f1);
plot(y)

A_chirp = A1 + y' + 0.4*randn(size(A1)); % random noise
figure
wiggle(A_chirp)

%%
%Kalman

q = 0.005;
xf=zeros(1501,100);
kk=zeros(1501,100);
pp=zeros(1501,100);

for j=1:100
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(A_chirp(:,j),q); % kalman filter
end

wiggle(xf)

%SNR
res=xf-A_chirp;
snrKalman=10*log10(sum(sum(xf.^2))/sum(sum(res.^2))) % snr of kalman
%%
[nt nx] = size(A2);

dt1 = 0.001;
fs=1/dt1;
% N=l000;
f = [-nt/2:nt/2-1]*fs/1000;

dd2 = fftshift(abs(fft(A_chirp(:,:))));

for i = 1 : 1 : nt
f_d2(i,:) = max(dd2(i,:));
end

figure
subplot(223);
plot(f,normalize(f_d2) ,'LineWidth',1.5,'color','k')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

title(' Data')

%

dd = fftshift(abs(fft(xf(:,:))));

for i = 1 : 1 : nt
f_dd(i,:) = max(dd(i,:));
end

subplot(224);
plot(f,normalize(f_dd) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])
% ylim([0 220])

title('Kalman Filter')

subplot(121);
wiggle(A_chirp)
title('data + chirp +gaussian noise')


subplot(122);
wiggle(xf)
title('Kalman')
%%
f1 = 40; f2 = 400;dt=.002;
[out] = ormsbybandpassfilter(A_chirp,f1,f2,dt);
figure
wiggle(out)

%%
% Chirp
% Chirp + Gaussian noise 40%
t = 0 : 0.01 : 1.50;
t1 = 0.05;
f1 = 30;
y = chirp(t,f0,t1,f1);
plot(y)
noise = 0.8*rand(size(t));
A_C =  y + noise ; % random noise

[A_K,KK,P] = Kalman1D(A_C, 0.03);

subplot(411)
plot(y)
title("Chirp")
ylim([-2 2])

subplot(412)
plot(noise)
title("Noise")
ylim([-2 2])

subplot(413)
plot(A_C)
title("Chirp + Gaussian Noise")
ylim([-2 2])

subplot(414)
plot(A_K)
title("Kalman")
ylim([-2 2])









