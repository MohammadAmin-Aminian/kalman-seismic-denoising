data = altreadsegy('<local path>/shot1.sgy');
[nt nx]=size(data)
shot = agc(data,100,1,0);
imagesc(shot);
colormap gray
colorbar
%%
q = 0.005;
xf=zeros(nt,nx);
kk=zeros(nt,nx);
pp=zeros(nt,nx);

for j=1:nx
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(shot(:,j),q); % kalman filter
end

figure
imagesc(xf);
colormap gray
colorbar
%SNR
res=xf-shot;
snrKalman=10*log10(sum(sum(xf.^2))/sum(sum(res.^2))) % snr of kalman
%%
dt1 = 0.001;
fs=1/dt1;
% N=l000;
f = [-nt/2:nt/2-1]*fs/1000;

dd2 = fftshift(abs(fft(shot(:,:))));

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

subplot(221);
imagesc(shot)
title('Data')
colormap gray
xlabel('Offset (m)');
ylabel('Two  Way Time (ms)');

subplot(222);
imagesc(xf)
title('Kalman')
colormap gray
xlabel('Offset (m)');
ylabel('Two  Way Time (ms)');


dd1 = fftshift(abs(fft(D1(:,:))));

for i = 1 : 1 : nt
f_dd1(i,:) = max(dd1(i,:));
end

plot(f,normalize(f_dd1) ,'LineWidth',1.5,'color','r')
%%
figure

plot(f,(f_d2) ,'LineWidth',1.5,'color','b')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])

hold on

plot(f,(f_dd) ,'LineWidth',1.5,'color','g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])

plot(f,(f_dd1) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 250])

legend("Data" , "Kalman" , "FX")
%%
  lf = 15;        % Operator length
  flow = 1;       % Min freq to process
  fhigh= 120;      % Max freq to process
  mu = 0.01;      % trade-off parameter

% Apply FX filter

  D1 = fx_decon(D,dtsec,lf,mu,flow,fhigh);
  figure
  imagesc(D1)
  title('FX')

colormap gray
xlabel('Offset (m)');
ylabel('Two  Way Time (ms)');
