data = altreadsegy('<local path>/shot20.sgy');

shot = agc(data,100,1,0);

figure

imagesc(shot)

colormap gray

[nt nx] = size(shot)


title("Data")
xlabel('Offset (m)');
ylabel('time (ms)');
%%

A4=fx_decon(shot,.001,8,0,0,200); % FX filter


figure

imagesc(A4)

colormap gray

title("FX")
xlabel('Offset (m)');
ylabel('time (ms)');
%%


q = 0.005;

xf=zeros(nt,nx);
kk=zeros(nt,nx);
pp=zeros(nt,nx);
for j=1:nx
    [xf(:,j),kk(:,j),pp(:,j)]=Kalman1D(shot(:,j),q); % kalman filter
end
figure
imagesc(xf)
colormap gray

title("Kalman")
xlabel('Offset (m)');
ylabel('time (ms)');
%%


[nt nx] = size(shot);

dt1 = 0.001;
fs=1/dt1;
% N=l000;
f = [-nt/2:nt/2-1]*fs/1000;

dd2 = fftshift(abs(fft(shot(:,:))));

for i = 1 : 1 : nt
f_d2(i,:) = max(dd2(i,:));
end

figure
plot(f,(f_d2) ,'LineWidth',1.5,'color','k')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 400])

hold on

dd = fftshift(abs(fft(A4(:,:))));

for i = 1 : 1 : nt
f_dd(i,:) = max(dd(i,:));
end

plot(f,(f_dd) ,'LineWidth',1.5,'color','r')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 400])

hold on

mm2 = fftshift(abs(fft(xf(:,:))));

for i = 1 : 1 : nt
myfil(i,:) = max(mm2(i,:));
end
plot(f,(myfil) ,'LineWidth',1.5,'color','g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 400])
%%



  [S,k,f] = fk_spectra(shot,0.001,12.5,6);
  
  figure
  
  imagesc(k,f,S);

  title("FK Data")
  
  ylim([0 150])


  [S,k,f] = fk_spectra(xf,0.001,12.5,6);
  
  figure
  
  imagesc(k,f,S);

  title("FK Kalman")
  ylim([0 150])

