% Earth Quake Kalman 
fileID = fopen('<local path>/IU.DAV.00.BHZ.M.2021-09-05T152300.019538.txt');
C = textscan(fileID, '%s %f');
Data = C{1,2};
Time = C{1,1};

%%
Q = 0.005;
[Y,KK,P]= Kalman1D(Data, Q );

figure
subplot(211)
plot(Data)
xlabel("Time (s)")
ylabel("Amplitude (m)")
ylim([-40000 40000])


subplot(212)
plot(Y)
xlabel("Time (s)")
ylabel("Amplitude (m)")
ylim([-40000 40000])

%%

dt = 0.025;
dx = 12.5;
[S1,k1,f1] = fk_spectra(Data',dt,dx,6);
imagesc(k1,f1,S1)
xlabel("Frequency (Hz)")
ylabel("Wave Number (1/km)")
% xlim([-0.01 0.01])

dt = 0.025;
dx = 12.5;
[S1,k1,f1] = fk_spectra(Y,dt,dx,6);
figure
imagesc(k1,f1,S1)
xlabel("Frequency (Hz)")
ylabel("Wave Number (1/km)")
% xlim([-0.01 0.01])
%%



[p,f] = pspectrum(Data);

f = f*1000;

plot(f,normalize(p,'range') ,"LineWidth",1.5,'color','k')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
xlim([0 250])



[p1,f1] = pspectrum(Y);
f1 = f1*1000

hold on

plot(f1,normalize(p1,'range') ,"LineWidth",1.5,'color','r');
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
xlim([0 250])
ylim([0 0.1])

legend("Data","Kalman")