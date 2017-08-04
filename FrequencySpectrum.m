close all
clc
load Measurement_32dc_36V_228Ohm;
signalNameLong='Voltage signal for dutyCy=0.32;Vout=36V;Rload=228\Omega';
Vin=(Measurement_32dc_36V_228Ohm.Y(3).Data);

signalName='Vin';
t=(Measurement_32dc_36V_228Ohm.X.Data);

signalDuration=60;
Fs=ceil(length(t)/max(t)/1000)*1000;
xDC=24;
desiredFs=50;
[y1] = decimate(Vin,10,5);
y=decimate(y1,10,5);
Ty=t(1:100:length(t));
figure (1);
plot(t,Vin,'-- ')
hold on
plot(Ty,y,'.-')
hold off
legend('Original','Resampled using decimate')
% ylim([-1.2 1.2])

sig=y(1:1000);
t1=Ty(1:1000);
% WHT.m creates a chirp signal and perform Wigner-Hough transform on it.

% t1=0:0.001:1; 
% N1 = length(t1);
% sig = chirp(t1,100,1,20);
sig = hilbert(sig);
% Signal in time domain
figure(2);
subplot(221);
plot(t1,real(sig));
xlabel('time (s)');
ylabel('Amplitude (V)');
title('signal in time')
sig_n = awgn(sig,-3);
subplot(222)
plot(t1, real(sig_n));
title('signal+WGN')
xlabel('time (s)');
ylabel('Amplitude (V)');
% Wigner-Ville Distribution
[tfr,t,f] = wv(sig);
% f = f * (N1-1)/2/N1;
t = t * 1/1000;
[F, T] = meshgrid(f, t);
subplot(223);
mesh(F, T, abs(tfr));
xlabel('frequency (Hz)');
ylabel('time (s)')
zlabel('Abs W-V distribution')
title('Wigner-Ville distribution');
subplot(224);
contour(f,t,abs(tfr));
xlabel('frequency (Hz)');
ylabel('time (s)');
% Hough Transform
[ht, rho, theta] = hough(tfr, f, t);
figure (3);
mesh(theta*180/pi, rho, abs(ht)); 
xlabel('theta'); ylabel('rho'); 
% figure; image(theta*(180/pi), rho, abs(ht));
% xlabel('theta'); ylabel('rho'); 
