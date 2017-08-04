function [uUniform, uGaussian,t]=generateSignal(xDC, fs, signalDuration, noiseAmplitude, fHar1, fHar2,fHar3,fHar4, AHar1, AHar2, AHar3, AHar4)
%function [uUniform, uGaussian]=generateSignal(xDC, fs, signalLength, noiseAmplitude, hormonic1, harmonic2, harmonic3)
%fs sampling frequency - e.g. 20kHz =20000
%signalLength=60 %seconds;
%xDC=230; the Dc component of the signal%Volts;
% noiseAmplitude=10;%Volt
%fs=20000; % samples per second
Ts=1/fs; %seconds per sample
%signalDuration=10 %seconds - when we stop the signal
signalLength=signalDuration*fs; %number of points 
t=(0:signalLength-1)*Ts; %  t = (0:Ts:signalDuration-Ts);             % seconds
% decimationFactor=1000;

Har1=AHar1*sin(2*pi*fHar1*t)';
Har2=AHar2*sin(2*pi*fHar2*t)';
Har3=AHar3*sin(2*pi*fHar3*t)';
Har4=AHar4*sin(2*pi*fHar4*t)';
uUniform=xDC+(-1+2.*rand(signalLength,1))*noiseAmplitude+Har1+Har2+Har3+Har4; 
% yUniformDecimated=decimate(uUniform,decimationFactor);
uGaussian=xDC+(-1+2*randn(signalLength,1))*noiseAmplitude+Har1+Har2+Har3+Har4;

% yGaussianDecimated=decimate(uGaussian,decimationFactor);
% % % figure;
% % % title 'original DC signal with noise and harmonic components'
% % % subplot 211
% % % plot(t,uUniform)
% % % % axis([0 fs*3/decimationFactor 0 xDC+10])
% % % grid on
% % % xlabel 'Time (sec)',ylabel 'Original uUniform'
% % % % subplot 222
% % % % plot(0:length(uUniform)/decimationFactor-1,yUniformDecimated)
% % % % axis([0 fs*3/decimationFactor 0 xDC+10])
% % % % grid on
% % % % xlabel 'Sample number',ylabel 'Decimated uUniform' 
% % % subplot 212
% % % plot(t,uGaussian)
% % % % axis([0 fs*3/decimationFactor 0 xDC+10])
% % % grid on
% % % xlabel 'Sample number',ylabel 'Original uGaussian'
% subplot 224
% plot(0:length(uUniform)/decimationFactor-1,yGaussianDecimated)
% axis([0 fs*3/decimationFactor 0 xDC+10])
% grid on
% xlabel 'Sample number',ylabel 'Decimated uGaussian'
% 

% subplot(2,2,1);
% plot(uUniform,0:(t-1));
% title('noisy signal with xDC and uniformly distributed noise of amplitude \noiseAmplitude');
% xlabel('time (sec)')