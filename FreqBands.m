function []=FreqBands(u,Fs, meanU)
%OLD: function []=FreqBands(xDC,Fs, noiseAmplitude,samplingTime)
%% Frequency domain PQ analysis - derivation of meaningful parameters
%The frequency-domain representation of a signal carries information about
%the signal's magnitude and phase at each frequency. 
%This is why the output of the FFT computation is complex. A complex number, x,
%has a real x_r, and an imaginary part, x_i, such that
%x = x_r + i*x_i.
%The magnitude of x is computed as  sqrt{(x_r^2+x_i^2)}, and the phase of x
%is computed as  arctan{(x_i/x_r)}.
%in MATLAB, functions abs and angle, respectively, get the magnitude and phase of any complex number.

%********************************************
%Steps to follow in our analysis as they were discussed in the 20160818 Weekly Meeting:
% 1. FFT of each Tw of the signal u(t)
% 2. FFT of each Hamming Window*u(t) over each Tw
% 3. Compare the 2 above (e.g. plot on same figure)
% 4. Determine the frequency bandwidths
% 5. Compute energy density E(omega)
% 6. Compute an energy density for each frequency bandwidths observed from the FFTs above

%***************************************************************************
%Since the fft gives us the frequency representation of the signal,
% we want to look for the maximum (to get the frequencies with the maximum energy),
%and since the fft is a complex signal, we will take the absolute value first - look at the magnitude.
%The index will correspond to the position in the fft vector wich corresponds to 
%the normalized frequency with maximum energy.
%Last, if our signal has an offset (that means that there is a DC component),
%then we may want to get rid of that offset before taking the fft so that we do not get a max at the origin representing the DC component.



[maxValue,indexMax] = max(abs(fft(u-meanU))); %indexMax is the index where the max fft value can be found

%to get from indexMax to the actual frequency of interest,
% we need to know the length NFFT of the fft (same as the length of your signal),
NFFTMax = length(u);
% and the sampling frequency Fs. The signal frequency will then be:
frequencyMaxEnergy = indexMax * Fs / NFFTMax;

f = Fs/2*linspace(0,1,NFFTMax/2+1); % is the frequency vector associated with the (one-sided) power spectral density (PSD)
% For a given base frequency F (e.g. the classical AC grid is 50Hz), then,

fHar=f; % f/F Standardize to harmonic reference

HarOrder=[];
%The n_th_ harmonic will be at the bin nearest the locations
for n=1:10
    if n*fix(fHar)==n*fHar; %as noted unless we fix the sampling rate and times precisely 
    %                     %there won't necessarily be a frequency bin at, say, 50 Hz identically.
    HarOrder=[HarOrder,n];
    end
end
HarOrder

L=length(u);
NFFT = 2^nextpow2(L);    
U = fft(u,NFFT)/L;
f = Fs/2*linspace(0,Fs/2,NFFT/2+1);
figure;
plot(U)
title 'FFT of u'
figure;
plot(f)
title 'f - frequencies'

figure;
title 'Frequency bands'
subplot(2,1,1)
plot((1:50),U(1:50));
xlabel('time (milliseconds)')

NFFT = 2^nextpow2(length(U)); % Next power of 2 from length of y
Y = fft(U-meanU,NFFT)/length(U);
%plot(Y)
f = Fs/2*linspace(0,1,NFFT/2+1); %frequency

subplot(2,1,2)
% Plot single-sided amplitude spectrum
title('Single-Sided Amplitude Spectrum of y(t)');
xlabel('Frequency (Hz)');
bar(f,2*abs(Y(1:NFFT/2+1)))
ylabel('Amplitude');
% U(1) = 0; % remove the DC component for better visualization

%The output of the FFT is a complex vector containing information about the frequency content of the signal.
%The magnitude tells us the strength of the frequency components relative to other components.
%The phase tells us how all the frequency components align in time.

%Plot the magnitude and the phase components of the frequency spectrum of the signal.
%The magnitude is conveniently plotted in a logarithmic scale (dB). 
%The phase is unwrapped using the unwrap function so that we can see a continuous function of frequency.
magnitudeU = abs(U);        % Magnitude of the FFT
phaseU = unwrap(angle(U));  % Phase of the FFT
% % % figure;
% % % subplot(2,1,1);
% % % plot(magnitudeU);
% % % title('Magnitude response of the DC noisy signal @230V DC');
% % % xlabel('Frequency in kHz');
% % % ylabel('dB');
% % % 
% % % subplot(2,1,2);
% % % plot(phaseU);
% % % title('Phase response of the DC noisy signal @230V DC');
% % % xlabel('Frequency in kHz');
% % % ylabel('radius');

%***************************************************
% helperFrequencyAnalysisPlot1(F,magnitudeU,phaseU,NFFT);
% 
% use the code of the function
% helperFrequencyAnalysisPlot1(F,Ymag,Yangle,NFFT,ttlMag,ttlPhase) from the
% Matlab example: 

% Plot helper function for the FrequencyAnalysisExample
% Copyright 2012 The MathWorks, Inc.

figure
subplot(2,1,1)
plot(F(1:NFFT/2)/1e3,20*log10(magnitudeU(1:NFFT/2)));
if nargin > 4 && ~isempty(ttlMag)
  tstr = {'Magnitude response of the noisy DC signal',ttlMag};
else
  tstr = {'Magnitude response of the noisy DC signal'};
end
title(tstr)
xlabel('Frequency in kHz')
ylabel('dB')
grid on;
axis tight 
subplot(2,1,2)
plot(F(1:NFFT/2)/1e3,phaseU(1:NFFT/2));
if nargin > 5
  tstr = {'Phase response of the DC noisy signal',ttlPhase};
else  
  tstr = {'Phase response of the DC noisy signal'};
end
title(tstr)
xlabel('Frequency in kHz')
ylabel('radians')
grid on;
axis tight
%this plot is not clear what has inside because the figure doesn't match with the ones above
%*******************************************************

%****************************************************
%Take time intervals of 1 s each (Tw=1 sec => No of sampling elements in Tw is Fs) and make 60 analysis (60 seconds <=> 60 consecutive Tw):
% Define the block parameter.  Average in a Tw=1sec=Fs columns over 1 row wide window. In 1 sec we have Fs sampling points
blockSize = [1, Fs];
NFFTblock = length(blockSize);

% % % Block process the signal to replace every element in the Tw element wide block by [...](e.g. the FFT) of the sampling elements in the block.
% % % 1: define the FFT fucntion for use by blockproc().
% % FFTFilterFunction = @(theBlockStructure) fft(theBlockStructure.data(:),NFFTblock);
% % 
% % %1: Do the actual (block FFT down to smaller size array).
% % FFTDownSignal = blockproc(u, blockSize, FFTFilterFunction); %aka FFT for each Tw
% % FFTDownSignal(1)=0; % remove the DC component for better visualization
% % NFFTDownSignal=length(FFTDownSignal)
% % FDownSignal = ((0:1/NFFTDownSignal:1-1/NFFTDownSignal)*Fs).'; %sampling step
% % lengthFDownSignal=length(FDownSignal)
% % magnitudeDownSignal=abs(FFTDownSignal);
% % phaseDownSignal=unwrap(angle(FFTDownSignal));
% % % print result per blocks

figure
subplot(2,1,1)

plot(FDownSignal(1:NFFTDownSignal/2)/1e3,20*log10(magnitudeDownSignal(1:NFFTDownSignal/2)));
if nargin > 4 && ~isempty(ttlMag)
  tstr = {'Tw Magnitude response of the noisy DC signal',ttlMag};
else
  tstr = {'Tw Magnitude response of the noisy DC signal'};
end
title(tstr)
xlabel('Frequency in kHz')
ylabel('dB')
grid on;
axis tight 
subplot(2,1,2)
plot(FDownSignal(1:NFFTDownSignal/2)/1e3,phaseU(1:NFFTDownSignal/2));
if nargin > 5
  tstr = {'Tw Phase response of the DC noisy signal',ttlPhase};
else  
  tstr = {'Tw Phase response of the DC noisy signal'};
end
title(tstr)
xlabel('Frequency in kHz')
ylabel('radians')
grid on;
axis tight

%% %imi da niste aberatii acum de aceea e comentat

% 1: define the FFT function for use by blockproc().
% FFTFilterFunction = @(theBlockStructure) fft(theBlockStructure.data(:)); 
% % 2: define the FFTHamming function for use by blockproc().
% FFTHammingFilterFunction = @(theBlockStructure) fft(theBlockStructure.data(:).*hamming(Fs));
% %*********************************************
% %**************************************************************************************
% % 1: Do the actual (block FFT down to smaller size array).
% blockFFTDownSignal = blockproc(u, blockSize, FFTFilterFunction); %aka FFT for each Tw


% % 1: Do the actual (block FFT of the Hamming window*u signal on each Tw down to smaller size array).
% blockFFTHammingDownSignal = blockproc(u, blockSize, FFTHammingFilterFunction); %aka FFT for each Tw

% figure;
% subplot(2,2,1);
% plot(u);
% title('DC signal 230V DC with random noise of amplitude 2V');
% xlabel('Sampling @ Fs=20kHz for 1 min');
% ylabel('u(volt)');
% 
% subplot(2,2,2);
% U=fft(u)
% U(1)=0;
% plot(fft(u'));
% title('FFT applied to the entire DC noisy signal');
% xlabel('Freq (Hz)');
% ylabel('u(volt)');
% % figure(2);
% subplot(2,2,3);
% plot(blockFFTDownSignal'); %see if the mean corresponds to the DC component
% title('FFT of the u signal applied directly to each Tw=1 sec');
% xlabel('freq (Hz)');
% ylabel('u(volt)');
% 
% subplot(2,2,4);
% plot(blockFFTHammingDownSignal'); %see if the mean corresponds to the DC component
% title('FFT of the u signal with Hamming window applied to each Tw=1 sec');
% xlabel('freq (Hz)');
% ylabel('u(volt)');
