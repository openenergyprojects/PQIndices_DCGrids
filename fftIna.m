% peak-to-peak variation of the waveform defined as xpp parameter in (Albu,
% 2010) eq. (3) @ pp. 1112
function [EnDominantFreq,dominantFreq]= fftIna(u, Fs, signalDuration, k, signalName)

%% we calculate the following steps in order to determine the relevalce of energy spread into the signal
% 1. Energy for the DC component  +/-5Hz
% 2. Energy of the dominating frequencies other than the DC component (takea bandwidth of +/-5Hz)
% 3. evaluate the rate between 2 and 1: 2/1
% 4. total energy
% 5. evaluate the rate between 1 & 4: 1/4
% 5. Calculate an energy coref En=(1-sum2)/4
            % signalName=strcat('FullSignal',signalName); - am incercat si varinata
            % full, dar componeneta DC e prea mare fata de zgomotul creat de riple sau
            % switching control actions 
L=length(u);

uNoMean=u-mean(u);

% NFFT = 2^nextpow2(L)   


% we need to normalize the fft to the number of 
                %sample points indicated when calling FFT in order to keep the magnitude of the signal
U = fft(uNoMean)/L;                
% U = fft(u)/L;  %iese doar 0 ce e sub amplitudinea componentei DC   - iar peaks va da eroare, caci are un vector cu 0 elemente            
                

f = Fs*(0:L-1)/L;
fNiquist=f(1:floor(L/2));%we look on half of the mirrored signal obtained from FFT call
UNiquist=U(1:floor(L/2));

%% get the dominated frequencies;
absUNiquist=abs(UNiquist); % we look only into the half mirroring signal
%find the peaks in the FFT of the signal of whose aplitude is grated than
%e.g. 1 V or another dynamic threshold (e.g. to be larger than the std(x))
%calculation of dynamic threshold to get the peaks
plot(fNiquist,UNiquist)
 
[pcs, positionPcs]=findpeaks(absUNiquist,'MinPeakProminence',std(absUNiquist));
%get the frequencies that correspond to these peaks
dominantFreq=f(positionPcs);

%% calculate the energy of the signal between dominant frequency bands
Energy=[];
EnergyAllSpectrum=[];
EnBeforeW1=1/(2*pi)*sum(absUNiquist(1:positionPcs(1)-1));
EnergyAllSpectrum(1:positionPcs(1)-1)=EnBeforeW1;
EnAfterWn=1/(2*pi)*sum(absUNiquist(positionPcs(end)+1:end));
for i=1:length(dominantFreq)-1
    for j=i+1
        absUFreqBand=absUNiquist(positionPcs(i):positionPcs(j)); %
        Energy(i)=1/(2*pi)*sum(absUFreqBand.^2);
        absUFreqBand=[];
        EnergyAllSpectrum(positionPcs(i):positionPcs(j))=Energy(i);
    end
end
EnergyAllSpectrum(positionPcs(end)+1:length(fNiquist))=EnAfterWn;
EnDominantFreq=sum(Energy)/sum(EnergyAllSpectrum);%aici am modificat in loc de Energy/sum(EnergyAllSpectrum)


%%%%%%%%%
%plot the spectogram
 % - Build figure.
 figure(k) ;  clf ;
 set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
subplot(2,2,1);
plot(fNiquist,abs(UNiquist));
% title ('FFT of u-mean(u)', 'fontsize', 9, 'fontweight', 'normal');
title ('FFT of u', 'fontsize', 9, 'fontweight', 'normal');
xlabel('frequency (Hz)', 'FontSize',8);
ylabel('Signal Amplitude (Volt)','FontSize',8);

subplot(2,2,2);
% sizeZoom=min(length(fNiquist),max(dominantFrequency)*signalDuration);
sizeZoom=min(length(fNiquist),300*signalDuration);
plot(fNiquist(1:sizeZoom),abs(UNiquist(1:sizeZoom)));
title('Zoom on FFT analysis on half mirror', 'fontsize', 9, 'fontweight', 'normal');
xlabel('frequency (Hz)','FontSize',8);
ylabel('Signal Amplitude (Volt)','FontSize',8);
% figure(7);
%now call the spectrogram
% spectrogram(x, window, noverlap, Nfft, Fs);
%define FFT parameters
des_df_Hz = 25;  %desired frequency resolution for the display, Hz
Nfft = round(Fs / des_df_Hz);  %general rule for FFT resolution
Nfft = 2*Nfft;  %double the bins to account for spreading due to windowing
Nfft = 2*round(0.5*Nfft);  %make Nfft an even number
window = Nfft;  %make your window the same length as your FFT
noverlap = round(0.95);  %overlap a lot to make the plot pretty
% subplot(2,2,3)
% spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
% title ('spectrogram of the signal', 'fontsize', 9, 'fontweight', 'normal');
% xlabel('frequency (kHz)', 'FontSize',8);
% ylabel('time (seconds)','FontSize',8);
subplot(2,2,3:4);
spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
title ('spectrogram', 'fontsize', 9, 'fontweight', 'normal');
ylabel('frequency (kHz)', 'FontSize',8);
xlabel('time (seconds)','FontSize',8);
view(90,90);
% suptitle('PQ indices for Frequency Domain analysis') ;

% - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'PQ indices for Frequency Domain analysis', 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
  
  
%save figure as a png file 
printK=strcat('-f',num2str(k));
% print(printK,strcat(signalName,strcat(num2str(k),'PQidicators4FreqDomain')),'-dpng');
set(gcf, 'PaperPositionMode','auto')   %# WYSIWYG

print('-dpng', '-r0', strcat(signalName,strcat(num2str(k),'PQidicators4FreqDomain')));               %# at screen resolution