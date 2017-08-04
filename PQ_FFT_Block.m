function [EnDominantFreqFFT]=PQ_FFT_Block(u,Fs)


%Take time intervals of 1 s each (Tw=1 sec => No of sampling elements in Tw is Fs) and make e.g. signalDuration=60 analysis (60 seconds <=> 60 consecutive Tw):

%***********************************************************************************
% Define the block parameter.  Average in a Tw=1sec=Fs columns over 1 row
% wide window. In 1 sec we have Fs sampling points
blockSize = [Fs,1];


% Block process the signal to replace every element in the Tw element wide block by [...](e.g. the mean) of the sampling elements in the block.
% 1: define the averaging function for use by blockproc().
fftFilterFunction = @(theBlockStructure) fftIna(theBlockStructure.data(:), Fs); %fft

%**************************************************************************************
% 1: Do the actual averaging (block average down to smaller size array).
EnDominantFreqFFT = blockproc(u, blockSize, fftFilterFunction); %aka our UU in the minutes of the meeting documnet

% Let's check the output size.
% [rows, columns] = size(blockAveragedDownSignal)

% 
% NP=length(0:1/Fs:1); %number of sample points in a time window interval of 1 second
 
% %% Frequency domain analysis
% ui=[];
% for blockIndex=1:signalDuration
%     ui(blockIndex)=u((blockIndex-1)*Fs:blockIndex*Fs)
% L=length(u);
% 
% uNoMean=u-mean(u);
% % NFFT = 2^nextpow2(L)   
% U = fft(uNoMean)/L; % we need to normalize the fft to the number of 
%                 %sample points indicated when calling FFT in order to keep the magnitude of the signal
%                 
% 
% f = Fs*(0:L-1)/L;
% fNiquist=f(1:floor(L/2));%we look on half of the mirrored signal obtained from FFT call
% UNiquist=U(1:floor(L/2));
% 
% %% get the dominated frequencies;
% absUNiquist=abs(UNiquist); % we look only into the half mirroring signal
% %find the peaks in the FFT of the signal of whose aplitude is grated than
% %e.g. 1 V or another dynamic threshold (e.g. to be larger than the std(x))
% %calculation of dynamic threshold to get the peaks
%  
% [pcs, positionPcs]=findpeaks(absUNiquist,'MinPeakProminence',std(absUNiquist));
% %get the frequencies that correspond to these peaks
% dominantFreq=f(positionPcs);
% 
% %% calculate the energy of the signal between dominant frequency bands
% Energy=[];
% EnergyAllSpectrum=[];
% EnBeforeW1=1/(2*pi)*sum(absUNiquist(1:positionPcs(1)-1));
% EnergyAllSpectrum(1:positionPcs(1)-1)=EnBeforeW1;
% EnAfterWn=1/(2*pi)*sum(absUNiquist(positionPcs(end)+1:end));
% for i=1:length(dominantFreq)-1
%     for j=i+1
%         absUFreqBand=absUNiquist(positionPcs(i):positionPcs(j)); %
%         Energy(i)=1/(2*pi)*sum(absUFreqBand.^2);
%         absUFreqBand=[];
%         EnergyAllSpectrum(positionPcs(i):positionPcs(j))=Energy(i);
%     end
% end
% EnergyAllSpectrum(positionPcs(end)+1:length(fNiquist))=EnAfterWn;
% EnDominantFreq=Energy/sum(EnergyAllSpectrum);
% 
% %% Wigner-Ville analysis
% 
% 
% 
% %% Plots for Time Domain
% %plot all PQ indices for the Time domain in a single figure: there are 12
% %subplots in total, including the plot of the signal itself over time
% 
% % nSeries=12; %number of subplots in the PQ time domain analysis
% % 
% %  % - Build figure.
% %  figure(1) ;  clf ;
% %  set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
% %     'Position', [0.1,0.1,0.6,0.6] ) ;
% %  % - Compute #rows/cols, dimensions, and positions of lower-left corners.
% %  nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
% %  rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
% %  colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
% %  rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
% %  % - Build subplots axes and plot data.
% % if length(signalName)<5 
% % ylabel1=signalName(1:4);
% % else
% %  ylabel1='x(t)'  ;
% % end;
% % ylabel2='p.u';
% % if isempty(strfind(signalName,'V'))==0
% %    ylabel1=strcat(ylabel1,'(volt)');
% % elseif isemply(strfind(signalName,'I'))==0
% %    ylabel1=strcat(ylabel1,'(Amps)');
% % end;
% %  for dId = 1 : nSeries
% %     rowId = ceil( dId / nCol ) ;
% %     colId = dId - (rowId - 1) * nCol ;
% %     axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
% %     switch dId
% %         case 1
% %             plot(t(1:10000),u(1:10000));
% %             title(signalName(1:min(length(signalName),10)), 'fontsize', 7, 'fontweight', 'normal');
% %             
% %         case 2
% %             plot(blockAveragedDownSignal); %see if the mean corresponds to the DC component
% %             title('Mean of the signal @ Tw=1 sec', 'fontsize', 7, 'fontweight', 'normal');
% %             xlabel('time (sec)', 'FontSize',8);
% %             ylabel(ylabel1, 'FontSize',8);
% %         case 3
% %             plot(blockRMSDownSignal);
% %             title('RMS of the signal @ Tw=1 sec', 'fontsize', 7, 'fontweight', 'normal');
% %             xlabel('time (sec)', 'FontSize',8);
% %             ylabel(ylabel1, 'FontSize',8);
% %         case 4
% %             plot(blockMedianDownSignal);
% %             title('Median=x50% of the signal @ Tw=1 sec', 'fontsize', 7, 'fontweight', 'normal');
% %             xlabel('time (sec)', 'FontSize',8);
% %             ylabel(ylabel1, 'FontSize',8);
% %         case 5
% %             plot(blockXppDownSignal);
% %             title('Peak-to-peak variation from median', 'fontsize', 7, 'fontweight', 'normal');
% %         case 6
% %             plot(blockXppStarDownSignal);
% %             title('Peak-to-peak variation from xDC', 'fontsize',7, 'fontweight', 'normal');
% %         case 7
% %             plot(blockX75DownSignal);
% %             title('75 percentile variation from median', 'fontsize', 7, 'fontweight', 'normal');
% %         case 8
% %             plot(blockXEDownSignal);
% %             title('RMS variation from xDC', 'fontsize', 7, 'fontweight', 'normal');
% %         case 9
% %             plot(blockXiPPDownSignal);
% %             title('peak-to-peak displacement factor', 'fontsize', 7, 'fontweight', 'normal');
% %         case 10
% %             plot(blockXi75DownSignal);
% %             title('75 percentile displacement factor', 'fontsize', 7, 'fontweight', 'normal');
% %         case 11
% %             plot(blockXiDownSignal);
% %             title('xi -asymmetry factor', 'fontsize', 7, 'fontweight', 'normal');
% %         case 12
% %             plot(blockXiRMSDownSignal);
% %             title('xiRMS -RMS variation displacement', 'fontsize', 7, 'fontweight', 'normal');
% % 
% %             
% %         otherwise
% %             disp('Unknown dId')
% %     end
% %     grid on ;
% %     xlabel('time (seconds)','FontSize',8);
% %     if dId<=4
% %         ylabel(ylabel1,'FontSize',8);
% %     else
% %         ylabel(ylabel2,'FontSize',8);
% %     end
% % 
% %  end
% %  % - Build title axes and title.
% %  axes( 'Position', [0, 0.95, 1, 0.05] ) ;
% %  set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
% % %  text( 0.5, 0, 'PQ indices for Time Domain analysis', 'FontSize', 14', 'FontWeight', 'Bold', ...
% % %       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
% % suptitle('PQ indices for Time Domain analysis');
% % %save figure as .png file   
% % print('-f1',strcat(signalName,'PQidicators4TimeDomain'),'-dpng');  
% % 
% % 
% % 
% % %% plots Frequency Domain parameters/(PQ indices)
% % 
% % %check the max, the min and the mean values. If they are 
% % % size(U)
% % % size(f)
% % figure (2);
% % % plot(f(1:1000),abs(U(1:1000)));
% % subplot(2,2,1);
% % plot(fNiquist,abs(UNiquist));
% % title ('FFT of u-mean(u)', 'fontsize', 9, 'fontweight', 'normal');
% % xlabel('frequency (Hz)', 'FontSize',8);
% % ylabel('Signal Amplitude (Volt)','FontSize',8);
% % 
% % subplot(2,2,2);
% % plot(fNiquist(1:3100*signalDuration),abs(UNiquist(1:3100*signalDuration)));
% % title('Zoom on FFT analysis on half mirror', 'fontsize', 9, 'fontweight', 'normal');
% % xlabel('frequency (Hz)','FontSize',8);
% % ylabel('Signal Amplitude (Volt)','FontSize',8);
% % % figure(7);
% % %now call the spectrogram
% % % spectrogram(x, window, noverlap, Nfft, Fs);
% % %define FFT parameters
% % des_df_Hz = 25;  %desired frequency resolution for the display, Hz
% % Nfft = round(Fs / des_df_Hz);  %general rule for FFT resolution
% % Nfft = 2*Nfft;  %double the bins to account for spreading due to windowing
% % Nfft = 2*round(0.5*Nfft);  %make Nfft an even number
% % window = Nfft;  %make your window the same length as your FFT
% % noverlap = round(0.95);  %overlap a lot to make the plot pretty
% % % subplot(2,2,3)
% % % spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
% % % title ('spectrogram of the signal', 'fontsize', 9, 'fontweight', 'normal');
% % % xlabel('frequency (kHz)', 'FontSize',8);
% % % ylabel('time (seconds)','FontSize',8);
% % subplot(2,2,3:4);
% % spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
% % title ('spectrogram reversed', 'fontsize', 9, 'fontweight', 'normal');
% % ylabel('frequency (kHz)', 'FontSize',8);
% % xlabel('time (seconds)','FontSize',8);
% % view(90,90);
% % suptitle('PQ indices for Frequency Domain analysis') ;
% % %save figure as a png file 
% % print('-f2',strcat(signalName,'PQidicators4FreqDomain'),'-dpng');