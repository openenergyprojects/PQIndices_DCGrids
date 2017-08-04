function [EnDominantFreq,dominantFreq, TimeDomParam]=PQ_DC_signal(u,Fs,xDC,t,signalDuration,signalName, Tw)


%Take time intervals of 1 s each (Tw=1 sec => No of sampling elements in Tw is Fs) and make 60 analysis (60 seconds <=> 60 consecutive Tw):

%***********************************************************************************
% Define the block parameter.  Average in a Tw=1sec=Fs columns over 1 row
% wide window. In 1 sec we have Fs sampling points
blockSize = [Fs,1];
perfectDCsignal=xDC*ones(Fs*signalDuration,1);
% u=u-perfectDCsignal;

%% In time domain calculate the following parameters:
%UU=mean(U) on each Tw; determine the mean value (in our case we have selected it 230 at the generation of the signal)
%RMS=sqrt(sum(x.^2)/n); determine the pure RMS value for each Tw interval, note
%that n=FS
%xE = RMS variation over the DC component, defined as xE=rms(x-xDC)/xDC---> eq.(6) in (Albu, 2010) @ pp. 1113
%Xipp= peak to peak displacement factor
%Xi75 = 75th percentile peak to peak displacement factor=> %XiY, Y% percentile displacement factor take => Xi75
%XiRMS = RMS variation displacement factor
%XiRMSstar = combined RMS displacement factor

% Block process the signal to replace every element in the Tw element wide block by [...](e.g. the mean) of the sampling elements in the block.
% 1: define the averaging function for use by blockproc().
meanFilterFunction = @(theBlockStructure) mean(theBlockStructure.data(:)); %mean2() calculates the mean of all elements in a matrix - in our case is a row vector
% 2: define the RMS function for use by blockproc()
rmsFilterFunction= @(theBlockStructure) rms(theBlockStructure.data(:));
% 3: define the median function for use in blockproc(): xDC,m=x50%
medianFilterFunction=@(theBlockStructure) prctile(theBlockStructure.data(:),50);
% 4: define the xpp function for use in blockproc(): check xpp.m file,
% where xpp is defined according to (Albu, 2010) - eq. (3) from pp 1112
xppFilterFunction=@(theBlockStructure) xpp(theBlockStructure.data(:)); %check xpp as defined in (Albu, 2010) - eq. (3) from pp 1112
% 5: define the xppStar function for use in blockproc(): check xppStar.m file,
% where xppStar is defined according to (Albu, 2010) - eq. (4) from pp 1112
xppStarFilterFunction=@(theBlockStructure) xppStar(theBlockStructure.data(:),xDC); %check xpp as defined in (Albu, 2010) - eq. (3) from pp 1112
% 6: define the x75 function for use in blockproc(): where x75 is defined as the 75% displacement percentile according to (Albu, 2010) - eq. (2) from pp 1112
x75FilterFunction=@(theBlockStructure) x75(theBlockStructure.data(:),xDC); %check x75% as defined in (Albu, 2010) - eq. (2) from pp 1112
% 7: define the xE function for use in blockproc(): where xE is the RMS
% variation according to (Albu, 2010) - eq. (6) from pp 1113
xEFilterFunction=@(theBlockStructure) rms(theBlockStructure.data(:)-xDC)/xDC; %check x75% as defined in (Albu, 2010) - eq. (2) from pp 1112
% 8: define the xiPP function for use in blockproc(): where xPP is the peak to peak displacement factor, according to (Albu, 2010) - eq. (8) from pp 1113
xiPPFilterFunction=@(theBlockStructure) xiPP(theBlockStructure.data(:)); 
% 9: define the xi75 function for use in blockproc(): where xi75 is defined as the 75% displacement percentile according to (Albu, 2010) - eq. (2) from pp 1112
xi75FilterFunction=@(theBlockStructure) xi75(theBlockStructure.data(:),xDC); 
% 10: define the xi function in blockproc(): where xi captures the asymmetry in the signal shape, eq.(9) from (Albu, 2010) @pp. 1113
xiFilterFunction=@(theBlockStructure) xi(theBlockStructure.data(:),xDC); 
% 11: define the xiRMS function in blockproc(): where xiRMS is the RMS variation displacement of the u signal, eq.(10) from (Albu, 2010) @pp. 1113
xiRMSFilterFunction=@(theBlockStructure) xiRMS(theBlockStructure.data(:),xDC); 

%**************************************************************************************
% 1: Do the actual averaging (block average down to smaller size array).
blockAveragedDownSignal = blockproc(u, blockSize, meanFilterFunction); %aka our UU in the minutes of the meeting documnet
blockAveragedDownSignalPerf = blockproc(perfectDCsignal, blockSize, meanFilterFunction);
% 2: do the actual RMS (block RMS value down to smaller size array).
blockRMSDownSignal = blockproc(u, blockSize, rmsFilterFunction); %aka our RMS in the minutes of the meeting documnet
blockRMSDownSignalPerf = blockproc(perfectDCsignal, blockSize, rmsFilterFunction);
% 3: do the actual xDC,m (block Median value down to smaller size array): aka our x,DC, m=x50% in (Ablu, 2010)
blockMedianDownSignal = blockproc(u, blockSize, medianFilterFunction);
blockMedianDownSignalPerf = blockproc(perfectDCsignal, blockSize, medianFilterFunction);
% 4: do the actual xpp (block xpp = peak to peak variation value down to smaller size array): aka our xpp% in (Ablu, 2010)
blockXppDownSignal = blockproc(u, blockSize, xppFilterFunction);
blockXppDownSignalPerf = blockproc(perfectDCsignal, blockSize, xppFilterFunction);
% 5: do the actual xppStar (block xppStar = peak to peak variation when instead of median we use xDC value down to smaller size array): aka our xpp* in (Ablu, 2010)
blockXppStarDownSignal = blockproc(u, blockSize, xppStarFilterFunction);
blockXppStarDownSignalPerf = blockproc(perfectDCsignal, blockSize, xppStarFilterFunction);
% 6: do the actual x75 (block x75 = 75% percentile displacement from xDC
blockX75DownSignal = blockproc(u, blockSize, x75FilterFunction);
blockX75DownSignalPerf = blockproc(perfectDCsignal, blockSize, x75FilterFunction);
% 7: do the actual xE (block xE = RMS variation over the xDC
blockXEDownSignal = blockproc(u, blockSize, xEFilterFunction);
blockXEDownSignalPerf = blockproc(perfectDCsignal, blockSize, xEFilterFunction);
% 8: do the actual xiPP (block)
blockXiPPDownSignal = blockproc(u, blockSize, xiPPFilterFunction);
blockXiPPDownSignalPerf = blockproc(perfectDCsignal, blockSize, xiPPFilterFunction);
% 9: do the actual xi75 (block)
blockXi75DownSignal = blockproc(u, blockSize, xi75FilterFunction);
blockXi75DownSignalPerf = blockproc(perfectDCsignal, blockSize, xi75FilterFunction);
% 10: do the actual xi (block)
blockXiDownSignal = blockproc(u, blockSize, xiFilterFunction);
blockXiDownSignalPerf = blockproc(perfectDCsignal, blockSize, xiFilterFunction);
% 11: do the actual xiRMS (block)
blockXiRMSDownSignal = blockproc(u, blockSize, xiRMSFilterFunction);
blockXiRMSDownSignalPerf = blockproc(perfectDCsignal, blockSize, xiRMSFilterFunction);

TimeDomParam.Mean=[blockAveragedDownSignal,blockAveragedDownSignalPerf];
TimeDomParam.RMS=[blockRMSDownSignal, blockRMSDownSignalPerf];
TimeDomParam.Median=[blockMedianDownSignal, blockMedianDownSignalPerf];
TimeDomParam.Xpp=[blockXppDownSignal, blockXppDownSignalPerf];
TimeDomParam.XppStar=[blockXppStarDownSignal, blockXppStarDownSignalPerf];
TimeDomParam.X75=[blockX75DownSignal, blockX75DownSignalPerf];
TimeDomParam.XE=[blockXEDownSignal, blockXEDownSignalPerf];
TimeDomParam.XiPP=[blockXiPPDownSignal, blockXiPPDownSignalPerf];
TimeDomParam.Xi75=[blockXi75DownSignal, blockXi75DownSignalPerf];
TimeDomParam.Xi=[blockXiDownSignal, blockXiDownSignalPerf];
TimeDomParam.XiRMS=[blockXiRMSDownSignal, blockXiRMSDownSignalPerf];
% Let's check the output size.
% [rows, columns] = size(blockAveragedDownSignal)

% 
% NP=length(0:1/Fs:1); %number of sample points in a time window interval of 1 second
 
%% Frequency domain analysis
% ui=[];

for blockIndex=1:signalDuration/Tw-1
%     ui(blockIndex)=u((blockIndex-1)*Fs+1:blockIndex*Fs);
    [EnDominantFreq, dominantFreq]= fftIna(u((blockIndex-1)*Fs*Tw+1:blockIndex*Fs*Tw), Fs, signalDuration,blockIndex,signalName);
end

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
%% Wigner-Ville analysis



%% Plots for Time Domain
%plot all PQ indices for the Time domain in a single figure: there are 12
%subplots in total, including the plot of the signal itself over time

nSeries=12; %number of subplots in the PQ time domain analysis

 % - Build figure.
 figure(1) ;  clf ;
 
 set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
 % - Compute #rows/cols, dimensions, and positions of lower-left corners.
 nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
 rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
 % - Build subplots axes and plot data.
if length(signalName)<5 
ylabel1=signalName(1:4);
else
 ylabel1='x(t)'  ;
end;
ylabel2='p.u';
if isempty(strfind(signalName,'V'))==0
   ylabel1=strcat(ylabel1,'(V)');
elseif isemply(strfind(signalName,'I'))==0
   ylabel1=strcat(ylabel1,'(A)');
end;
 for dId = 1 : nSeries
    rowId = ceil( dId / nCol ) ;
    colId = dId - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    switch dId
        case 1
            plot(t(1:10000),u(1:10000));
            title(signalName(1:min(length(signalName),10)), 'fontsize', 9, 'fontweight', 'normal');
            
        case 2
            plot(blockAveragedDownSignal); %see if the mean corresponds to the DC component
            hold on
            plot(blockAveragedDownSignalPerf,'g'); %see if the mean corresponds to the DC component
            hold on
            plot(moving(blockAveragedDownSignal,Tw),'r'); 
            hold off
            title('Mean of the signal @ Tw=1s', 'fontsize', 9, 'fontweight', 'normal');
            xlabel('time (s)', 'FontSize',8);
            ylabel(ylabel1, 'FontSize',8);
            
        case 3
            plot(blockRMSDownSignal);
            hold on
            plot(blockRMSDownSignalPerf,'g');
            hold on
            plot(moving(blockRMSDownSignal,Tw),'r')
            hold off
            title('RMS of the signal @ Tw=1s', 'fontsize', 9, 'fontweight', 'normal');
            xlabel('time (s)', 'FontSize',8);
            ylabel(ylabel1, 'FontSize',8);
        case 4
            plot(blockMedianDownSignal);
            hold on
            plot(blockMedianDownSignalPerf,'g');
            hold on
            plot(moving(blockMedianDownSignal,Tw),'r')
            hold off
            title('Median=x50% of the signal @ Tw=1s', 'fontsize', 9, 'fontweight', 'normal');
            xlabel('time (s)', 'FontSize',8);
            ylabel(ylabel1, 'FontSize',8);
        case 5
            plot(blockXppDownSignal);
            hold on
            plot(blockXppDownSignalPerf,'g');
            hold on
            plot(moving(blockXppDownSignal,Tw),'r')
            hold off
            title('Peak-to-peak variation from median', 'fontsize', 9, 'fontweight', 'normal');
        case 6
            plot(blockXppStarDownSignal);
            hold on
            plot(blockXppStarDownSignalPerf,'g');
            hold on
            plot(moving(blockXppStarDownSignal,Tw),'r')
            hold off
            title('Peak-to-peak variation from xDC', 'fontsize',9, 'fontweight', 'normal');
        case 7
            plot(blockX75DownSignal);
            hold on
            plot(blockX75DownSignalPerf,'g');
            hold on
            plot(moving(blockX75DownSignal,Tw),'r')
            hold off
            title('75 percentile variation from median', 'fontsize', 9, 'fontweight', 'normal');
        case 8
            plot(blockXEDownSignal);            
            hold on
            plot(blockXEDownSignalPerf,'g');            
            hold on
            plot(moving(blockXEDownSignal,Tw),'r')
            hold off
            title('RMS variation from xDC', 'fontsize', 9, 'fontweight', 'normal');
        case 9
            plot(blockXiPPDownSignal);            
            hold on
            plot(blockXiPPDownSignalPerf,'g');            
            hold on
            plot(moving(blockXiPPDownSignal,Tw),'r')
            hold off
            title('peak-to-peak displacement factor', 'fontsize', 9, 'fontweight', 'normal');
        case 10
            plot(blockXi75DownSignal);            
            hold on
            plot(blockXi75DownSignalPerf,'g');            
            hold on
            plot(moving(blockXi75DownSignal,Tw),'r')
            hold off
            title('75 percentile displacement factor', 'fontsize', 9, 'fontweight', 'normal');
        case 11
            plot(blockXiDownSignal);            
            hold on
            plot(blockXiDownSignalPerf,'g');            
            hold on
            plot(moving(blockXiDownSignal,Tw),'r')
            hold off
            title('xi -asymmetry factor', 'fontsize', 9, 'fontweight', 'normal');
        case 12
            plot(blockXiRMSDownSignal);            
            hold on
            plot(blockXiRMSDownSignalPerf,'g');            
            hold on
            plot(moving(blockXiRMSDownSignal,Tw),'r')
            hold off
            title('xiRMS -RMS variation displacement', 'fontsize', 9, 'fontweight', 'normal');
            
        otherwise
            disp('Unknown dId')
    end
    grid on ;
    xlabel('time (s)','FontSize',8);
    if dId<=4
        ylabel(ylabel1,'FontSize',8);
    else
        ylabel(ylabel2,'FontSize',8);
    end

 end
 % - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
%  text( 0.5, 0, 'PQ indices for Time Domain analysis', 'FontSize', 14', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
% suptitle('PQ indices for Time Domain analysis');

text( 0.5, 0, 'PQ indices for Time Domain analysis', 'FontSize', 12',  ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

%save figure as .png file   
% print('-f1',strcat(signalName,'PQidicators4TimeDomain'),'-dpng');  
set(gcf, 'PaperPositionMode','auto')   %# WYSIWYG
print('-dpng','-r0', strcat(signalName,'PQidicators4TimeDomain'));               %# at screen resolution



%% plots Frequency Domain parameters/(PQ indices)

%All the below code forms the fftIna.m function

%check the max, the min and the mean values. If they are 
% size(U)
% size(f)
% figure (2);
% % plot(f(1:1000),abs(U(1:1000)));
% subplot(2,2,1);
% plot(fNiquist,abs(UNiquist));
% title ('FFT of u-mean(u)', 'fontsize', 9, 'fontweight', 'normal');
% xlabel('frequency (Hz)', 'FontSize',8);
% ylabel('Signal Amplitude (Volt)','FontSize',8);
% 
% subplot(2,2,2);
% plot(fNiquist(1:3100*signalDuration),abs(UNiquist(1:3100*signalDuration)));
% title('Zoom on FFT analysis on half mirror', 'fontsize', 9, 'fontweight', 'normal');
% xlabel('frequency (Hz)','FontSize',8);
% ylabel('Signal Amplitude (Volt)','FontSize',8);
% % figure(7);
% %now call the spectrogram
% % spectrogram(x, window, noverlap, Nfft, Fs);
% %define FFT parameters
% des_df_Hz = 25;  %desired frequency resolution for the display, Hz
% Nfft = round(Fs / des_df_Hz);  %general rule for FFT resolution
% Nfft = 2*Nfft;  %double the bins to account for spreading due to windowing
% Nfft = 2*round(0.5*Nfft);  %make Nfft an even number
% window = Nfft;  %make your window the same length as your FFT
% noverlap = round(0.95);  %overlap a lot to make the plot pretty
% % subplot(2,2,3)
% % spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
% % title ('spectrogram of the signal', 'fontsize', 9, 'fontweight', 'normal');
% % xlabel('frequency (kHz)', 'FontSize',8);
% % ylabel('time (seconds)','FontSize',8);
% subplot(2,2,3:4);
% spectrogram(u, window, noverlap, Nfft, Fs,'yaxis');
% title ('spectrogram reversed', 'fontsize', 9, 'fontweight', 'normal');
% ylabel('frequency (kHz)', 'FontSize',8);
% xlabel('time (seconds)','FontSize',8);
% view(90,90);
% suptitle('PQ indices for Frequency Domain analysis') ;
% %save figure as a png file 
% print('-f2',strcat(signalName,'PQidicators4FreqDomain'),'-dpng');