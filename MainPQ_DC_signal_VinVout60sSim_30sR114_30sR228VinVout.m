%% Main: 
%PQ comparison analysis of experimental set-up and simulation for the
%same experiment in Simulink. We are looking only on voltage signals for
%this comparative analysis with an emphasis on the output signal from DCDC
%converter

close all;
clear all;
clc;

%load aquired experimental data from the lab test with a real boost
%converter connected on the input DC side through a three phase grid side
%inverter; The measurements were taken for the input and output voltage and current signals for the following cases:
%   (a) by changing the duty-cycle and keeping the load (otput resistance) constant
%   (b) by changing the load conneted to the DC bus through the boost
%   converter and varying again the duty-cylce

%% load experimental deta collected accoring to the explanation above
load Measurement_52dc_48V_228Ohm.mat
load Measurement_57dc_48V_114Ohm.mat
load VinVoutSim60sec.mat


tSim=VacIacVdcSim.time;
FsamplSim=ceil(length(tSim)/max(tSim)/1000)*1000; 
VinSim60sec=VacIacVdcSim.signals(1).values;
VoutSim60sec=VacIacVdcSim.signals(2).values;
%retreave the input and output voltages in separate variables
% VinSim60sec is the name of the variable for Vin siumlation composed for
% the signal
Vout228=(Measurement_52dc_48V_228Ohm.Y(4).Data)';%Vout228 is the output voltage for a load R=228Ohm, Vinref=24V; Voutref=48V, dc=0.52
Vin228=(Measurement_52dc_48V_228Ohm.Y(3).Data)';%Vin228 is the input voltage for a load R=228Ohm, Vinref=24V; Voutref=48V, dc=0.52
t228=(Measurement_52dc_48V_228Ohm.X.Data)'; %the time stamp of the signals Vout228/Vin228
Vout114=(Measurement_57dc_48V_114Ohm.Y(4).Data)';%Vout114 is the output voltage for a load R=114Ohm, Vinref=24V; Voutref=48V,dc=0.57
Vin114=(Measurement_57dc_48V_114Ohm.Y(3).Data)';%Vin114 is the input voltage for a load R=114Ohm, Vinref=24V; Voutref=48V,dc=0.57
t114=(Measurement_57dc_48V_114Ohm.X.Data)'; %the time stamp of the signals Vout114/Vin114
%% Keep only the first minute of the stored signals from the experiment and copy them in corresponding variables with '_60s' at the end
sD114=30; %first 30s from the signal with R=114 to be concatenated with
sD228=30; %31-60s from the signal with R=228;
signalDuration=sD114+sD228;%signal duration in seconds = sD1+sD2 signalDuration=60; 

Fsampl114=ceil(length(t114)/max(t114)/1000)*1000; %sampling frequency
Fsampl228=ceil(length(t228)/max(t228)/1000)*1000; %sampling frequency
Vin114_30s=Vin114(1:sD114*Fsampl114);
Vout114_30s=Vout114(1:sD114*Fsampl114);
Vin228_30s=Vin228(sD228*Fsampl228+1:(sD228+30)*Fsampl228);
Vout228_30s=Vout228(sD228*Fsampl228+1:(sD228+30)*Fsampl228);
Vin114_228_60s=[Vin114_30s;Vin228_30s];
t114_60s=t114(1:signalDuration*Fsampl114);
t228_60s=t228(1:signalDuration*Fsampl228);
Vout114_228_60s=[Vout114_30s;Vout228_30s];

%start numbering the figures to be stored for this small project
numFigures=1;
%Plot in a mirroing way the two pairs of Vin-Vout for the two load values
%R1=228Ohm and R2=114 Ohm - The plot will be limited to 1minute

% figure(1);
% subplot(2,2,1);
% plot(t114_60s,Vin114_60s);
% xlabel='time(s)';
% ylabel='Vin(V)';
% title ('Vin @R=114Ohm & dc=0.57');
% 
% subplot(2,2,2);
% plot(t114_60s,Vout114_60s);
% xlabel='time(s)';
% ylabel='Vout(V)';
% title ('Vout @R=114Ohm & dc=0.57');
% 
% subplot(2,2,3);
% plot(t228_60s,Vin228_60s);
% xlabel='time(s)';
% ylabel='Vin(V)';
% title ('Vin @R=228Ohm & dc=0.52');
% 
% subplot(2,2,4);
% plot(t228_60s,Vout228_60s);
% xlabel='time(s)';
% ylabel='Vout(V)';
% title ('Vout @R=228Ohm & dc=0.52');
% 
% numFigures=numFigures+1;

%% signal from real Lab set-up and Experiments - load the voltage and current signals for
% case  - set up of the dutyCycle=0.57, VoutRef=48 %V, resistiveLoad=114
% %Ohm, Vin=24V
Tw=10;%analysis window of 10 sec
%***************************************de aici back 4 ori


name='VinSim114_228_60s';
title(name, 'fontSize',10)
% signalDuration=60; %signal duration in seconds
signalName='VinSim_114_228Ohm_60sec';
% t=(Measurement_57dc_48V_114Ohm.X.Data)';
% Fsampl=ceil(length(t)/max(t)/1000)*1000;
xDC=22;
% xDCVinSim114_228=median(VinSim60sec);
% LB=min(VinSim60sec);
% UB=max(VinSim60sec);
% xDCVinSim114_228=(UB-LB)/2+LB;

[EnDominantFreqVout, dominantFreq, TimeDomParam114]=PQ_DC_signal_1(VinSim60sec(100:end),FsamplSim(100:end),xDC,tSim,signalDuration,signalName, Tw, numFigures);
numFigures=numFigures+4;
%*********************
name='VoutSim114_228_60s';
title(name, 'fontSize',10)
% signalDuration=60; %signal duration in seconds
signalName='VoutSim_114_228Ohm_60sec';
% t=(Measurement_57dc_48V_114Ohm.X.Data)';
% Fsampl=ceil(length(t)/max(t)/1000)*1000;
xDC=48;
% xDCVinSim114_228=median(VinSim60sec);
% LB=min(VinSim60sec);
% UB=max(VinSim60sec);
% xDCVinSim114_228=(UB-LB)/2+LB;

[EnDominantFreqVout, dominantFreq, TimeDomParam114]=PQ_DC_signal_1(VoutSim60sec,FsamplSim,xDC,tSim,signalDuration,signalName, Tw, numFigures);
numFigures=numFigures+4;
%*************************

name='Vout114_228_60s';
title(name, 'fontSize',10)
% signalDuration=60; %signal duration in seconds
signalName='Vin57_52dc_114_228Ohm_60sec';
% t=(Measurement_57dc_48V_114Ohm.X.Data)';
% Fsampl=ceil(length(t)/max(t)/1000)*1000;
xDC=22;
% xDCVin114_228=median(Vin114_228_60s);
% UB=max(Vin114_228_60s);
% LB=min(Vin114_228_60s);
% xDCVin114_228=(UB-LB)/2+LB;
[EnDominantFreqVout, dominantFreq, TimeDomParam114]=PQ_DC_signal_1(Vin114_228_60s,Fsampl114,xDC,t114_60s,signalDuration,signalName, Tw, numFigures);

numFigures=numFigures+4;
name='Vout114_228_60s';
title(name, 'fontSize',10)
% signalDuration=60; %signal duration in seconds
signalName='Vout57_52dc_114_228Ohm_60sec';
% t=(Measurement_57dc_48V_114Ohm.X.Data)';
% Fsampl=ceil(length(t)/max(t)/1000)*1000;
xDC=21;
% xDCVout114_228=median(Vout114_228_60s);
% UB=max(Vout114_228_60s);
% LB=min(Vout114_228_60s);
% xDCVout114_228=(UB-LB)/2+LB;
[EnDominantFreqVout, dominantFreq, TimeDomParam114]=PQ_DC_signal_1(Vout114_228_60s,Fsampl114,xDC,t114_60s,signalDuration,signalName, Tw, numFigures);






