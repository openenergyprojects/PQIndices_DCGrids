%Main
close all;
clear all;
clc;

%load aquired experimental data from the lab test with a real boost
%converter connected on the input DC side through a three phase grid side
%inverter; The measurements were taken for the input and output voltage and current signals for the following cases:
%   (a) by changing the duty-cylce and keeping the load (otput resistance) constant
%   (b) by changing the load conneted to the DC bus through the boost
%   converter and varying again the duty-cylce
load VacIacVdcOutWKS21s.mat;
load Measurement_52dc_48V_228Ohm;
Vout228=(Measurement_52dc_48V_228Ohm.Y(4).Data)';


%% signal from simulated Lab set-up and Experiments - load the voltage and current signals for
%1st case  - set up of the dutyCycle=0.48, VoutRef=48 %V, resistiveLoad=228 %Ohm
Tw=10;%analysis window of 10 sec

name='Voltage sim_sign for dutyCy=0.48;Vout=48V;Rload=114\Omega';

title(name, 'fontSize',10)
signalDuration=60; %signal duration in seconds
% Vin=(Measurement_52dc_48V_228Ohm.Y(3).Data)';
% signalName='Vin52dc_48V_228Ohm_10sec';
% t=(Measurement_52dc_48V_228Ohm.X.Data)';
% Fsampl=ceil(length(t)/max(t)/1000)*1000;
% xDC=24;
% u=Vin(1:signalDuration*Fsampl);
% [EnDominantFreqVin]=PQ_DC_signal(u,Fsampl,xDC,t,signalDuration,signalName, Tw);
Vout228VectorSim=(VacIacVdc2.signals(2).values(:,1))';
Vout228_sim=(repmat(Vout228VectorSim(150:end),1,3))';
xDC=48;
% t=(VacIacVdc2.time)',1,10);
% Fsampl=10000;
t=(Measurement_52dc_48V_228Ohm.X.Data)';
Fsampl=ceil(length(t)/max(t)/1000)*1000;
% Fsampl=10000;
signalName='Vout_sim_60s';
[EnDominantFreqVout,dominantFreq, TimeDomParam228]=PQ_DC_signal(Vout228_sim(1:signalDuration*Fsampl),Fsampl,xDC,t(1:signalDuration*Fsampl),signalDuration,signalName, Tw);

%% commented 20170228 *********************************
% % % % % % name='Voltage signal for dutyCy=0.57;Vout=48V;Rload=114\Omega';
% % % % % % 
% % % % % % title(name, 'fontSize',10)
% % % % % % signalDuration=60; %signal duration
% % % % % % % Vin=(Measurement_57dc_48V_114Ohm.Y(3).Data)';
% % % % % % signalName='Vin57dc_48V_114Ohm_10sec';
% % % % % % t=(Measurement_57dc_48V_114Ohm.X.Data)';
% % % % % % Fsampl=ceil(length(t)/max(t)/1000)*1000;
% % % % % % xDC=24;
%%**********************************************

% u=Vin(1:signalDuration*Fsampl);
% [EnDominantFreqVin]=PQ_DC_signal(u,Fsampl,xDC,t,signalDuration,signalName, Tw);

% % % % % Vout114=(Measurement_57dc_48V_114Ohm.Y(4).Data)';
% % % % % xDC=48;
% % % % % signalName='Vout57dc_48V_114Ohm_10sec';
% % % % % [EnDominantFreqVout, dominantFreq, TimeDomParam114]=PQ_DC_signal(Vout114(1:signalDuration*Fsampl),Fsampl,xDC,t,signalDuration,signalName, Tw);

% Iin=(Measurement_32dc_36V_228Ohm.Y(1).Data)';
% [EnDominantFreqVin]=PQ_DC_signal(Iin(1:signalDuration*Fsampl),Fsampl,xDC,t,signalDuration,signalName, Tw);
% 
% Iout=(Measurement_32dc_36V_228Ohm.Y(2).Data)';
% 
% Vin=(Measurement_37dc_36V_114Ohm.Y(3).Data)';
% xDC=24;
% signalName='Vin_37dc_36V_114Ohm_10sec';
% t=Measurement_37dc_36V_114Ohm.X.Data;
% [EnDominantFreqVin]=PQ_DC_signal(Vin(1:signalDuration*Fsampl),Fsampl,xDC,t,signalDuration,signalName, Tw);
% 
% Vout=(Measurement_37dc_36V_114Ohm.Y(4).Data)';
% signalName='Vout_37dc_36V_114Ohm_10sec';
% xDC=36;
% [EnDominantFreqVout]=PQ_DC_signal(Vout(1:signalDuration*Fsampl),Fsampl,xDC,t,signalDuration,signalName, Tw);




