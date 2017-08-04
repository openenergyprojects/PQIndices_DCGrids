%% Load my experimental data
%load each subset of aquired experimental data from the lab test with a real boost
%converter connected on the input DC side through a three phase grid side
%inverter; The measurements were taken for the input and output voltage and current signals for the following cases:
%   (a) by changing the duty-cylce and keeping the load (otput resistance)
%   constant: two files: "Measurement_32dc_36V_228Ohm" and
%   "Measurement_52dc_48V_228Ohm;"
%   (b) by changing the load conneted to the DC bus through the boost
%   converter and varying again the duty-cylce: we have another two files:
%   "Measurement_37dc_36V_114Ohm" and "Measurement_57dc_48V_114Ohm"

%% first set of data:
dataSet1.Name='Voltage and current signals for dutyCy=0.32;Vout=36V;Rload=228\Omega';
% title(dataSet1.Name, 'fontSize',10) %plot super title

dataSet1.t= Measurement_32dc_36V_228Ohm.X.Data;
dataSet1.Iin=Measurement_32dc_36V_228Ohm.Y(1).Data;
dataSet1.Vin=Measurement_32dc_36V_228Ohm.Y(3).Data;
dataSet1.Iout=Measurement_32dc_36V_228Ohm.Y(2).Data;
dataSet1.Vout=Measurement_32dc_36V_228Ohm.Y(4).Data;

%% second set of data:
dataSet2.Name='Voltage and current signals for dutyCy=0.52;Vout=48V;Rload=228\Omega';

dataSet2.t= Measurement_52dc_48V_228Ohm.X.Data;
dataSet2.Iin=Measurement_52dc_48V_228Ohm.Y(1).Data;
dataSet2.Vin=Measurement_52dc_48V_228Ohm.Y(3).Data;
dataSet2.Iout=Measurement_52dc_48V_228Ohm.Y(2).Data;
dataSet2.Vout=Measurement_52dc_48V_228Ohm.Y(4).Data;


%% third set of data:

dataSet3.Name='Voltage and current signals for dutyCy=0.37;Vout=36V;Rload=114\Omega';

dataSet3.t= Measurement_37dc_36V_114Ohm.X.Data;
dataSet3.Iin=Measurement_37dc_36V_114Ohm.Y(1).Data;
dataSet3.Vin=Measurement_37dc_36V_114Ohm.Y(3).Data;
dataSet3.Iout=Measurement_37dc_36V_114Ohm.Y(2).Data;
dataSet3.Vout=Measurement_37dc_36V_114Ohm.Y(4).Data;

%% fourth set of data:
load Measurement_57dc_48V_114Ohm;

dataSet4.Name='Voltage and current signals for dutyCy=0.57;Vout=48V;Rload=114\Omega';

dataSet4.t= Measurement_57dc_48V_114Ohm.X.Data;
dataSet4.Iin=Measurement_57dc_48V_114Ohm.Y(1).Data;
dataSet4.Vin=Measurement_57dc_48V_114Ohm.Y(3).Data;
dataSet4.Iout=Measurement_57dc_48V_114Ohm.Y(2).Data;
dataSet4.Vout=Measurement_57dc_48V_114Ohm.Y(4).Data;


%for each set we do the same analysis, however me may be interested in
%plotting on the same figure some of the data (e.g. to notice the effect of changing the load on the output voltage)

 t       = dataSet1.t ;
 nSeries = length( dataSet1, 2 ) ;
 % - Build figure.
 figure() ;  clf ;
 set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
 % - Compute #rows/cols, dimensions, and positions of lower-left corners.
 nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
 rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
 % - Build subplots axes and plot data.
 for dId = 1 : nSeries
    rowId = ceil( dId / nCol ) ;
    colId = dId - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    plot( t, data(:,dId), 'b' ) ;
    grid on ;
    xlabel( '\theta(t) [rad]' ) ;  ylabel( 'Anomaly [m]' ) ;
    title( sprintf( 'Time series %d', dId )) ;    
 end
 % - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'My Nice Title', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;