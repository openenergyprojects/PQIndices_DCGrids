% peak-to-peak variation of the waveform defined as xpp parameter in (Albu,
% 2010) eq. (3) @ pp. 1112
function Xpp= xpp(x)
Xmax=max(x);
Xmin=min(x);
Xmedian=prctile(x,50);
Xpp=(Xmax-Xmin)/Xmedian;