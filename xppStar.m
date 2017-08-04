% peak-to-peak variation of the waveform defined as xpp parameter in (Albu,
% 2010) eq. (3) @ pp. 1112
function XppStar= xppStar(x, xDC)
Xmax=max(x);
Xmin=min(x);

XppStar=(Xmax-Xmin)/xDC;