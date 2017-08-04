% peak-to-peak variation of the waveform defined as xpp parameter in (Albu,
% 2010) eq. (3) @ pp. 1112
function Xi75= xi75(x, xDC)
% Xmedian=prctile(x,50);
x75Plus=prctile(x,75);
x75Minus=prctile(x,25); %100-75%
Xi75=(x75Plus)/(x75Minus);