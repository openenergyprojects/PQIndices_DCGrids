% peak-to-peak variation of the waveform defined as xpp parameter in (Albu,
% 2010) eq. (3) @ pp. 1112
function Xi= xi(x, xDC)
Xmedian=prctile(x,50);
xPlus=x(x>=Xmedian);
N=length(xPlus);
xMinus=x(x<=Xmedian);
M=length(xMinus);
Xi=sqrt(sum((xPlus/xDC).^2)/N)/sqrt(sum((xMinus/xDC).^2)/M);