% Time specifications:
Fs = 20000;                       % samples per second
 dt = 1/Fs;                       % seconds per sample
  StopTime = 10;                    % seconds
  t = (0:dt:StopTime-dt);             % seconds
%  t1 = (0:dt:.25);
%  t2 = (.25:dt:.50);
%  t3 = (.5:dt:.75);
%  t4 = (.75:dt:1);

%  %get a full-length example of each signal component
%  x1 = (10)*sin(2*pi*50*t);
%  x2 = (15)*sin(2*pi*150*t);
%  x3 = (12)*sin(2*pi*350*t);
% %  x4 = (10)*sin(2*pi*400*t);
% 
%   %construct a composite signal
%   x = zeros(size(t));
%   I = find((t >= t1(1)) & (t <= t1(end)));
%   x(I) = x1(I);
%   I = find((t >= t2(1)) & (t <= t2(end)));
%   x(I) = x2(I);
%   I = find((t >= t3(1)) & (t <= t3(end)));
%   x(I) = x3(I);
%   I = find((t >= t4(1)) & (t <= t4(end)));
%   x(I) = x4(I);
% 
%   NFFT = 2 ^ nextpow2(length(t));     % Next power of 2 from length of y
%    Y    = fft(x, NFFT);
%    f    = Fs / 2 * linspace(0, 1, NFFT/2 + 1);
%    figure;
%    plot(f(1:200), 2 * abs( Y( 1:200) ) );

   T = 0:.001:1;
   spectrogram(x,10,9);
   ylabel('Frequency');
    axis(get(gcf,'children'), [0, 1, 1, 100]);