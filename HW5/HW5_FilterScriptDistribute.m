clear
close all
load FilterData_HW5.mat

% --- Visualize data in file ---
figure;plot(tvec,yn);grid
title('Data')
xlabel('Time (s)');ylabel('Amplitude')

figure;plot(tvec,yn);grid
title('Data, zoomed in')
axis([0.195,0.205,-30,30])
xlabel('Time (s)');ylabel('Amplitude')

% Frequency analysis
n = length(tvec);
fs = n/max(tvec);
fvec = [1:n/2-1]*fs/n; % frequency vector
Yn = abs(fft(yn));
figure;plot(fvec*1e-3,Yn(2:n/2));grid
xlabel('Frequency (kHz)');ylabel('Magnitude')
title('Magnitude FFT of data: the dominant frequencies correspond to the spikes.')
axis([0 50,0,4.1e5])
% ---------------------

% --- Lowpass (Butterworth) filter response ---
wn = .1;                      % define normalized cutoff frequency
[bu,au] = butter(6,wn,'low'); % create lowpass filter characteristics
FilSig = filter(bu,au,yn);    % apply filter to data
% -----------------------------------

% --- Apply filter to data ---
% [bu,au] = butter(6,wn,'low'); % create lowpass filter characteristics
% FilSig = filter(bu,au,yn);    % apply filter to data

figure;
subplot(211);plot(tvec,yn);grid
axis([0.195,0.205,-30,30])
title('Raw data')
xlabel('Time (s)');ylabel('Amplitude')
subplot(212);plot(tvec,FilSig);grid
title('Filtered data')
xlabel('Time (s)');ylabel('Amplitude')
axis([0.195,0.205,-30,30])
% -----------------------------------

% Frequency analysis
n = length(tvec);
fs = n/max(tvec);
fvec = [1:n/2-1]*fs/n; % frequency vector
FYn = abs(fft(FilSig));
figure;plot(fvec*1e-3,FYn(2:n/2));grid
xlabel('Frequency (kHz)');ylabel('Magnitude')
title('Magnitude FFT of Filtered Data')
axis([0 50,0,4.1e5])
% ---------------------

