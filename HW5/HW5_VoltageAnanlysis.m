function HW5_VoltageAnanlysis()
clc;
load('FilterData_HW5.mat')

% Get the indices where the time is less than or equal to 0.2 seconds
iNoise = find(tvec<=.2);

% Isolate the noise times and voltages
tNoise = tvec(iNoise);
Noise = yn(iNoise);

% Isolate the signal times and voltages
tSignal = tvec(iNoise(end):end);
Signal = yn(iNoise(end):end);

% Find RMS of the noise and signal voltage values
rmsNoise = RMS(Noise);
rmsSignal = RMS(Signal);
fprintf('RMS Signal: %5.3f\n',rmsSignal)
fprintf('RMS Noise: %5.3f\n',rmsNoise)

% Compute the signal to noise ratio
SNR = 20*log10(rmsSignal/rmsNoise);
fprintf('SNR: %5.3f\n',SNR)

% Filter signal and analyze
FilSig = analyze_data(tvec,yn);

% Isolate the noise times and voltages
Filt_tNoise = tvec(iNoise);
Filt_Noise = FilSig(iNoise);

% Isolate the signal times and voltages
Filt_tSignal = tvec(iNoise(end):end);
Filt_Signal = FilSig(iNoise(end):end);

% Find RMS of the noise and signal voltage values
rmsFiltNoise = RMS(Filt_Noise);
rmsFiltSignal = RMS(Filt_Signal);
fprintf('RMS Filtered Signal: %5.3f\n',rmsFiltSignal)
fprintf('RMS Filtered Noise: %5.3f\n',rmsFiltNoise)

% Compute the signal to noise ratio
FiltSNR = 20*log10(rmsFiltSignal/rmsFiltNoise);
fprintf('SNR: %5.3f\n',FiltSNR)
end

function rms = RMS(val)

N = length(val);
sq_val = val .^2;
ave = sum(sq_val) ./ N;

rms = sqrt(ave);
end

function FilSig = analyze_data(tvec,yn)
% --- Visualize data in file ---
figure;plot(tvec,yn);grid
title('Data')
xlabel('Time (s)');ylabel('Amplitude')

figure;plot(tvec,yn);grid
title('Raw Data Zoomed In')
axis([0.195,0.205,-30,30])
xlabel('Time (s)');ylabel('Amplitude')

% Frequency analysis
n = length(tvec);
fs = n/max(tvec);
fvec = [1:n/2-1]*fs/n; % frequency vector
Yn = abs(fft(yn));
figure;plot(fvec*1e-3,Yn(2:n/2));grid
xlabel('Frequency (kHz)');ylabel('Magnitude')
title('Magnitude FFT of Raw Data')
axis([0 50,0,4.1e5])
% ---------------------

% --- Lowpass (Butterworth) filter response ---
wn = .1;                      % define normalized cutoff frequency
[bu,au] = butter(6,wn,'low'); % create lowpass filter characteristics
FilSig = filter(bu,au,yn);    % apply filter to data
% -----------------------------------

figure;
subplot(211);plot(tvec,yn);grid
axis([0.195,0.205,-30,30])
title('Raw Data')
xlabel('Time (s)');ylabel('Amplitude')
subplot(212);plot(tvec,FilSig);grid
title('Filtered Data')
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

figure;plot(tvec,FilSig);grid
title('Filtered Data Zoomed In')
axis([0.195,0.205,-30,30])
xlabel('Time (s)');ylabel('Amplitude')
end

