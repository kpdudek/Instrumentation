close all; clc; clear;

% load('Forced_12.mat');
% load('Forced_100.mat');
% load('FreeDecay_LowerRelease.mat')
load('FreeDecay_UpperRelease.mat')

% Frequency
fs = 1/250000;

% Square the data to make all values positive
filtData = data.^2;

% Stays constant
trigThresh = 4;

% Must be adjusted based on pulse height
echoThresh = .0025;

% Increase if plot is moving right
% Decrease if plot is moving left
% Motor Speed 12: 9426
% Motor Speed 100: 9418
% Free Decay Upper: 9431 - start at 200 pulses
% Free Decay Lower: 9435
numberSamplesPerPulse = 9431;

% Arrays to store the data and a time variable to get absolute time
time = [];
dist = [];
prev_time = 0;

% Find the first trigger
firstEchoTime = min(find(filtData>trigThresh));

% The loop checks if the (current trigger time + pulse time) is in the
% range of the data set
while firstEchoTime+numberSamplesPerPulse <= length(filtData)
% Isolate the current trigger time plus the number of samples in a pulse to work with 
pulse=filtData(firstEchoTime:firstEchoTime+numberSamplesPerPulse);

% Uncomment for debugging the samples per pulse
% plot(pulse)
% hold on

% Go through the indices after the trigger width until the value is above
% the echo voltage threshold
foundPulse = 0;
sampleIndex = 350;
while ~foundPulse
    if sampleIndex > numberSamplesPerPulse
        disp("FUCK NO PULSE FOUND")
        sampleIndex = 350;
        break
    end
    
    if(pulse(sampleIndex)>echoThresh)
        foundPulse=1;
    else
        sampleIndex=sampleIndex+1;
    end
end
travelTime = sampleIndex*fs;
soundSpeed=343;

% Get the travel distance two ways. Using the speed of sound and our
% calibration
travelDistance_SoS = travelTime*soundSpeed/2;
travelDistance_Cal = (sampleIndex*fs+5.71e-4)/5.83e-3;

% Make sure the distances arent radically different
if abs(travelDistance_Cal-travelDistance_SoS) > .1
    disp('FUCK')
end

% Store the absolute time, and the travel distance
pulseTime = numberSamplesPerPulse * fs;
time = [time, prev_time + pulseTime]; %prev_time+travelTime];
dist = [dist, travelDistance_Cal];


% plot(sampleIndex,pulse(sampleIndex),'or')
% xlim([0 1200])
% ylim([0 .04])
% hold off
%%%%%% UPDATE LOOP VARIABLES %%%%%%%
% Add the samples per pulse
firstEchoTime = firstEchoTime+numberSamplesPerPulse;
prev_time = prev_time+pulseTime;
end

% Plot distance vs time
figure()
plot(time,dist)
xlabel('Time (s)')
ylabel('Distance From Sensor (m)')
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Decay Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xVal = [2.075,2.565,3.056,3.546,3.999,4.451,4.942];
yVal = [.3593,.4225,.4657,.4966,.5178,.5377,.5528];
plot(xVal,yVal,'ro')

figure('Name','Exponential Fit')
ln_env = -1*log(yVal);
plot(xVal,ln_env,'b.','MarkerSize',12)
hold on 
lnCoeff = polyfit(xVal,ln_env,1);
linDist = polyval(lnCoeff,xVal);
plot(xVal,linDist,'r','LineWidth',2)

xlabel('Time (s)')
ylabel('$-ln(Distance)$','Interpreter','latex')

fprintf('Linear slope: %5.5f\n',lnCoeff(1))

Td = mean([2.565-2.075, 3.056-2.565, 3.546-3.056, 3.999-3.546, 4.451-3.999, 4.942-4.451]);
fprintf('Td: %5.5f\n',Td)

wd = 2*pi / Td;
fprintf('Wd: %5.5f\n',wd)







