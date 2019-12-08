close all; clc; clear;

% load('Forced_12.mat');
load('Forced_100.mat');
% load('FreeDecay_LowerRelease.mat')
% load('FreeDecay_UpperRelease.mat')

% Frequency
fs = 1/250000;

% Square the data to make all values positive
filtData = data.^2;

% Stays constant
trigThresh = 4;

% Must be adjusted based on pulse height
echoThresh = .012;

% Increase if plot is moving right
% Decrease if plot is moving left
% Motor Speed 12: 9426
% Motor Speed 100: 9418
% Free Decay Upper: 9431 - start at 200 pulses
% Free Decay Lower: 9435
numberSamplesPerPulse = 9418;

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
plot(pulse)
hold on

% Go through the indices after the trigger width until the value is above
% the echo voltage threshold
foundPulse = 0;
sampleIndex = 400;
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


plot(sampleIndex,pulse(sampleIndex),'or')
xlim([0 1200])
ylim([0 .04])
hold off
xlabel('Samples')
ylabel('Distance From Sensor (m)')
%%%%%% UPDATE LOOP VARIABLES %%%%%%%
% Add the samples per pulse
firstEchoTime = firstEchoTime+numberSamplesPerPulse;
prev_time = prev_time+pulseTime;
end

% Plot it so we can see what the last frame looks like. It should start
% exactly at the beginning of a trigger
% figure()
% plot(pulse)


% for i = 1:length(dist)
%     if dist(i) > .8
%         dist(i) = .8;
%     end
% end

% Plot distance vs time
figure()
plot(time,dist)
xlabel('Time (s)')
ylabel('Distance From Sensor (m)')














% for iData = 2:length(filtData)-1
%     if (filtData(iData)==-1) && ((filtData(iData-1)~=0) || (filtData(iData+1)~=0))
%         filtData(iData)= mean([filtData(iData-1),filtData(iData+1)]);
%     end
% end
% 
% % plot the filtered data vs index
% 
% plot(filtData)
% 
% % Loop through the data and get time between signal going high and then to
% % low. This loop counts the time for both the trigger pulse and the echo
% % pulse and stores them in the same vector
% PulseTime = [];
% PulseHeight = [];
% time = 0;
% timeEnd = 0;
% trig = 0;
% pulseB = 0;
% for iData = 1:length(filtData)
%     if (filtData(iData) > 0) && (trig == 0) && (filtData(iData-1) == -1)
%         time = iData * fs;
%         
%         pulseB = iData;
%         trig = 1;
%     end
%     
%    if (filtData(iData) == -1) && (trig == 1) %&& (filtData(iData+1) == 0)
%        timeEnd = iData * fs;
%        PulseTime(1,end+1) = time;
%        PulseTime(2,end+1) = timeEnd;
%        
%        pulVolts = filtData(pulseB:iData);
%        PulseHeight(end+1) = sum(pulVolts);
%        time = 0;
%        trig = 0;
%    end 
% end
% len = 1:length(PulseTime(1,:));
% 
% end
% 
% 
% 
% 
% 
% % 
% % times = struct('Time',[]);
% % dataCount = 1;
% % timeCount = 0;
% % trig = 0;
% % for iData = 1:length(filtData)
% %     if filtData(iData) ~= 0 && trig == 0
% %         trig = 1;
% %     elseif (trig == 1) && (filtData(iData) == 0)
% %         trig = 2;
% %     elseif trig == 1
% %         timeCount = timeCount + fs;
% %     elseif trig == 2
% %         timeCount = timeCount + fs;
% %     elseif (trig == 2) && (filtData(iData) ~= 0)
% %         times(dataCount).time = timeCount;
% %         dataCount = dataCount + 1;
% %         trig = 3;
% %     elseif (trig == 3) && (filtData(iData) == 0)
% %         timeCount = 0;
% %         trig = 0;
% %     end
% % end