function [PulseTime,PulseHeight,len] = analyzeSensorData(data,iden)
figure('Name',iden)

fs = 1/250000;
% time = 1:fs:length(data);
thresh = .01;

% Square the data to make all values positive and then set the noise to
% zero
filtData = data.^2;
filtData(find(filtData < thresh)) = -1;

echoThresh=4;
numberSamplesPerPulse = 9417;%1738000-1729000;
firstEchoTime = min(find(filtData>echoThresh));
while 1
pulse=filtData(firstEchoTime:firstEchoTime+numberSamplesPerPulse);
plot(pulse)
foundPulse = 0;
sampleIndex = 200;
while ~foundPulse
    if(pulse(sampleIndex)>-0.5)
        foundPulse=1;
    else
        sampleIndex=sampleIndex+1;
    end
end
travelTime = sampleIndex*fs;
soundSpeed=343;
travelDistance=travelTime/soundSpeed;
travelDistance=(sampleIndex*fs+5.71e-4)/5.83e-3
plot(pulse)
firstEchoTime=firstEchoTime+numberSamplesPerPulse;
end
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