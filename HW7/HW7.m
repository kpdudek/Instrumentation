% Get all that nasty shit outta here
close all; clear; clc;

% Insert that shit into my workspace
load('HW7_Problem1.mat')

% % Massage that shit
% ssVal = Voltage(find(t > .004));
% initVal = find(t<0);
% Voltage(initVal) = 1;

% Plot that shit
figure('Name','Raw Data')
plot(t,Voltage)
ylim([.8 2.6])
xlim([-1 5].*.001)
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')

% Display that useful shit
disp(stepinfo(Voltage,t))

% % Find those spikey bois
% [peakVals,locs] = findpeaks(-1*Voltage);

% % Linearize that decaying shit
% tLin = t(locs(1:3));
% linVals = -log(Voltage(locs(1:3))) ./ tLin;
% figure('Name','Linear')
% plot(tLin,linVals,'.','MarkerSize',12)

% coeff = polyfit(tLin,linVals,1);
% linfit = polyval(coeff,tLin);
% hold on; plot(tLin,linfit,'r','LineWidth',2)
% fprintf('Slope of the line: %5.3e\n',coeff(1))

K = 1;
A = 1;
wd = 14920;
% z = sqrt(1-(wd/wn)^2);

F = @(x,t) K*A - (K*A)*exp(-x(1)*x(2)*t)*((x(1)/sqrt(1-x(1)^2))*sin(wd*t)*cos(wd*t));
x0 = [1,1];

x = lsqcurvefit(F,x0,t,Voltage)








