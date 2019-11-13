% Get all that nasty shit outta here
close all; clear; clc;

% Insert that shit into my workspace
load('HW7_Problem1.mat')

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

% Find those spikey bois
[peakVals,locs] = findpeaks(-1.*Voltage);
hold on
plot(t(locs(1:5)),-1.*peakVals(1:5),'ro')
error = .1 * 1;
error = ones(1,length(t)) .* error;
plot(t,2+error,'b--',t,2-error,'b--')
plot(.001363,2.1,'g.','MarkerSize',16)
legend('Signal','Peaks','Upper Bound','Lower Bound','Settling Time')

% Linearize that decaying shit
tLin = t(locs(1:5));
linVals = -log(Voltage(locs(1:5)));
figure('Name','Linear')
plot(tLin,linVals,'.','MarkerSize',12)
xlabel('$Time (s)$','Interpreter','latex')
ylabel('$-ln(y_{env})$','Interpreter','latex')

% Linear fit
coeff = polyfit(tLin,linVals,1);
linfit = polyval(coeff,tLin);
hold on; plot(tLin,linfit,'r','LineWidth',2)
fprintf('Slope of the line: %5.3e\n',coeff(1))


% Second order parameters
wd = pi/(4.21*10^-4);
wn = 7460.1;
z = sqrt(abs(1-(wd/wn)^2));
fprintf('Zeta = %5.10f\n',z)

% Error
dError = abs((Voltage-2)/(1-2)) .* 100;
figure()
plot(t,dError,'LineWidth',2)
xlim([0 max(t)])
ylim([0 105])
xlabel('Time (s)')
ylabel('Dynamic Error (%)')










