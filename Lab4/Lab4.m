function Lab4()
close all
clc

thermistor()

end

function out = read_data(suffix)
out = [];
for iTemp = 1:5
    temp = {'0','15','30','45','60'};
    str = sprintf('%s%s.csv',temp{iTemp},suffix);
    data = csvread(str,23);
    data = data(:,2);
    out = [out,data];
end
end

function Ro = Volt_to_Resistance(Vo,Vi,R)
Ro = (R./(Vo./Vi)) - R;
end

function thermistor()
% Read data from csv files
data = read_data('cel');

% Get averages and standard deviations
avgData = mean(data);
stdData = std(data);

% Lab specific values
temp = [0,15,30,45,60]; % Water bath temps
Vapplied = 1.5; %V - voltage divider applied voltage
R = 980; %ohms - Voltage divider constant resistance

% Plot average voltage vs temp
figure('Name','Thermistor Volt')
plot(temp,avgData,'.','MarkerSize',12)
xlabel('Temperature ($^{\circ}C)$','Interpreter','latex')
ylabel('Voltage ($V$)','Interpreter','Latex')

% Convert average voltages to resistances
Rtherm = Volt_to_Resistance(avgData,Vapplied,R);
Rmeas = [7.357,3.511,1.818,1.006,.561].*1000; %ohms

% Plot calculated and measured resistance vs temp
figure('Name','Thermistor Resistance')
plot(temp,Rtherm,'.','MarkerSize',12)
hold on
plot(temp,Rmeas,'.','MarkerSize',12)
xlabel('Temperature ($^{\circ}C)$','Interpreter','latex')
ylabel('Voltage ($V$)','Interpreter','Latex')
legend('Voltage Divider','Measured')

% Getting data into linear form
lnR = log(Rtherm);
invT = 1./(temp+273); % convert temp to degrees K
figure('Name','Ln resistance vs 1/T')
plot(lnR,invT,'.','MarkerSize',12)
xlabel('$ln(R)$','Interpreter','latex')
ylabel('$\frac{1}{T}(^{\circ}C)$','Interpreter','latex')

end

function 


