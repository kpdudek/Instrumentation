function Lab4()
close all
clc

%thermistor()
RTD()
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

function Syx = standard_error_fit(y_meas,y_fit)
nu = length(y_meas)-2;
Syx = sqrt(sum((y_fit-y_meas).^2)/nu);
%disp({'Syx: ',Syx})
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

function RTD()
% Read data from csv files
rtdData = read_data('rtd');
rtd10Data = read_data('rtd_10');
rtdDataVd = read_data('rtd_vd');
rtd10DataVd = read_data('rtd_vd_10');

% Averages and standard deviations
ave_rtdData = mean(rtdData);
ave_rtd10Data = mean(rtd10Data ./ (10^.5));

ave_rtdDataVd = mean(rtdDataVd);
ave_rtd10DataVd = mean(rtd10DataVd ./ (10^.5));

std_rtdData = std(rtdData);
std_rtd10Data = std(rtd10Data);
std_rtdDataVd = std(rtdDataVd);
std_rtd10DataVd = std(rtd10DataVd);


analyze_data(ave_rtdData,std_rtdData,'RTD with Bridge')
analyze_data(ave_rtd10Data,std_rtd10Data,'RTD with Bridge (10dB)')
analyze_data(ave_rtdDataVd,std_rtdDataVd,'RTD with Voltage Divider')
analyze_data(ave_rtd10DataVd,std_rtd10DataVd,'RTD with Voltage Divider (10dB)')

end

function analyze_data(volt,std,iden)
fprintf('\n---  %s | Data  ---\n',iden)
temp = [0,15,30,45,60];
tv_3 = 3.182;
tv_99 = 1.96;

% Perform linear regression
coeff = polyfit(temp,volt,1);
volt_fit = coeff(2) + coeff(1).*temp;
fprintf('a0 = %5.3e\na1 = %5.3e\n',coeff(2),coeff(1))

% Standard error of the fit
Syx = standard_error_fit(volt,volt_fit);
fprintf('Syx = %5.3e\n',Syx)

%%% Uncertainty
% Given
U_adc = .5 * (20/2^16);
UT_acc = 1;
UT_res = .05;
% Data specific
Up_fit = tv_3 * Syx;
Up_mean = tv_99.* (std./sqrt(100));

U_volt = sqrt(Up_fit^2 + Up_mean.^2 + U_adc^2);
U_temp = sqrt(UT_acc^2 + UT_res^2);

uncert = sqrt((U_volt/coeff(1)).^2 + U_temp^2);
fprintf('Up Fit = %5.3e\n',Up_fit)
fprintf('Up Mean = %5.3e\n',Up_mean)
fprintf('Total U Voltage = %5.3e\n',U_volt)
fprintf('Total U Temperature = %5.3e\n',U_temp)
fprintf('Uncertainty = %5.3f\n',uncert)

% Plot calibration
figure('Name',iden)
plot(temp,volt,'.','MarkerSize',12)
hold on
plot(temp,volt_fit,'r')
legend(iden,'Linear Fit')
xlabel('Temperature ($^{\circ}C$)','Interpreter','latex')
ylabel('Voltage (V)')
% xlim([min(temp)*.9 max(temp)*1.1])
% ylim([min(volt_fit)*.9 max(volt_fit)*1.1])

% Plot uncertainty
temp_func = (volt - coeff(2)) ./ coeff(1);
figure('Name',sprintf('%s, uncert',iden))
errorbar(volt,temp_func,uncert,'.','MarkerSize',12)
hold on
legend(iden)
ylabel('Temperature ($^{\circ}C$)','Interpreter','latex')
xlabel('Voltage (V)')
% xlim([min(volt)*.9 max(volt)*1.1])
% ylim([min(temp_func)*.9 max(temp_func)*1.1])

end


