function SphereDrag_Lab3()
close all
clear
clc

[load_coeff,Syx_load] = load_cell_calibration();
[press_coeff,Syx_press] = pressure_transducer_calibration();

load_coef_lab = [-6.3,53393.1]; % N , N/V
Syx_load_lab = 5*10^-6; % V
fprintf('\nLinear fit for load cell:\n')
fprintf('Mine | a0 : %5.5e a1 : %5.5e\n',load_coeff(2),load_coeff(1))
fprintf('Lab  | a0 : %5.5e a1 : %5.5e\n',load_coef_lab(1),load_coef_lab(2))

fprintf('\nStandard error of the fit for load cell:\n')
fprintf('Mine: %5.5e\n',Syx_load)
fprintf('Lab : %5.5e\n',Syx_load_lab)

press_coef_lab = [-75.13,2696.01]; % N/m^2 , N/m^2*V
Syx_press_lab = .052201;
fprintf('\nLinear fit for pressure transducer:\n')
fprintf('Mine | a0 : %5.5e a1 : %5.5e\n',press_coeff(2),press_coeff(1))
fprintf('Lab  | a0 : %5.5e a1 : %5.5e\n',press_coef_lab(1),press_coef_lab(2))

fprintf('\nStandard error of the fit for pressure transducer:\n')
fprintf('Mine: %5.5e\n',Syx_press)
fprintf('Lab : %5.5e\n',Syx_press_lab)


end

function [coeff,Syx] = load_cell_calibration()
% Measured lab values
loads = [0:.1:1,.9:-.1:0]; % Increase to 1kg and then decrease without double counting 1
force = loads.*9.8;

voltIncrease = [.1192,.1366,.1527,.1691,.1850,.2020,.2173,.2355,.2661,.2863,.3046]./1000;
voltDecrease = [.2862,.2680,.2497,.2330,.2131,.1948,.1766,.1583,.1398,.1200]./1000;
volts = [voltIncrease,voltDecrease];

% Plotting the raw data
figure('Name','Load Cell Calibration')
plot(force,volts,'.','MarkerSize',12)
xlabel('Force ($kg$)','Interpreter','latex')
ylabel('Voltage ($mV$)','Interpreter','latex')
hold on

% Linear fit on the data
coeff = polyfit(force,volts,1);
volt_fit = coeff(2) + coeff(1).*force;
plot(force,volt_fit)
% fprintf('\nLoad Cell Calibration\n')
% fprintf('a0 : %5.5e\na1 : %5.5e\n',coeff(2),coeff(1))

% Standard error of the fit
Syx = standard_error_fit(volts,volt_fit);

end

function [coeff,Syx] = pressure_transducer_calibration()
%Measured lab values
water_height = [.03,.94,2.4,3.4,4.6,5.9,7.3,8.8,7.3,5.85,4.5,3.4,2.4,.94,.04];
volt = [.002,.09,.239,.338,.456,.584,.762,.875,.726,.583,.455,.340,.239,.089,.002]./1000;
volt_std = [.017,.078,.077,.082,.081,.083,.072,.073,.071,.084,.084,.080,.079,.079,.017./1000];

% Plotting raw data
figure('Name','Pressure Transducer Calibration')
plot(water_height,volt,'.','MarkerSize',12)
xlabel('Dynamic Pressure ($\frac{N}{m^2}$)','Interpreter','latex')
ylabel('Voltage ($mV$)','Interpreter','latex')
hold on

% Linear fit on the raw data
coeff = polyfit(water_height,volt,1);
volt_fit = coeff(2) + coeff(1).*water_height;
plot(water_height,volt_fit)
% fprintf('\nPressure Transducer Calibration\n')
% fprintf('a0 : %5.5e\na1 : %5.5e\n',coeff(2),coeff(1))

% Standard error of the fit
Syx = standard_error_fit(volt,volt_fit);
end

function Syx = standard_error_fit(y_meas,y_fit)
nu = length(y_meas)-2;
Syx = sqrt(sum((y_fit-y_meas).^2)/nu);
%disp({'Syx: ',Syx})
end