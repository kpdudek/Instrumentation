function SphereDrag_Lab3()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SphereDrag_Lab3()
%
% Coefficients are stored in the vector as follows
%       coeffs = [ a0 , a1 ]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc

% Perform linear regression on the calibration data sets
[load_coeff,Syx_load] = load_cell_calibration();
[press_coeff,Syx_press] = pressure_transducer_calibration();

% Compare the lab displayed calibration parameters against mine
load_coef_lab = [-6.3,53393.1]; % N , N/V
Syx_load_lab = 5*10^-6; % V
fprintf('\nLinear fit for load cell:\n')
fprintf('Mine | a0 : %5.5e a1 : %5.5e\n',load_coeff(2),load_coeff(1))
fprintf('Lab  | a0 : %5.5e a1 : %5.5e\n',load_coef_lab(1),load_coef_lab(2))

fprintf('\nStandard error of the fit for load cell:\n')
fprintf('Mine: %5.5e\n',Syx_load)
fprintf('Lab : %5.5e\n',Syx_load_lab)

fprintf('-----------------------------------------------------------------------')

press_coef_lab = [-75.13,2696.01]; % N/m^2 , N/m^2*V
Syx_press_lab = .052201;
fprintf('\nLinear fit for pressure transducer:\n')
fprintf('Mine | a0 : %5.5e a1 : %5.5e\n',press_coeff(2),press_coeff(1))
fprintf('Lab  | a0 : %5.5e a1 : %5.5e\n',press_coef_lab(1),press_coef_lab(2))

fprintf('\nStandard error of the fit for pressure transducer:\n')
fprintf('Mine: %5.5e\n',Syx_press)
fprintf('Lab : %5.5e\n',Syx_press_lab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sphere drag analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spheres = {'Smooth','Low Rough','High Rough'};
sphere_diam = [.1014476,.1302,.10304]; %m
rho_water = 1.193; %997.59; % kg / m^3 for water
rho_air = 1.193; % kg / m^3
mu = 15.32 * 10^-6; % m^2 / s for air

manometer_smooth = [1.62,2.4,3.4,3.95,4.55,5.2,5.85,6.5,7.35,8.0,8.75,8.0,...
    7.25,6.5,5.8,5.15,4.5,3.9,3.4,2.9,2.4]';
volt_smooth = [.1341,.1433,.1540,.16,.1660,.1717,.1777,.1830,.1870,.1830,.1550,...
    .1825,.1872,.1837,.1780,.1720,.1661,.16,.1541,.1489,.1438]'./1000;
manometer_low = [1.62,2.4,3.4,3.95,4.55,5.2,5.9,6.6,7.3,8.05,8.7,8.0,...
    7.25,6.5,5.8,5.15,4.5,4.0,3.4,2.9,2.4]';
volt_low = [.1323,.1423,.1548,.16,.1650,.1695,.1733,.1730,.1717,.1680,...
    .1407,.1670,.1720,.1745,.1740,.1701,.1646,.1595,.1537,.1484,.1433]'./1000;
manometer_high = [1.62,2.4,3.4,4.0,4.2,5.2,5.9,6.6,7.3,8.0,8.7,8.0,7.3,...
    6.55,5.80,5.20,4.5,4.0,3.4,2.9,2.4]';
volt_high = [.1288,.1356,.1447,.1499,.1553,.161,.1677,.175,.1822,.19,...
    .1975,.19,.1822,.175,.1683,.162,.1557,.15,.1453,.1406,.1367]'./1000;

% Smooth sphere
[Cd_smooth,Re_smooth,D_smooth,vel_smooth,press_smooth] = sphere_drag(rho_water,rho_air,mu,spheres{1},sphere_diam(1),press_coeff,load_coeff,volt_smooth,manometer_smooth);

% Low rough sphere
[Cd_low,Re_low,D_low,vel_low,press_low] = sphere_drag(rho_water,rho_air,mu,spheres{2},sphere_diam(2),press_coeff,load_coeff,volt_low,manometer_low);

% High rough sphere
[Cd_high,Re_high,D_high,vel_high,press_high] = sphere_drag(rho_water,rho_air,mu,spheres{3},sphere_diam(3),press_coeff,load_coeff,volt_high,manometer_high);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncertainty Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ucd_smooth = cd_uncertainty(rho_air,vel_smooth,sphere_diam(1),D_smooth,press_smooth);
Ucd_low = cd_uncertainty(rho_air,vel_low,sphere_diam(2),D_low,press_low);
Ucd_high = cd_uncertainty(rho_air,vel_high,sphere_diam(3),D_high,press_high);


% Plotting the Coefficient of drag vs reynolds number for increasing and
% Decreasing tests
figure('Name','Re vs Cd Experimental')
hold on
errorbar(Re_smooth(1:11),Cd_smooth(1:11),Ucd_smooth(1:11),'ob');errorbar(Re_smooth(12:end),Cd_smooth(12:end),Ucd_smooth(12:end),'*k')
errorbar(Re_low(1:11),Cd_low(1:11),Ucd_low(1:11),'og');errorbar(Re_low(12:end),Cd_low(12:end),Ucd_low(12:end),'*r')
errorbar(Re_high(1:11),Cd_high(1:11),Ucd_high(1:11),'oc');errorbar(Re_high(12:end),Cd_high(12:end),Ucd_high(12:end),'*m')
legend('Increasing Smooth','Decreasing Smooth','Increasing Low','Decreasing Low','Increasing High','Decreasing High','Location','northeast')
xlabel('$Re$','Interpreter','latex')
ylabel('$C_d$','Interpreter','latex')
xlim([.9*min(Re_smooth) 1.1*max(Re_low)])

% Plotting the coefficient of drag vs reynolds number for increasing
% Plotting the computed theoretical value
[Cd_Theory] = compute_theoretical_Cd(Re_smooth(1:11),Re_low(1:11),Re_high(1:11));

figure('Name','Re vs Cd Theory');

errorbar(Re_smooth(1:11),Cd_smooth(1:11),Ucd_smooth(1:11),'ob')
hold on
errorbar(Re_low(1:11),Cd_low(1:11),Ucd_low(1:11),'og')
errorbar(Re_high(1:11),Cd_high(1:11),Ucd_high(1:11),'oc')

plot(Re_smooth(1:11),Cd_Theory(1,:),'r','LineWidth',2)
plot(Re_low(1:11),Cd_Theory(2,:),'y','LineWidth',2)
plot(Re_high(1:11),Cd_Theory(3,:),'k','LineWidth',2)

set(gca,'XScale','log', 'YScale','log')
xlim([.9*min(Re_smooth) 1.1*max(Re_low)])
grid on
legend('Increasing Smooth','Increasing Low','Increasing High','Theory Smooth','Theory Low','Theory High','Location','southeast')
xlabel('$Re$','Interpreter','latex')
ylabel('$C_d$','Interpreter','latex')


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
xlabel('Force ($N$)','Interpreter','latex')
ylabel('Voltage ($V$)','Interpreter','latex')
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
volt = [.002,.09,.239,.338,.456,.584,.762,.875,.726,.583,.455,.340,.239,.089,.002];
volt_std = [.017,.078,.077,.082,.081,.083,.072,.073,.071,.084,.084,.080,.079,.079,.017];

% Plotting raw data
figure('Name','Pressure Transducer Calibration')
plot(water_height,volt,'.','MarkerSize',12)
xlabel('Dynamic Pressure ($in H_{2}O$)','Interpreter','latex')
ylabel('Voltage ($V$)','Interpreter','latex')
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

function out = get_sensor_value(coeffs,volt)

out = (volt - coeffs(2)) ./ coeffs(1);

end

function Pa = inH2O_to_Pa(in)
Pa = 249.08890833333 .* in;
end

function [Cd,Re,force,velocity,press] = sphere_drag(rho_water,rho_air,mu,sphere,diam,press_coeffs,load_coeffs,man_load_volt,man_read)
fprintf('-----------------------------------------------------------------------')
str = sprintf('%s.csv',sphere);
sphere_dat = csvread(str,1);

load_volt = sphere_dat(:,1);
load_volt_std = sphere_dat(:,2);
press_volt = sphere_dat(:,3);
press_volt_std = sphere_dat(:,4);
Re_lab = sphere_dat(:,5);
Cd_lab = sphere_dat(:,6);

% Take the voltage readings and convert them to a pressure using the
% calibration coefficients. Then convert in H2O to Pa
press = inH2O_to_Pa(get_sensor_value(press_coeffs,press_volt));

% Compute the velocity of the flow using the dynamic pressure that was just
% calculated and the density of water.
velocity = sqrt((2.*press)./rho_water); % m/s?

% Convert voltage readings from load cell into forces
force = get_sensor_value(load_coeffs,load_volt);
S = (pi * (diam^2)) / 4;
Cd = (2.*force) ./ (rho_air .* (velocity.^2) .* S);
error_Cd = (abs(Cd_lab - Cd) ./ Cd_lab) .* 100;
fprintf('\nCd percent error for: %s\n',sphere)
for iCd = 1:length(error_Cd)
    fprintf('%.2f ',error_Cd(iCd))
end
fprintf('\n')

% Compute the Reynolds number
Re = (rho_air .* velocity .* diam) ./ (mu);
error_Re = (abs(Re_lab - Re) ./ Re_lab) .* 100;
fprintf('\nReynolds number percent error for: %s\n',sphere)
for iRe = 1:length(error_Cd)
    fprintf('%.2f ',error_Re(iRe))
end
fprintf('\n')


%%%  Compare the manual readings to the computer ones
man_press = inH2O_to_Pa(man_read);
man_force = get_sensor_value(load_coeffs,man_load_volt);
man_velocity = sqrt((2.*man_press)./rho_water);

man_Cd = (2.*man_force) ./ (rho_air .* (man_velocity.^2) .* S);
man_Re = (rho_air .* man_velocity .* diam) ./ (mu);

man_error_Cd = (abs(Cd_lab - man_Cd) ./ Cd_lab) .* 100;
man_error_Re = (abs(Re_lab - man_Re) ./ Re_lab) .* 100;

fprintf('\nMANUAL Cd percent error for: %s\n',sphere)
for iCd = 1:length(man_error_Cd)
    fprintf('%.2f ',man_error_Cd(iCd))
end
fprintf('\n')

fprintf('\nMANUAL Re percent error for: %s\n',sphere)
for iRe = 1:length(man_error_Re)
    fprintf('%.2f ',man_error_Re(iRe))
end
fprintf('\n')

end

function Cd_Theory = compute_theoretical_Cd(Re_smooth,Re_low,Re_high)

cd = @(Re) (24./Re) + ((2.6 .* (Re./5)) ./ ((1 + (Re./5).^1.52))) + ((.411 .* ((Re./263000).^-7.94)) ./ (1+((Re./263000).^-8))) + ((Re.^.8) ./ (461000));

Cd_Theory = zeros(3,length(Re_smooth));

Cd_Theory(1,:) = cd(Re_smooth);
Cd_Theory(2,:) = cd(Re_low);
Cd_Theory(3,:) = cd(Re_high);

end

function Ucd = cd_uncertainty(rho,vel,diam,drag,press)

Ucal = .00002236; %m
S = (pi * diam^2) / 4;
Us = sqrt(((pi * diam) / 2) * Ucal);
Ud = .70569; % N
Uu = .5.*sqrt(2./(rho.*press)).*inH2O_to_Pa(1.527);

Ucd = ones(size(press));
Ucd(:) = sqrt(((2./(rho.*(vel.^2).*S)).*Ud).^2 + (((4.*drag)./(rho.*(vel.^3).*S)).*Uu).^2 + (((2.*drag)./(rho.*(vel.^2).*(S.^2))).*Us).^2);
%Ucd(:) = mean(Ucd);
end















