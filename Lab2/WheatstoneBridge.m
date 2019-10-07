function WheatstoneBridge()
clear; clc;

question1_2();

question3_4();

end

function question1_2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DRB_balanced = 10002; % K ohm
Vac_balanced = -.199; % mV

V_applied = 15; %V

DRB_decreasing = [9982,9960,9940,9920,9900,9800,9700,9601,9501,...
    9401,8902,8403,7704,7405,6908]; % ohm
Vac_decreasing = [-8.75,-16.955,-24.412,-32.151,-39.541,-77.534,-115.982,...
    -154.487,-194.142,-233.665,-437.113,-652.084,-879.422,-1119.482,-1372.65]; % mV

DRB_increasing = [10022,10042,10062,10082,10102,10202,10301,10401,...
    10501,10601,11100,11599,12098,12597,13094]; % ohm
Vac_increasing = [5.811,13.208,20.476,27.809,35.156,71.877,108.084,...
    144.251,179.731,215.285,387.251,551.084,707.481,856.785,999.264]; %mV

DRB_Full_Range = [DRB_decreasing,DRB_increasing];
[DRB_Full_Range,order] = sort(DRB_Full_Range);
Vac_Full_Range = [Vac_decreasing,Vac_increasing];
Vac_Full_Range = Vac_Full_Range(order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot -3000 ohm to 3000 ohm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_R1 = DRB_Full_Range - DRB_balanced;
deltaR1_over_R1 = delta_R1./DRB_balanced;
Vac_over_V = (Vac_Full_Range./1000)./V_applied;

full_fig = figure('Name','Full Delta R');
full_ax = axes(full_fig);
plot(full_ax,deltaR1_over_R1,Vac_over_V,'.','MarkerSize',13)
hold on

%%% Theory
Vac_over_V_non_lin = deltaR1_over_R1 ./ (4+2*deltaR1_over_R1);
plot(full_ax,deltaR1_over_R1,Vac_over_V_non_lin,'LineWidth',1)
hold on

Vac_over_V_lin = deltaR1_over_R1 ./ 4;
plot(full_ax,deltaR1_over_R1,Vac_over_V_lin,'LineWidth',1)
hold on

legend('Experimental','Nonlinear Theory','Linear Theory','Location','northwest')
xlabel('$\frac{\Delta R1}{R_{10}}$','Interpreter','latex','FontSize',14)
ylabel('$\frac{V_{ac}}{V}$','Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot -100 ohm to 100 ohm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_R1_linear = delta_R1(find(abs(delta_R1) < 100)); 
Vac_Linear = Vac_Full_Range(find(abs(delta_R1) < 100)); 

linear_deltaR1_over_R1 = delta_R1_linear./DRB_balanced;
linear_Vac_over_V = (Vac_Linear./1000)./V_applied;

lin_fig = figure('Name','Linear Regime');
plot(linear_deltaR1_over_R1,linear_Vac_over_V,'.','MarkerSize',13)
hold on

%%% Theory
linear_Vac_over_V_non_lin = linear_deltaR1_over_R1 ./ (4+2*linear_deltaR1_over_R1);
plot(linear_deltaR1_over_R1,linear_Vac_over_V_non_lin)
hold on

linear_Vac_over_V_lin = linear_deltaR1_over_R1 ./ 4;
plot(linear_deltaR1_over_R1,linear_Vac_over_V_lin)
hold on

legend('Experimental','Nonlinear Theory','Linear Theory','Location','northwest')
xlabel('$\frac{\Delta R1}{R_{10}}$','Interpreter','latex','FontSize',14)
ylabel('$\frac{V_{ac}}{V}$','Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute linearity error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linear_Vac = Vac_over_V_lin .* V_applied;
disp('Linear Vac')
disp(linear_Vac)
disp('Vac')
disp(Vac_Full_Range./1000)
percent_error = (((Vac_Full_Range./1000) - linear_Vac)./linear_Vac).*100;
disp('Linearity Error:')
disp(percent_error)
end

function question3_4()
F = 2.125;
V = 5; % V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quarter bridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quarter_micrometer = in2m(abs([.760,.733,.717,.691]-.781)); %in
quarter_unamp = [-.321,-.103,-.0299,-.0472]; %mV
quarter_amp = [-.1731,-.1522,-.1432,-.1351]; %V
quarter_zero_unamp = -5.94; %mV
quarter_zero_amp = -.477; %V

quarter_unamp_strain = (4.*(quarter_unamp-quarter_zero_unamp)./1000) ./ (F * V);
quarter_amp_strain = (4.*(quarter_amp-quarter_zero_amp)./100) ./ (F * V);
disp('Quarter bridge strain:')
fprintf('Unamplified | %.3e %.3e %.3e %.3e\n',quarter_unamp_strain(1),quarter_unamp_strain(2),quarter_unamp_strain(3),quarter_unamp_strain(4))
fprintf('Amplified   | %.3e %.3e %.3e %.3e\n\n',quarter_amp_strain(1),quarter_amp_strain(2),quarter_amp_strain(3),quarter_amp_strain(4))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Half bridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
half_micrometer = in2m(abs([.760,.738,.717,.696]-.781)); %in
half_unamp = [-3.159,-3.086,-2.852,-2.746]; %mV
half_amp = [-.466,-.453,-.432,-.414]; %V
half_zero_unamp = -3.720; %mV
half_zero_amp = -.521; %V

half_unamp_strain = (2.*(half_unamp-half_zero_unamp)./1000) ./ (F * V);
half_amp_strain = (2.*(half_amp-half_zero_amp)./100) ./ (F * V);
disp('Half bridge strain:')
fprintf('Unamplified | %.3e %.3e %.3e %.3e\n',half_unamp_strain(1),half_unamp_strain(2),half_unamp_strain(3),half_unamp_strain(4))
fprintf('Amplified   | %.3e %.3e %.3e %.3e\n\n',half_amp_strain(1),half_amp_strain(2),half_amp_strain(3),half_amp_strain(4))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Theoretical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = .0317;
m = .2325;
L = .2585;
b = .02579;
h = .00347;
c = h/2;
P = [.05,.1,.15,.2] .* 9.8; %Kg
E = 69 * 10^9;
I = (b * h^3) / 12;

%%% Strain
theory_strain = (P * (L - a) * c) ./ ( E .* I);
disp('Theory Strain:')
fprintf('%.3e %.3e %.3e %.3e\n\n',theory_strain(1),theory_strain(2),theory_strain(3),theory_strain(4))

%%% Deflections
y = (-P ./ (E .* I)) .* (((L .* m.^2)./2)-(m.^3./6));
disp('Theory Deflection:')
fprintf('%.3e %.3e %.3e %.3e\n\n',y(1),y(2),y(3),y(4))

disp('Measured Deflections:')
fprintf('Quarter Micrometer | %.3e %.3e %.3e %.3e\n',quarter_micrometer(1),quarter_micrometer(2),quarter_micrometer(3),quarter_micrometer(4))
fprintf('Half Micrometer    | %.3e %.3e %.3e %.3e\n\n',half_micrometer(1),half_micrometer(2),half_micrometer(3),half_micrometer(4))

end

function m = in2m(in)
m = in .* .0254;
end