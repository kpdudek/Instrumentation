function WheatstoneBridge()

question1();

question2();

question3();

end

function question1()

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
plot(full_ax,deltaR1_over_R1,Vac_over_V,'.','MarkerSize',12)
hold on

%%% Theory
Vac_over_V_non_lin = deltaR1_over_R1 ./ (4+2*deltaR1_over_R1);
plot(full_ax,deltaR1_over_R1,Vac_over_V_non_lin)
hold on

Vac_over_V_lin = deltaR1_over_R1 ./ 4;
plot(full_ax,deltaR1_over_R1,Vac_over_V_lin)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot -100 ohm to 100 ohm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_R1_linear = delta_R1(find(abs(delta_R1) < 100)); 
Vac_Linear = Vac_Full_Range(find(abs(delta_R1) < 100)); 

linear_deltaR1_over_R1 = delta_R1_linear./DRB_balanced;
linear_Vac_over_V = (Vac_Linear./1000)./V_applied;

lin_fig = figure('Name','Linear Regime');
plot(linear_deltaR1_over_R1,linear_Vac_over_V,'.','MarkerSize',12)
hold on

%%% Theory
linear_Vac_over_V_non_lin = linear_deltaR1_over_R1 ./ (4+2*linear_deltaR1_over_R1);
plot(linear_deltaR1_over_R1,linear_Vac_over_V_non_lin)
hold on

linear_Vac_over_V_lin = linear_deltaR1_over_R1 ./ 4;
plot(linear_deltaR1_over_R1,linear_Vac_over_V_lin)
hold on
end

function question2()



end

function question3()
F = 2.125;
V = 5; % V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quarter bridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quarter_micrometer = in2m([.760,.733,.717,.691]); %in
quarter_unamp = [-.321,-.103,-.0299,-.0472]; %mV
quarter_amp = [-.1731,-.1522,-.1432,-.1351]; %V
quarter_zero_unamp = -.320; %mV
quarter_zero_amp = -.1735; %V

quarter_unamp_strain = (4.*(quarter_unamp-quarter_zero_unamp)./1000) ./ (F * V);
quarter_amp_strain = (4.*(quarter_amp-quarter_zero_amp)./100) ./ (F * V);
disp('Quarter bridge strain')
disp(quarter_unamp_strain)
disp(quarter_amp_strain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Half bridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
half_micrometer = in2m([.760,.738,.717,.696]); %in
half_unamp = [-3.159,-3.086,-2.852,-2.746]; %mV
half_amp = [-.466,-.453,-.432,-.414]; %V
half_zero_unamp = -3.720; %mV
half_zero_amp = -.521; %V

half_unamp_strain = (4.*(half_unamp-half_zero_unamp)./1000) ./ (F * V);
half_amp_strain = (4.*(half_amp-half_zero_amp)./100) ./ (F * V);
disp('Half bridge strain')
disp(half_unamp_strain)
disp(half_amp_strain)

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
disp('theory_strain')
disp(theory_strain)

%%% Deflection
y = (-P ./ (E .* I)) .* (((L .* m.^2)./2)-(m.^3./2));
disp('Deflections')
disp(y)
disp(quarter_micrometer)
disp(half_micrometer)
end

function m = in2m(in)
m = in .* .0254;
end