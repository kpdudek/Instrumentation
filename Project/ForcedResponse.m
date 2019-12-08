close all; clear; clc;

% Motor knob setting converted to RPM, then to radians/sec
motorSpeed = [12,20,25,28:8:100] .* 2.67 .* (2*pi) ./ 60;

% Response frequency -- NOT USED
wResponse = [1.0438,1.3087,1.9564,2.2622,3.8406,3.7879,5.5970,5.6818,5.3004,5.3191,5.3060,5.3191,5.3191];

% Distance PP
distPP = [.0824,.1455,.3723,.2110,.1986,.2099,.1536,.1637,.0114,.0119,.0112,.0117,.0105];


figure('Name','Freq Response')
plot(motorSpeed,distPP,'b')
hold on

reducedPeak = .7 * max(distPP);
plot(motorSpeed,ones(size(motorSpeed)).*reducedPeak,'r')

xlabel('Frequency ($\frac{rad}{s}$)','Interpreter','latex')
ylabel('Amplitude ($\frac{X_{plate}PP}{x_{Base}}$)','Interpreter','latex')


wr = 6.99;
wu = 7.57;
wl = 6.30;
Q_fwhm = wr/(wu-wl);
zeta = 1/(2*Q_fwhm);
fprintf('Q FWHM: %5.5f\n',Q_fwhm)
fprintf('Zeta: %5.5f\n', zeta)


figure()

plot(motorSpeed,distPP,'b')
hold on

wn = 6.99;
z = .09084;
M_theory = @(w) 1 ./ ((1-(w/wn).^2).^2 + (2.*z.*(w/wn)).^2).^.5;

Mt = M_theory(motorSpeed).*.08;
plot(motorSpeed,Mt,'g')

xlabel('Frequency ($\frac{rad}{s}$)','Interpreter','latex')
ylabel('Amplitude ($\frac{X_{plate}PP}{x_{Base}}$)','Interpreter','latex')
legend('Experimental','Theoretical')
