close all; clear; clc;

load('Data.mat')
m = .290; %kg

% Displacement response
figure('Name','Displacement')
force = (HW8Springconstantdata(:,1)./ 1000 ).* 9.8;
delta = HW8Springconstantdata(:,2) ./ 100;
plot(force,delta,'b.','MarkerSize',12)
hold on;
coeff = polyfit(force,delta,1);
linFit = polyval(coeff,force);
plot(force,linFit,'r','LineWidth',2)

xlabel('Force (N)')
ylabel('Displacement (m)')

K = abs(1/coeff(1));
wnat = sqrt(K/m);
fprintf('Spring constant: %5.5f\n',K)
fprintf('Natural Frequency: %5.5f\n',wnat)


% Reponse Data
figure('Name','Response Data')
omega = HW8responsedata(:,1);
Vpp = HW8responsedata(:,2);

Vpp_max = max(Vpp);

plot(omega,Vpp,'bx-','MarkerSize',12)
hold on
plot(omega,ones(size(Vpp)).*(.7*Vpp_max),'r')

xlabel('Frequency ($\frac{rad}{s}$)','Interpreter','latex')
ylabel('$V_{pp}$ (V)','Interpreter','latex')

wr = 15.235;
wu = 16.829;
wl = 13.236;
Q_fwhm = wr/(wu-wl);
zeta = 1/(2*Q_fwhm);
fprintf('Q FWHM: %5.5f\n',Q_fwhm)
fprintf('Zeta: %5.5f\n', zeta)

% Free Data
figure('Name','Free Decay')
time = HW8freedecaydata(:,1);
volt = HW8freedecaydata(:,2);

plot(time,volt,'b')
hold on 
peakT = [.172,.566,.98,1.4];
peakV = [.425,.435,.479,.513];
plot(peakT,peakV,'ro','MarkerSize',8)

figure('Name','Exponential Fit')
ln_env = -1*log(peakV);
plot(peakT,ln_env,'b.','MarkerSize',12)
hold on 
lnCoeff = polyfit(peakT,ln_env,1);
linVolt = polyval(lnCoeff,peakT);
plot(peakT,linVolt,'r')

fprintf('Linear slope: %5.5f\n',lnCoeff(1))


Td = mean([.566-.172, .98-.566, 1.4-.98]);
fprintf('Td: %5.5f\n',Td)

wd = 2*pi / Td;
fprintf('Wd: %5.5f\n',wd)












