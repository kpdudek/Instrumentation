close all; clear; clc;

wn = 180000;
K = .1; %V/g
o = .9; %V
a = 7; %g

% For A
zA = .05;
% For B
zB = .3;


M = @(w,z) 1./sqrt((1-(w./wn).^2).^2 + (2.*z.*(w./wn)).^2);
phi = @(w,z) -1.*atan((2.*z.*(w./wn))./(1-(w./wn).^2));

w = 0:2:400000;

magA = M(w,zA);
phiA = phi(w,zA);

magB = M(w,zB);
phiB = phi(w,zB);

% MAGNITUDE PLOT
figure()
plot(w,magA,w,magB,'LineWidth',2)
xlabel('Frequency (w)')
ylabel('M(w)')
hold on
plot(w,ones(1,length(w)).*1.41,'r--',w,ones(1,length(w)).*.708,'r--')
legend('A-1','B-1')

% PHASE PLOT
figure()
plot(w,phiA,w,phiB,'LineWidth',2)
xlabel('Frequency (w)')
ylabel('$\Phi(w)$','Interpreter','latex')
legend('A-1','B-1')

% TIME LAG
input = @(t,w) 7.*sin(w*t);
output = @(t,w,z) K.*a.*M(w,z).*sin(w*t+phi(w,z))+o;

t = 0:.0001:.005;

accelA = input(t,280000);
VoltA = output(t,280000,zA);

accelB = input(t,260000);
VoltB = output(t,260000,zB);

figure('Name','A')
plot(t,accelA,t,VoltA)
xlabel('Time (s)')
legend('Input','Output')

figure('Name','B')
plot(t,accelB,t,VoltB)
xlabel('Time (s)')
legend('Input','Output')

% ERROR PLOT
figure()
errorA = abs(magA-1);
errorB = abs(magB-1);

plot(w,errorB,'LineWidth',2)
legend('B-1')






