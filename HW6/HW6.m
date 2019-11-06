clear; close all; clc

t = 0:.01:2*pi;
tau = .15;

T = @(t) 115 + 12*sin(2*t);
V = @(t) 57.48*sin(2*t -.29) + 580;

figure('Name','Input Signal')
plot(t,T(t),'LineWidth',2)
xlabel('Time (s)')
ylabel('Temperature')

figure('Name','Output Signal')
plot(t,V(t),'LineWidth',2);%,t,ones(size(t)).*580,[pi/4,pi/4],[520,640])
xlabel('Time (s)')
ylabel('Voltage')


I = normalize(T(t));
V = normalize(V(t));

figure('Name','Normalized')
plot(t,I,t,V,'LineWidth',2)
hold on
xlabel('time (s)')
ylabel('Magnitude')
legend('Intput','Output')

t_I = t(find(I == max(I)));
t_V = t(find(V == max(V)));

t_lag = -t_V + t_I;
disp(t_lag)