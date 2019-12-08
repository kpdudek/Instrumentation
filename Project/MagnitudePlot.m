clear; close all; clc;

Freq=0:0.01:3;
D1=.010;
D2=5;
D3=2;
D4=1;
D5=0.707;
D6=0.5;
D7=0.3;
D8=0.2;
D9=0;

M1=1./((((1-Freq.^2).^2)+(2.*D1.*Freq).^2).^(1/2));

figure();
plot(Freq,M1,'b');
