clear; clc;

zeta1 = .09084;
zeta2 = .01078;

k = 77.7;
m = (126.7 + 296.2)/1000;

C1 = zeta1 * 2 * sqrt(k * m);
C2 = zeta2 * 2 * sqrt(k * m);

fprintf('C1: %5.5f\n',C1)
fprintf('C2: %5.5f\n',C2)

