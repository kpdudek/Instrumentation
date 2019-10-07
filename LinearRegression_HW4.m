function LinearRegression_HW4()
clear; clc; close all;
temp = [.5,9.8,21.2,39.7,76.4,99.2];
volt = [.02,.45,1.1,1.4,3.3,5];

[a0,a1] = lin_reg(volt,temp);
lin_temp = a0 + a1.*volt;

% Raw Data
figure()
hold on
plot(volt,temp,'.','MarkerSize',12)
plot(volt,lin_temp,'r','LineWidth',2)
legend('Raw Data','Linear Fit','Location','northwest')
xlim([-.2 5.5])
ylim([-5 110])
xlabel('Voltage (mV)')
ylabel('Temperature (^{\circ}C)')
end

function [a0,a1] = lin_reg(x,y)

N = length(x);

x_sq = zeros(size(x));
for iX = 1:N
    x_sq(iX) = x(iX)^2;
end

B = N*sum(x_sq) - (sum(x))^2;
a0 = (sum(x_sq)*sum(y) - sum(x)*sum((x.*y))) / B;
a1 = (N*sum((x.*y)) - sum(x)*sum(y)) / B;

f = polyfit(x,y,1);
disp([a0,f(2);a1,f(1)])

%%% Standard error of the fit
fit_error = zeros(size(y));
y_fit = a0 + a1.*x;
for iY = 1:N
    fit_error(iY) = (y_fit(iY) - y(iY))^2;
end
nu = N-2;
Syx = sqrt(sum(fit_error)/nu);
disp({'Syx: ',Syx})

end