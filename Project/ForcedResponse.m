close all; clear; clc;

motorSpeed = [12,20,25,28:8:100] .* (2*pi) ./ 60;
wResponse = [1.0438,1.3087,1.9564,2.2622,3.8406,3.7879,5.5970,5.6818,5.3004,5.3191,5.3060,5.3191,5.3191];

distPP = [.0824,.1455,.3723,.2110,.1986,.2099,.1536,.1637,.0114,.0119,.0112,.0117,.0105] ./ .5796;

figure('Name','Freq Response')
plot(motorSpeed,distPP,'b')
hold on

reducedPeak = .7 * max(distPP);
plot(motorSpeed,ones(size(motorSpeed)).*reducedPeak,'r')

xlabel('Frequency ($\frac{rad}{s}$)','Interpreter','latex')
ylabel('Magnitude')
