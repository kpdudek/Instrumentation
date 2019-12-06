close all; clc; clear;

FreeLow = load('FreeDecay_LowerRelease.mat');
FreeHigh = load('FreeDecay_UpperRelease.mat');
Forced100 = load('Forced_100.mat');


analyzeSensorData(Forced100.data,'F100');
% analyzeSensorData(FreeLow.data,'FreeLow');
% analyzeSensorData(FreeHigh.data,'FreeHigh');

