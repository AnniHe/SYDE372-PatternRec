clearvars; clc; close all;
%Model estimation 

load('lab2_1.mat');
load('lab2_2.mat');

dataAOne = OneDimData(a, 5, 1, 0);
dataAOne = dataAOne.getNormPdf();
dataAOne.plotAny();

%This is wack dawg
dataBOne = OneDimData(b, 0, 0, 1);
dataBOne = dataBOne.getNormPdf();
dataBOne.plotAny();

dataBTwo = OneDimData(b, 0, 0, 1);
dataBTwo = dataBTwo.getExpPdf()
dataBTwo.plotAny();

dataAThree = OneDimData(a, 0, 0, 0);
dataAThree = dataAThree.getUniPdf();
dataAThree.plotAny();

dataBThree = OneDimData(b, 0, 0, 0);
dataBThree = dataBThree.getUniPdf();
dataBThree.plotAny();

% ------------------------------
% 1D - Non-parametric Estimation
% -------------------------------
sigma1 = .1;
sigma2 = .4;
Parzen.plot(a',sigma1);
Parzen.plot(b',sigma1);
Parzen.plot(a',sigma2);
Parzen.plot(b',sigma2);

% ------------------------------
% 2D - Non-parametric Estimation
% -------------------------------
sigma = sqrt(400);
Parzen.plot(al,sigma);
Parzen.plot(bl,sigma);
Parzen.plot(cl,sigma);


