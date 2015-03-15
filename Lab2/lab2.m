clearvars; clc; close all;

load('lab2_1.mat');
load('lab2_2.mat');
dataAOne = OneDimData(a, 5, 1, 0);
dataBOne = OneDimData(b, 1, 0, 1);

% -----------------------------------
% 2A - Parametric Estimation Gaussian
% Question 2.1
% -----------------------------------
dataAOne.plotGaussian();
dataBOne.plotGaussian(); %problem with this is that variance is set to 0
% -----------------------------------
% 2B - Parametric Estimation Exponential
% Question 2.2
% -----------------------------------
dataAOne.plotExponential();
dataBOne.plotExponential();
% -------------------------------
% 3C - Parametric Estimation Uniform
% Question 2.3
% -------------------------------
dataAOne.plotUniform();
dataBOne.plotUniform();

dataBOne = dataBOne.getNormPdf();
dataBOne.plotAny();




%This is wack dawg
%this doesn't work because the following line returns a bunch of NaN:
%OneD.pdf_true = normpdf(OneD.range, OneD.mu_true, OneD.sigma_true); 


dataBTwo = OneDimData(b, 0, 0, 1);
dataBTwo = dataBTwo.getExpPdf()
dataBTwo.plotAny();

dataAThree = OneDimData(a, 0, 0, 0);
dataAThree = dataAThree.getUniPdf();
dataAThree.plotAny();

dataBThree = OneDimData(b, 0, 0, 0);
dataBThree = dataBThree.getUniPdf();
dataBThree.plotAny();

% -------------------------------
% 1D - Non-parametric Estimation
% Question 2.4
% -------------------------------
sigma1 = .1;
sigma2 = .4;
muA = 5;
range = [0:.1:10];
Parzen.plot(range, a', sigma1);
Plot.plotNormPdf(range, muA, 1, 'r');
Plot.applyCase(1);
Parzen.plot(range, a', sigma2);
Plot.plotNormPdf(range, muA, 1, 'r');
Plot.applyCase(1);

Parzen.plot(range, b',sigma1);
Plot.plotExponentialPdf(range, 1, 'r');
Plot.applyCase(1);
Parzen.plot(range, b',sigma2);
Plot.plotExponentialPdf(range, 1, 'r');
Plot.applyCase(1);

% -------------------------------
% 2D - Non-parametric Estimation
% Question 3.2
% -------------------------------
Parzen.plot2D(al,bl,cl,20);



