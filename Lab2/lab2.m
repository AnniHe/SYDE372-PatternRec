clearvars; clc; close all;

load('lab2_1.mat');
load('lab2_2.mat');
dataAOne = OneDimData(a, 5, 1, 0);
dataBOne = OneDimData(b, 1, 0, 1);

% -----------------------------------
% Dataset A - Parametric Estimation Gaussian
% Question 2.1
% -----------------------------------
figure();
hold on;
dataAOne.plotEstimatedGaussian(); 
dataAOne.plotTrueGaussian();
hold off;
% -----------------------------------
% Dataset B - Parametric Estimation Gaussian
% Question 2.1
% -----------------------------------
figure();
hold on;
dataBOne.plotEstimatedGaussian();
dataBOne.plotTrueExponential();
hold off;
% -----------------------------------
% Dataset A - Parametric Estimation Exponential
% Question 2.2
% -----------------------------------
figure();
hold on;
dataAOne.plotEstimatedExponential(); %%%%why is this sooo low
dataAOne.plotTrueGaussian();
hold off;
% -----------------------------------
% Dataset B - Parametric Estimation Exponential
% Question 2.2
% -----------------------------------
figure();
hold on;
dataBOne.plotEstimatedExponential(); %%%%why is this sooo low
dataBOne.plotTrueExponential();
hold off;
% -----------------------------------
% Dataset A - Parametric Estimation Uniform
% Question 2.3
% -----------------------------------
figure();
hold on;
dataAOne.plotEstimatedUniform();
dataAOne.plotTrueGaussian();
hold off;
% -----------------------------------
% Dataset B - Parametric Estimation Uniform
% Question 2.3
% -----------------------------------
figure();
hold on;
dataBOne.plotEstimatedUniform();
dataBOne.plotTrueExponential();
hold off;

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



