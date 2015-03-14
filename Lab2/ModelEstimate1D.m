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
muA = 5;
range = [0:.1:10];
Parzen.plot(a',sigma1);
Plot.plotNormPdf(range, muA, sigma1, 'r');
axis tight;
Parzen.plot(b',sigma1);
Parzen.plot(a',sigma2);
Parzen.plot(b',sigma2);

% ------------------------------
% 2D - Non-parametric Estimation
% Question 3.2 - INCOMPLETE
% -------------------------------
sigma = sqrt(400);

%Create a matrix window with a Gaussian shape
%http://dali.feld.cvut.cz/ucebna/matlab/toolbox/images/fspecial.html
window = fspecial('gaussian', [25 25], sigma);

%Resolution (step size)
%Determines the spatial step between PDF estimates
lowx = min([min(al(:,1)), min(bl(:,1)), min(cl(:,1))]);
lowy = min([min(al(:,2)), min(bl(:,2)), min(cl(:,2))]);
highx = max([max(al(:,1)), max(bl(:,1)), max(cl(:,1))]);
highy = max([max(al(:,2)), max(bl(:,2)), max(cl(:,2))]);
step = 1;
resolution = [step (lowx) (lowy) (highx) (highy)];

[pdfA, xA, yA] = Parzen.parzen2D(al, resolution, window);
[pdfB, xB, yB] = Parzen.parzen2D(bl, resolution, window);
[pdfC, xC, yC] = Parzen.parzen2D(cl, resolution, window);

%Create a grid
grid = zeros(length(xA)-1,length(yA)-1);

%Split grid into discrete squares
xVals = [lowx:step:highx];
yVals = [lowy:step:highy];

%classify using ML
for i = 1:xA
    for j = 1:yA
        grid(i,j) = ML.classify(i,j, pdfA, pdfB, pdfC);
    end
end

contour(xVals,yVals,grid',3,'k');

% Parzen.plot(al,sigma);
% Parzen.plot(bl,sigma);
% Parzen.plot(cl,sigma);


