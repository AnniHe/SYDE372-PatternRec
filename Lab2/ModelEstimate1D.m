clear all; clc; close all;
%Model estimation 

load('lab2_1.mat');

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

% m_a = mean(a);
% m_b = mean(b);
% 
% var_a = var(a);
% var_b = var(b);
% 
% range_a = [min(a): 0.1:max(a)];
% pdf_a = normpdf(range_a, m_a, var_a);
% pdf_a_true = normpdf(range_a, 5, 1);
% 
% lambda_b = inv(m_b);
% range_b = [min(b): 0.1:max(b)];
% pd_b = fitdist(b', 'Exponential');
% pdf_b = pdf(pd_b, range_b);
% 
% 
% figure(1);
% plot(range_a, pdf_a, 'b');
% hold on;
% plot(range_a, pdf_a_true, 'r');
% title('Gaussian Parametric Estimation: Data A');
% legend('Sample Distribution', 'True Distribution');
% 
% figure(2);
% plot(range_b, pdf_b, 'b');
% hold on;
% plot(range_b, pdf_b, 'r');
% title('Gaussian Parametric Estimation: Data B');
% legend('Sample Distribution', 'True Distribution');





