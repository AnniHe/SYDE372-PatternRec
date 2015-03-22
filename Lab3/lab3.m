% -------------------
% Setup
% -------------------
clearvars; clc; 
image = readim('cloth.im') ;
imagesc(image) ;
colormap(gray) ;
load feat.mat
aplot(f2);





% -------------------
% Unlabeled Clustering
% -------------------
Clustering.doKMeans(f32);