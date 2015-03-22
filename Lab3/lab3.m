% -------------------
% Setup
% -------------------
clearvars; clc; close all;
image = readim('cloth.im') ;
imagesc(image) ;
colormap(gray) ;
load feat.mat
% figure;
% aplot(f2);
% figure;
% aplot(f8);
% figure;
% aplot(f32);





% ---------------------
% Unlabeled Clustering
% ---------------------
prototypes = Clustering.pickRandomPrototypes(f32);
Clustering.doKMeans(f32, prototypes);
Clustering.dofuzzyKMeans(f32, prototypes);