%three feature matrices
%find mean and covariance matrix
%
clear vars; clc; close all;

load feat.mat;

%for each image, find mean, covariance and whitening transform matrices 
[f2_whiten,f2_mean,f2_cov] = whitenTransform(f2);
[f8_whiten,f8_mean,f8_cov] = whitenTransform(f8);
[f32_whiten,f32_mean,f32_cov] = whitenTransform(f32);

f2_classified = MICD.classify(f2([1:2],:)', f2_mean, f2_whiten);
f2_confuse = confusionmat( f2t([3],:), f2_classified)

f8_classified = MICD.classify(f8([1:2],:)', f8_mean, f8_whiten);
f8_confuse = confusionmat( f8t([3],:), f8_classified)

f32_classified = MICD.classify(f32([1:2],:)', f32_mean, f32_whiten);
f32_confuse = confusionmat( f32t([3],:), f32_classified)