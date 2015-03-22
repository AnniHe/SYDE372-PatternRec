clear vars; clc; close all;

load feat.mat;

%for each image, find mean, covariance and whitening transform matrices 
[f2_whiten,f2_mean,f2_cov] = whitenTransform(f2);
[f8_whiten,f8_mean,f8_cov] = whitenTransform(f8);
[f32_whiten,f32_mean,f32_cov] = whitenTransform(f32);

%do classification and compute confusion matrices
f2_classified = MICD.classify(f2([1:2],:)', f2_mean, f2_whiten);
f2_confuse = confusionmat( f2t([3],:), f2_classified)
f2_error = 1 - trace(f2_confuse)/sum(f2_confuse(:))

f8_classified = MICD.classify(f8([1:2],:)', f8_mean, f8_whiten);
f8_confuse = confusionmat( f8t([3],:), f8_classified)
f8_error = 1 - trace(f8_confuse)/sum(f8_confuse(:))

f32_classified = MICD.classify(f32([1:2],:)', f32_mean, f32_whiten);
f32_confuse = confusionmat( f32t([3],:), f32_classified)
f32_error = 1 - trace(f32_confuse)/sum(f32_confuse(:))

%multiple textures
multi_class = [reshape(multf8(:,:,1),[],1) reshape(multf8(:,:,2),[],1)];
multi_f8 = MICD.classify( multi_class, f8_mean, f8_whiten);

cimage = reshape(multi_f8, size(multf8(:,:,1), 1), size(multf8(:,:,1), 2));
imagesc(cimage);
colormap(gray);


