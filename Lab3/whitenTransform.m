function [set_whiten, set_mean, set_cov] = whitenTransform(data)

set_mean = [];
set_cov = [];
set_whiten = [];

for i = 1:10
    im_mean = mean(data([1:2],[15*(i-1)+1:15*i+i])');
    im_cov = cov(data([1],[15*(i-1)+1:15*i+i]),data([2],[15*(i-1)+1:15*i+i]));
    set_mean = [set_mean im_mean'];
    set_cov = [set_cov im_cov];
    [e_vec, e_diag] = eig(im_cov);
    whiten = e_vec^(1/2) * e_diag';
    set_whiten = [set_whiten whiten];
end

end