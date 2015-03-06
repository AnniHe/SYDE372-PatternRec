% function MAP_2classes(s1, s2, mu1, mu2, p1, p2)
%     q0 = inv(s1) - inv(s2)
%     q1 = 2 * (transpose(mu2)* inv(s2) - transpose(mu1)* inv(s1))
%     q2 = transpose(mu1) * inv(s1) * mu1 - transpose(mu2)* inv(s2) * mu2
%     q3 = -1 * log(p2/p1)
%     q4 = -1 * log(det(s1)/det(s2))
% 
%     syms x1 x2;
%     ezplot([x1 x2] * q0 * [x1; x2] + q1 * [x1;x2] + q2 + 2 * q3 + q4 == 0, [-100, 100]);
% end