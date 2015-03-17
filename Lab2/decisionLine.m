function y = decisionLine(point1, point2,x)

mid = [(point1(1) + point2(1))/2 (point1(2) + point2(2))/2]';

slope = -1/((point2(2) - point1(2))/ (point2(1) - point1(1)));

interc = mid(2) - slope * mid(1);

y = slope * x + interc;

end