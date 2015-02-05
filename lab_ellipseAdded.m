function main

    clear;
    clf();

    %CASE 1:

    Na = 200;
    muA = [5; 10];
    covarA = [8 0; 0 4];
    colourA = 'm'

    Nb = 200;
    muB = [10; 15];
    covarB = [8 0; 0 4];
    colourB = 'b'

    %CASE 2:
    Nc = 100;
    muC = [5; 10];
    covarC = [8 4; 4 40];
    colourC = 'r'

    Nd = 200;
    muD = [15; 10];
    covarD = [8 0; 0 8];
    colourD = 'g'

    Ne = 150;
    muE = [10; 5];
    covarE = [10 -5; -5 20;]
    colourE = 'c'

    %Generate bivariate normal distribution with mean vector and covariance
    %matrix
    dataA = gen_normal_distribution(Na, muA, covarA);
    figure(1);
    scatter(dataA(:,1), dataA(:,2), 5, colourA, 'filled');
    hold on;

    dataB = gen_normal_distribution(Nb, muB, covarB);
    scatter(dataB(:,1), dataB(:,2), 5, colourB, 'filled');
    hold on;

    figure(2);
    dataC = gen_normal_distribution(Nc, muC, covarC);
    scatter(dataC(:,1), dataC(:,2), 5, colourC, 'filled');
    hold on;
    dataD = gen_normal_distribution(Nd, muD, covarD);
    scatter(dataD(:,1), dataD(:,2), 5, colourD, 'filled');
    hold on;
    dataE = gen_normal_distribution(Ne, muE, covarE);
    scatter(dataE(:,1), dataE(:,2), 5, colourE, 'filled');
    hold on;

    %Find whitening transform W -> unit std.dev. contour
    [eigen_vec_A, eigen_diag_A] = get_eigenvalues(covarA);

    %plot ellipse
    figure(1);
    thetaA = atan2(eigen_vec_A(3),eigen_vec_A(1));
    draw_ellipse(muA(1), muA(2), thetaA, sqrt(eigen_diag_A(1)), sqrt(eigen_diag_A(4)), colourA);
    hold on;

    [eigen_vec_B, eigen_diag_B] = get_eigenvalues(covarB);
    thetaB = atan2(eigen_vec_B(3),eigen_vec_B(1));
    draw_ellipse(muB(1), muB(2), thetaB, sqrt(eigen_diag_B(1)), sqrt(eigen_diag_B(4)), colourB);
    hold on;

    figure(2);
    [eigen_vec_C, eigen_diag_C] = get_eigenvalues(covarC);
    thetaC = atan2(eigen_vec_C(3),eigen_vec_C(1));
    draw_ellipse(muC(1), muC(2), thetaC, sqrt(eigen_diag_C(1)), sqrt(eigen_diag_C(4)), colourC);
    hold on;

    [eigen_vec_D, eigen_diag_D] = get_eigenvalues(covarD);
    thetaD = atan2(eigen_vec_D(3),eigen_vec_D(1));
    draw_ellipse(muD(1), muD(2), thetaD, sqrt(eigen_diag_D(1)), sqrt(eigen_diag_D(4)), colourD);
    hold on;

    [eigen_vec_E, eigen_diag_E] = get_eigenvalues(covarE);
    thetaE = atan2(eigen_vec_E(3),eigen_vec_E(1));
    draw_ellipse(muE(1), muE(2), thetaE, sqrt(eigen_diag_E(1)), sqrt(eigen_diag_E(4)), colourE);
    hold on;

    whiten_A = eigen_diag_A^(-1/2) * transpose(eigen_vec_A);
    whiten_B = eigen_diag_B^(-1/2) * transpose(eigen_vec_B);
    whiten_C = eigen_diag_C^(-1/2) * transpose(eigen_vec_C);
    whiten_D = eigen_diag_D^(-1/2) * transpose(eigen_vec_D);
    whiten_E = eigen_diag_E^(-1/2) * transpose(eigen_vec_E);

    muA_whitened = whiten_A * muA;
    muB_whitened = whiten_B * muB;
    muC_whitened = whiten_C * muC;
    muD_whitened = whiten_D * muD;
    muE_whitened = whiten_E * muE;

    figure(1)
    %draw unit standard deviation contours of class A and B
    draw_circle(muA_whitened, colourA);
    draw_circle(muB_whitened, colourB);

    %MED: plot MED decision boundary with original means
    syms xA xB;
    ezplot(transpose(muB - muA)*[xA; xB] + 1/2 * (transpose(muA) * muA - transpose(muB) * muB) == 0, [-100, 100]);
    %GED: plot GED decision boundary with Whitening-transformed means
    syms xA_w xB_w;
    ezplot(transpose(muB_whitened - muA_whitened)*[xA_w; xB_w] + 1/2 * (transpose(muA_whitened) * muA_whitened - transpose(muB_whitened) * muB_whitened) == 0, [-100, 100]);

    %MAP:
    countAB = size(dataA) + size(dataB);
    pA = size(dataA)/countAB;
    pB = size(dataB)/countAB;
    MAP_2classes(covarA, covarB, muA, muB, pA, pB);

    %NN:
    %For every point in A, find shortest distance to a B point, find mid
    %distance

    figure(2)
    %draw unit standard deviation contours of class C, D and E
    draw_circle(muC_whitened, colourC);
    draw_circle(muD_whitened, colourD);
    draw_circle(muE_whitened, colourE);
    %do we find 3 decision boundaries of A/B, A/C and B/C because we have
    %three classes?
    %syms xC xD xE;
    %explot();

    %testing meshgrid stuff

    %draw meshgrid
    %     xMin = min(dataA(1));
    %     xMax = min(dataA(2));
    %     yMin = max(dataA(1));
    %     yMax = max(dataA(2));


    %     plot(meshgrid(-5:5,-5:5));
    %
    %     xrange = [-8 8];
    % yrange = [-8 8];
    % % step size for how finely you want to visualize the decision boundary.
    % inc = 0.1;
    %
    % % generate grid coordinates. this will be the basis of the decision
    % % boundary visualization.
    % (meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2)));

    %-----------------
    %
    %-----------------
    for i = 1:size(MED_distances)
    MED_Class = MED_Class()
    end

    print('help')
end

function MAP_2classes(s1, s2, mu1, mu2, p1, p2)
    q0 = inv(s1) - inv(s2)
    q1 = 2 * (transpose(mu2)* inv(s2) - transpose(mu1)* inv(s1))
    q2 = transpose(mu1) * inv(s1) * mu1 - transpose(mu2)* inv(s2) * mu2
    q3 = -1 * log(p2/p1)
    q4 = -1 * log(det(s1)/det(s2))

    syms x1 x2;
    ezplot([x1 x2] * q0 * [x1; x2] + q1 * [x1;x2] + q2 + 2 * q3 + q4 == 0, [-100, 100]);
end

function gauss = gen_normal_distribution(N, mu, covar)
    rand = randn(N,2);
    R = chol(covar);
    gauss = repmat(transpose(mu), length(rand), 1) + rand*R;
end

function [V,D] = get_eigenvalues(matrix)
    [V,D] = eig(matrix);
end

function cir = draw_circle(mu, colour)
    x = mu(1);
    y = mu(2);

    ang =0:0.01:2*pi;
    xp = 1*cos(ang);
    yp = 1*sin(ang);
    plot(x + xp, y + yp, colour);
    axis([-50 50 -50 50]);
    axis('square');
    hold on;
end

function draw_ellipse(x,y,theta,a,b, colour)

if nargin<5, error('Too few arguments to Plot_Ellipse.'); end;

np = 100;
ang = [0:np]*2*pi/np;
pts = [x;y]*ones(size(ang)) + [cos(theta) -sin(theta); sin(theta) cos(theta)]*[cos(ang)*a; sin(ang)*b];
plot( pts(1,:), pts(2,:), colour);

end

%Returns an array of [minimum, the class that was classified]
function MED_Class(point, meanArgs)
classDistances=[];

%compare mean A and mean B
for i = 1:length(meanArgs)
distance = sum(meanArgs(i) - point).^2
classDistances =[classDistances distance] %literally just adds another column
end

[minDistance, classType] = min(classDistances)
end


% gridSize -
% varargs - can take in multiple args, like varargs in java
%           random gaussian distrubution points
function prepareGrid (gridSize, varargs)
% Input gridsize and the (2D) data sets
minX_val=[];
minY_val=[];
maxX_val=[];
maxY_val=[];

for i = 1:length(varargs)
xVals = (:, 1);
yVals = (:, 2);

%bullshit matlab syntax
minX_val = [minX_vals min(xVals)];
maxX_val = [minX_vals max(xVals)];
minY_val = [minY_vals min(yVals)];
maxY_val = [maxY_vals max(yVals)];
end
end

% function med(mean1, mean2)
%     x1 = [0:1:100]
%     x2 = transpose(mean1-mean2) * x - 1/2*(transpose(mean1) * mean1 - mean2 * transpose(mean2)
%
