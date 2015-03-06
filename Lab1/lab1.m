function main

    clear;
    clf();

    gridSize = 0.1;
    k = 5;
    
    %CASE 1:

    Na = 200;
    muA = [5; 10];
    covarA = [8 0; 0 4];
    colourA = 'm';

    Nb = 200;
    muB = [10; 15];
    covarB = [8 0; 0 4];
    colourB = 'b';

    %CASE 2:
    Nc = 100;
    muC = [5; 10];
    covarC = [8 4; 4 40];
    colourC = 'r';

    Nd = 200;
    muD = [15; 10];
    covarD = [8 0; 0 8];
    colourD = 'g';

    Ne = 150;
    muE = [10; 5];
    covarE = [10 -5; -5 20;];
    colourE = 'c';

    %Generate bivariate normal distribution with mean vector and covariance
    %matrix
    dataA = gen_normal_distribution(Na, muA, covarA);
    dataATest = gen_normal_distribution(Na, muA, covarA);
    figure(1);
    scatter(dataA(:,1), dataA(:,2), 5, colourA, 'filled');
    hold on;
    figure(3);
    scatter(dataA(:,1), dataA(:,2), 5, colourA, 'filled');
    hold on;

    dataB = gen_normal_distribution(Nb, muB, covarB);
    dataBTest = gen_normal_distribution(Nb, muB, covarB);
    figure(1);
    scatter(dataB(:,1), dataB(:,2), 5, colourB, 'filled');
    hold on;
    figure(3);
    scatter(dataB(:,1), dataB(:,2), 5, colourB, 'filled');
    hold on;


    figure(2);
    dataC = gen_normal_distribution(Nc, muC, covarC);
    dataCTest = gen_normal_distribution(Nc, muC, covarC);
    scatter(dataC(:,1), dataC(:,2), 5, colourC, 'filled');
    hold on;
    dataD = gen_normal_distribution(Nd, muD, covarD);
    dataDTest = gen_normal_distribution(Nd, muD, covarD);
    scatter(dataD(:,1), dataD(:,2), 5, colourD, 'filled');
    hold on;
    dataE = gen_normal_distribution(Ne, muE, covarE);
    dataETest = gen_normal_distribution(Ne, muE, covarE);
    scatter(dataE(:,1), dataE(:,2), 5, colourE, 'filled');
    hold on;
    
    figure(4);
    scatter(dataC(:,1), dataC(:,2), 5, colourC, 'filled');
    hold on;
    scatter(dataD(:,1), dataD(:,2), 5, colourD, 'filled');
    hold on;
    scatter(dataE(:,1), dataE(:,2), 5, colourE, 'filled');
    hold on;

    

    %Find whitening transform W -> unit std.dev. contour
    [eigen_vec_A, eigen_diag_A] = get_eigenvalues(covarA);

    %plot ellipse
    figure(1);
    thetaA = atan2(eigen_vec_A(3),eigen_vec_A(1));
    Draw.draw_ellipse(muA(1), muA(2), thetaA, sqrt(eigen_diag_A(1)), sqrt(eigen_diag_A(4)), colourA);
    hold on;

    [eigen_vec_B, eigen_diag_B] = get_eigenvalues(covarB);
    thetaB = atan2(eigen_vec_B(3),eigen_vec_B(1));
    Draw.draw_ellipse(muB(1), muB(2), thetaB, sqrt(eigen_diag_B(1)), sqrt(eigen_diag_B(4)), colourB);
    hold on;
    
    figure(3);
    Draw.draw_ellipse(muA(1), muA(2), thetaA, sqrt(eigen_diag_A(1)), sqrt(eigen_diag_A(4)), colourA);
    hold on;
    Draw.draw_ellipse(muB(1), muB(2), thetaB, sqrt(eigen_diag_B(1)), sqrt(eigen_diag_B(4)), colourB);
    hold on;

    figure(2);
    [eigen_vec_C, eigen_diag_C] = get_eigenvalues(covarC);
    thetaC = atan2(eigen_vec_C(3),eigen_vec_C(1));
    Draw.draw_ellipse(muC(1), muC(2), thetaC, sqrt(eigen_diag_C(1)), sqrt(eigen_diag_C(4)), colourC);
    hold on;

    [eigen_vec_D, eigen_diag_D] = get_eigenvalues(covarD);
    thetaD = atan2(eigen_vec_D(3),eigen_vec_D(1));
    Draw.draw_ellipse(muD(1), muD(2), thetaD, sqrt(eigen_diag_D(1)), sqrt(eigen_diag_D(4)), colourD);
    hold on;

    [eigen_vec_E, eigen_diag_E] = get_eigenvalues(covarE);
    thetaE = atan2(eigen_vec_E(3),eigen_vec_E(1));
    Draw.draw_ellipse(muE(1), muE(2), thetaE, sqrt(eigen_diag_E(1)), sqrt(eigen_diag_E(4)), colourE);
    hold on;
    
    figure(4);
    Draw.draw_ellipse(muC(1), muC(2), thetaC, sqrt(eigen_diag_C(1)), sqrt(eigen_diag_C(4)), colourC);
    hold on;
    Draw.draw_ellipse(muD(1), muD(2), thetaD, sqrt(eigen_diag_D(1)), sqrt(eigen_diag_D(4)), colourD);
    hold on;
    Draw.draw_ellipse(muE(1), muE(2), thetaE, sqrt(eigen_diag_E(1)), sqrt(eigen_diag_E(4)), colourE);
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

    %draw unit standard deviation contours
    figure(1)
    Draw.draw_circle(muA_whitened, colourA);
    Draw.draw_circle(muB_whitened, colourB);
    figure(2)
    Draw.draw_circle(muC_whitened, colourC);
    Draw.draw_circle(muD_whitened, colourD);
    Draw.draw_circle(muE_whitened, colourE);
    figure(3)
    Draw.draw_circle(muA_whitened, colourA);
    Draw.draw_circle(muB_whitened, colourB);
    figure(4)
    Draw.draw_circle(muC_whitened, colourC);
    Draw.draw_circle(muD_whitened, colourD);
    Draw.draw_circle(muE_whitened, colourE);

    AB_means = [muA muB];
    AB_tmeans = [muA_whitened muB_whitened];
    AB_transforms = [whiten_A whiten_B];
    AB_probs = [Na/(Na+Nb) Nb/(Na+Nb)];
    AB_covars = [covarA covarB];
    
    CDE_means = [muC muD muE];
    CDE_tmeans = [muC_whitened muD_whitened muE_whitened];
    CDE_transforms = [whiten_C whiten_D whiten_E];
    CDE_probs = [Nc/(Nc+Nd+Ne) Nd/(Nc+Nd+Ne) Ne/(Nc+Nd+Ne)];
    CDE_covars = [covarC covarD covarE];
    
    %Grid Helper
    [xVals_AB, yVals_AB, grid_AB] = Plotter.prepareGrid(gridSize, dataA, dataB);
    [xVals_AB_Test, yVals_AB_Test, grid_AB_Test] = Plotter.prepareGrid(gridSize, dataATest, dataBTest);
    [xVals_CDE, yVals_CDE, grid_CDE] = Plotter.prepareGrid(gridSize, dataC, dataD, dataE);
    [xVals_CDE_Test, yVals_CDE_Test, grid_CDE_Test] = Plotter.prepareGrid(gridSize, dataCTest, dataDTest, dataETest);

    figure(1);
    MED_AB = Plotter.plot_MED(grid_AB, xVals_AB, yVals_AB, AB_means);
    GED_AB = Plotter.plot_GED(grid_AB, xVals_AB, yVals_AB, AB_tmeans, AB_transforms);
    MAP_AB = Plotter.plot_MAP(grid_AB, xVals_AB, yVals_AB, AB_means, AB_covars, AB_probs);
    figure(2);
    MED_CDE = Plotter.plot_MED(grid_CDE, xVals_CDE, yVals_CDE, CDE_means);
    GED_CDE = Plotter.plot_GED(grid_CDE, xVals_CDE, yVals_CDE, CDE_tmeans, CDE_transforms);
 	MAP_CDE = Plotter.plot_MAP(grid_CDE, xVals_CDE, yVals_CDE, CDE_means, CDE_covars, CDE_probs);
    figure(3);
    NN_AB  = Plotter.plot_nn(1, grid_AB_Test, xVals_AB_Test, yVals_AB_Test, dataA, dataB);
    KNN_AB = Plotter.plot_nn(k, grid_AB_Test, xVals_AB_Test, yVals_AB_Test, dataA, dataB);
    figure(4);
    NN_CDE  = Plotter.plot_nn(1, grid_CDE_Test, xVals_CDE_Test, yVals_CDE_Test, dataC, dataD, dataE);
 	KNN_CDE = Plotter.plot_nn(k, grid_CDE_Test, xVals_CDE_Test, yVals_CDE_Test, dataC, dataD, dataE);

    ErrorAnalysis.confusion(xVals_AB, yVals_AB, MED_AB, 'MED', dataA, dataB);
    ErrorAnalysis.confusion(xVals_AB, yVals_AB, GED_AB, 'GED', dataA, dataB);
    ErrorAnalysis.confusion(xVals_AB, yVals_AB, MAP_AB, 'MAP', dataA, dataB);
    ErrorAnalysis.confusion(xVals_AB_Test, yVals_AB_Test, NN_AB, 'NN', dataATest, dataBTest);
    ErrorAnalysis.confusion(xVals_AB_Test, yVals_AB_Test, KNN_AB, 'KNN', dataATest, dataBTest);
    
    ErrorAnalysis.confusion(xVals_CDE, yVals_CDE, MED_CDE, 'MED', dataC, dataD, dataE);
    ErrorAnalysis.confusion(xVals_CDE, yVals_CDE, GED_CDE, 'GED', dataC, dataD, dataE);
    ErrorAnalysis.confusion(xVals_CDE, yVals_CDE, MAP_CDE, 'MAP', dataC, dataD, dataE);
    ErrorAnalysis.confusion(xVals_CDE_Test, yVals_CDE_Test, NN_CDE, 'NN', dataCTest, dataDTest, dataETest);
    ErrorAnalysis.confusion(xVals_CDE_Test, yVals_CDE_Test, KNN_CDE,'kNN', dataCTest, dataDTest, dataETest);

end

function gauss = gen_normal_distribution(N, mu, covar)
    rand = randn(N,2);
    R = chol(covar);
    gauss = repmat(transpose(mu), length(rand), 1) + rand*R;
end

function [V,D] = get_eigenvalues(matrix)
    [V,D] = eig(matrix);
end