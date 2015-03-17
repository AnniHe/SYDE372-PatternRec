%SEQUENTIAL CLASSIFIER

%let a and b represent data points in classes a and b. let j = 1
%randomly select one point from a and one point from b
%create a discriminant G using MED with the two points as prototypes
%using all of the data in a and b, work out the confusion matrix entries
%remove all the points from a and b that are classified correctly as B and A
%if a and b still contain points, repeat process
%once all points have been correctly, discriminant is good. 
% j = j+1

clear vars; clc; close all;
gridSize = 1;
xrange = 0:800;

load('lab2_3.mat');

% [xVals, yVals, grid] = Plotter.prepareGrid(gridSize, a,b);
ab = vertcat(a,b);

%add a class identifier for the class a and b. This is could be done
%better, but we only have two classes which works fine for now
a_true = [repmat(a,1) ones(size(a,1),1)];
twos = repmat(2,1,size(b,1))';
b_true = [repmat(b,1) twos];
ab_true = vertcat(a_true, b_true);

numDiscriminants = 1:1;
aveErrorRate = [];
minErrorRate = [];
maxErrorRate = [];
stddevErrorRate = [];

for i = numDiscriminants
    errorRate = [];
    
    for numClassifiers = 1:10
        figure();
        scatter(a(:,1), a(:,2), 5, 'b', 'filled');
        hold on;
        scatter(b(:,1), b(:,2), 5, 'g', 'filled');
        hold on;
        numPairs = 0;
        axis([50 550 0 450]);
        
        a_classify = repmat(a, 1);
        b_classify = repmat(b, 1);
        
        j = 1;
        naB = 1;
        nbA = 1;
        while j <= i || (naB ~= 0 && nbA ~= 0)
%         while naB ~= 0 || nbA ~= 0
            pointa = datasample(a_classify,1,1);
            pointb = datasample(b_classify,1,1);
            ab_classify = vertcat(a_classify, b_classify);
            MED_res = MED(ab_classify, [pointa' pointb']);      
            confMat = confusionMatrix(MED_res, a_classify, b_classify);
            newDiscrim = [ab_classify MED_res]; 
            MED_line = decisionLine(pointa, pointb, xrange);
            
            naB = confMat(3);
            nbA = confMat(2);
            
            %get error of the classifier at the moment        
            probError = classifierError(ab_true, newDiscrim);
            if naB == 0
                %get the correctly classified points. i.e. the intersection between
                %the new discriminant output and the true class data
                b_intersect = intersect(newDiscrim, b_true, 'rows');
                b_remove = b_intersect(:,[1,2]);  

                %subtract the correctly labelled samples from the "to-classify"
                %list. find indices of to-delete elements, remove them by indices
                [X, ind_b] = ismember(b_remove, b_classify, 'rows');  
                b_classify([ind_b(ind_b ~= 0)],:) = [];
                
                scatter(b_remove(:,1), b_remove(:,2), 5, 'm', 'filled');
            elseif nbA == 0
                a_intersect = unique(intersect(newDiscrim, a_true, 'rows'), 'rows');
                a_remove = a_intersect(:,[1,2]);
                [X, ind_a] = ismember(a_remove, a_classify, 'rows');
                a_classify([ind_a(ind_a ~= 0)],:) = [];
                scatter(a_remove(:,1), a_remove(:,2), 5, 'k', 'filled');
            end
            numPairs = numPairs + 1;
           
            %check whether we are done with this discriminant
            if naB == 0 || nbA ==0
                plot(xrange, MED_line, 'r'); hold on;
                j = j+1;
            end
            
        end
        if (j-1) == i
            errorRate = [errorRate probError];
        else
            numClassifiers = numClassifiers - 1;
        end
    end
    
    aveErrorRate = [aveErrorRate mean(errorRate)];
    minErrorRate = [minErrorRate min(errorRate)];
    maxErrorRate = [maxErrorRate max(errorRate)];
    stddevErrorRate = [stddevErrorRate std(errorRate)];
end

figure();
plot(numDiscriminants, aveErrorRate, 'b'); hold on;
% plot(numDiscriminants, aveErrorRate, 'om'); hold on;
plot(numDiscriminants, minErrorRate, 'r'); hold on;
% plot(numDiscriminants, minErrorRate, 'or'); hold on;
plot(numDiscriminants, maxErrorRate, 'g'); hold on;
% plot(numDiscriminants, maxErrorRate, 'og'); hold on;
plot(numDiscriminants, stddevErrorRate, 'k'); hold on;
% plot(numDiscriminants, stddevErrorRate, 'ok'); hold on;
legend('Average', 'Minimum','Maximum', 'Std.Dev.');

