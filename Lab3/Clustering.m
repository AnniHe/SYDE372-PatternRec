classdef Clustering
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        K = 10;
    end
    
    methods (Static)
        function doKMeans(data)
            
            nrSamples = length(data);
            %Cluster prototype means
            prototypes = pickRandomPrototypes(data)
            
            %Calculate the distances
            classified = classifyMED(data, prototypes);
            
            %iterate over and over again
        end
    end
    
end

function CP = pickRandomPrototypes(data)
    
    % gets a row of 10 columns with random indices from data set range
    perms = randperm(size(data,2), Clustering.K);
    %Cluster Prototypes
    CP = zeros(length(perms),2);
    
    for i = 1:length(CP)
        CP(i,:) = [data(1,perms(i)) data(2, perms(i))];
    end
end

function classified = classifyMED(data, CP)
    % Let's make it less intuitive for the students YAY
    % Use transpose to have x's in col1, y's in col 2
    data = data.'
    clz = length(CP);
    distances = zeros(clz, 1);
    
    points = length(data)
    classified = zeros(points, 1);
    
    %Each sample Point (xi,yi)
    for i = 1:points
        %Each cluster prototype
        for j = 1:clz
            distances(j) = sqrt((data(i, 1) - CP(j, 1))^2 + (data(i, 2) - CP(j, 2))^2);
        end
        [~, index] = min(distances);
        classified(i) = index;
    end
end

