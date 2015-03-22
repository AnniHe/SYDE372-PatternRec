classdef Clustering
    %Explanations
    % classified: 160 elements (for each sample point), value corresponds
    % to a number 1-10 according to which class prototype is the closest.
    % p0 - initial prototype selection - this is random
    % pn - converted prototypes - after many iterations
    
    properties (Constant)
        K = 10;
        p0;
    end
    
    methods (Static)
        function doKMeans(data, p0)
            
            nrSamples = length(data);
            %Cluster prototype means
%             Clustering.p0 = pickRandomPrototypes(data);
%             p0 = Clustering.p0;
            
            %Calculate the distances
            classified = classifyMED(data, p0);
            prevClassified = [];
            
            %iterate over and over again
            %Now iteratively go create new mean class prototypes
            while(~isequal(prevClassified, classified))
                %compute new cluster prototype means
                prevClassified = classified;
                pn = pickMeanPrototypes(data,classified)
                classified = classifyMED(data, pn);
                
                
                %compare is classifed the same as before or not?
                %if yes, retry with new calulated means
                %if no, we are done
                
            end
            
            Plotter.plotKMeans(data, p0, pn)
            
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
        
        function dofuzzyKMeans(data, p0)
        end
    end
    
end



function CP = pickMeanPrototypes(data, classified)
    %classified - count all the instances of 1's 2's and their index in the
    %-160 range
    CP = zeros(Clustering.K, 2);
    %Sum the number of occurances of an
    %Each Cluster type 1-10
    for i = 1:Clustering.K
        shizzles = find(classified == i);
        x = data(1, shizzles);
        y = data(2, shizzles);
        CP(i,:) = [mean(x) mean(y)];

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



