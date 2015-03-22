classdef MICD

    methods (Static)
        function class = classify(point, prototypes, weight)
            class = zeros(size(points,1),1);
            for i = 1:length(points)
                classDistances = [];
                for j = 1:length(prototypes)
                    m = prototypes(:,j) * weight;
                    distance = sqrt((m(1)-points(i,1))^2 + (m(2)-points(i,2))^2);
                    classDistances = [classDistances distance];
                end
                [~, classType] = min(classDistances);
                class(i) = classType;
            end
        end
    end
    
end

