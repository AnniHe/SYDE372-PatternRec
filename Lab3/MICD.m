classdef MICD

    methods (Static)
        function class = classify(points, prototypes, weights)
            class = zeros(size(points,1),1);
            for i = 1:size(points,1)
                classDistances = [];
                for j = 1:size(prototypes,2)
                    j_weight = j*2;
                    point_t = weights(:,j_weight-1:j_weight) * points(i,:)';
                    prototype_t = weights(:,j_weight-1:j_weight) * prototypes(:,j);
                    
                    distance = sqrt((prototype_t(1)-point_t(1))^2 + (prototype_t(2)-point_t(2))^2);
                    classDistances = [classDistances distance];
                end
                [~, classType] = min(classDistances);
                class(i) = classType;
            end
        end
    end
    
end