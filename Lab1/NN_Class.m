function classType = NN_Class(k, point, varargin)

    minDistances = [];
    %for A,B,C,D,E
    for i = 1:length(varargin)
        classDistances = [];
        %each row
        for j = 1:length(varargin{i})
            %[row col] = ind2sub(size(varargin{i}), j);
            distance = sqrt((varargin{i}(j,1) - point(1))^2 + (varargin{i}(j,2) - point(2)).^2);
            classDistances =[classDistances distance];
        end
        
        if k > 1
            tmp = sort(classDistances);
            sorted = tmp(1:k);
            avgDistance = mean(sorted);
            minDistances = [minDistances avgDistance];
        else
            closest = min(classDistances);
            minDistances = [minDistances closest];
            
        end
            
    end
    %classType index should be 1 or 2
    [dist, classType] = min(minDistances);

    
end