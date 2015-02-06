function classType = MED_Class(point, meanArgs)
classDistances=[];

    %compare mean A and mean B
    for i = 1:length(meanArgs)
        m = meanArgs(:,i);
        distance = sqrt((m(1)-point(1))^2 + (m(2)-point(2))^2);
        classDistances =[classDistances distance];
    end

    [minDistance, classType] = min(classDistances);
end