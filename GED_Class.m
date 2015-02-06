function classType = GED_Class(point, meanArgs, transformArgs)
    classDistances=[];

    for i = 1:length(meanArgs)
        shit = i*2;
        pt_transformed = transformArgs(:, shit -1 : shit) * point;
        m = meanArgs(:,i);
        distance = sqrt((m(1)-pt_transformed(1))^2 + (m(2)-pt_transformed(2))^2);
        classDistances =[classDistances distance];
    end

    [minDistance, classType] = min(classDistances);
end