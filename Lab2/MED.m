function classified = MED(points, prototypes)

classified = zeros(size(points,1),1);
    for i = 1:length(points)
        classDistances = [];
        for j = 1:length(prototypes)
            m = prototypes(:,j);
            distance = sqrt((m(1)-points(i,1))^2 + (m(2)-points(i,2))^2);
            classDistances =[classDistances distance];
        end
        [minDistance, classType] = min(classDistances);
        classified(i) = classType;
    end
end