function classType = MAP_Class(pt, meanArgs, covarArgs, probArgs)
    classVals = [];
    for i = 1:length(meanArgs)
        shit = i*2;
        disp(i);
        val = (probArgs(i)* exp( (-1/2) * transpose(pt - meanArgs(i))*inv(covarArgs(:,shit-1:shit))*(pt - meanArgs(i))))/((2*pi)*sqrt(det(covarArgs(:,shit-1:shit))));
        classVals = [classVals val];
    end

    [minVal, classType] = max(classVals);
end