function classType = MAP_Class(pt, meanArgs, covarArgs, probArgs)
    classVals = [];
    for i = 1:length(meanArgs)
        shit = i*2;
%         disp(i);
%         val = (probArgs(i)* exp( (-1/2) * transpose(pt - meanArgs(i))*inv(covarArgs(:,shit-1:shit))*(pt - meanArgs(i))))/((2*pi)*sqrt(det(covarArgs(:,shit-1:shit))));
        Sigma = covarArgs(:,shit-1:shit);
        Mu = meanArgs(:,i);
        Prob = probArgs(i);
        %Gause 2D equation * Posterior Probability
        val = 1/(sqrt(2*pi)*det(Sigma)) * exp(-0.5*(pt-Mu)'*inv(Sigma)*(pt-Mu)) * Prob;
        classVals = [classVals val];
    end

    [minVal, classType] = max(classVals(:));
end