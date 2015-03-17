function error = classifierError(trueSet, classifyingSet)

%classified results
classifiedRes = repmat(trueSet, 1);

for index = 1:size(classifyingSet, 1)
    pointRow = classifyingSet(index, [1,2]);
    [X, ind_c] = ismember(pointRow, classifiedRes(:,[1,2]), 'rows');
    classifiedRes(ind_c, 3) = classifyingSet(index, 3);
end
confMat = confusionmat(trueSet(:,3),classifiedRes(:,3));
error = 1 - trace(confMat)/sum(confMat(:));

end