function confM = confusionMatrix(classifierRes, varargin)
    dataset = [];
    for i = 1:length(varargin)
        new = [varargin{i}, ones(size(varargin{i},1), 1)*i];

        if length(dataset) == 0
            dataset = new;
        else
            dataset = vertcat(dataset,new);
        end
    end
    z = dataset(:,3);
    confM = confusionmat(z, classifierRes);
end