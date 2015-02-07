classdef ErrorAnalysis

    methods (Static)
        function confusion(xVals, yVals, z_classifier, type, varargin)

            dataset = [];
            
            %build args for vercat
             for i = 1:length(varargin)

                new = [varargin{i}, ones(length(varargin{i}), 1)*i];

                if length(dataset) == 0
                    dataset = new;
                else
                    dataset = vertcat(new, dataset);
                end

             end

            x = dataset(:,1);
            y = dataset(:,2);
            z = dataset(:,3);

            classifierResult = griddata(xVals,yVals, z_classifier, x,y, 'nearest');
            
            confusionmatrix = confusionmat(z, classifierResult);
            disp(type);
            disp(confusionmatrix);
            
            %total error points/total points
            Perror = 1 - trace(confusionmatrix)/sum(confusionmatrix(:));
            disp(Perror);
            
            
        end
    end
end