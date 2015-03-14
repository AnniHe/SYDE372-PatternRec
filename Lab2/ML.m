classdef ML
    methods (Static)
        %return 1,2,3 based on pdf that is most likely
        %x,y are points on the grid
        function class = classify(x, y, varargin)
            probs = [];
            for i = 1:length(varargin)
                probs = [probs varargin{i}(x,y)];
            end
            [~,class] = max(probs);
        end
    end
end

