classdef Plotter
   methods(Static)
       
       function [xVals, yVals, classGrid] = prepareGrid (gridSize, varargin)
            offset = 2;

            % Input gridsize and the (2D) data sets
            minX_val=[];
            minY_val=[];
            maxX_val=[];
            maxY_val=[];

            %iterate through each dataset and get the bounds
            for i = 1:length(varargin)
                xVals = varargin{i}(:, 1);
                yVals = varargin{i}(:, 2);
                minX_val = [minX_val min(xVals)];
                maxX_val = [minX_val max(xVals)];
                minY_val = [minY_val min(yVals)];
                maxY_val = [maxY_val max(yVals)];
            end

            %return a nice grid
            xVals = min(minX_val) - offset: gridSize : max(maxX_val) + offset;
            yVals = min(minY_val) - offset: gridSize : max(maxY_val) + offset;

            classGrid = zeros(length(yVals), length(xVals));
        end
      
        function classifier = plot_MED(grid, xVals, yVals, meanArgs, fig)
        classifier = repmat(grid, 1);
            for i = 1:numel(classifier)
                [row col] = ind2sub(size(classifier), i);
                classIndex = MED_Class([xVals(col);yVals(row)], meanArgs);
                classifier(i) = classIndex;
            end
            figure(fig)
            contour(xVals, yVals, classifier);
            hold on;
        end

        function classifier = plot_GED(grid, xVals, yVals, mean_transformed_Args, mat_transforms)
            classifier = repmat(grid, 1);
            for i = 1:numel(classifier)
                [row col] = ind2sub(size(classifier), i);
                classIndex = GED_Class([xVals(col);yVals(row)], mean_transformed_Args, mat_transforms);
                classifier(i) = classIndex;       
            end
%             figure(f);
            contour(xVals, yVals,classifier);
            hold on;
        end

        function classifier = plot_MAP(grid, xVals, yVals, meanArgs, covarArgs, probsArgs)
            classifier = repmat(grid, 1);
            for i = 1:numel(classifier)
                [row col] = ind2sub(size(classifier), i);
                classIndex = MAP_Class([xVals(col);yVals(row)],meanArgs, covarArgs, probsArgs);
                classifier(i) = classIndex;
            end
        end
        
        function plot_nn(k, grid, xVals, yVals, fig, varargin)

            classifier = repmat(grid, 1);
            for i = 1:numel(classifier)
                [row col] = ind2sub(size(classifier), i);
                classIndex = NN_Class(k,[xVals(col);yVals(row)], fig,varargin{:});
                classifier(i) = classIndex; 
            end
            figure(fig);
            contour(xVals,yVals,classifier);
            hold on;
        end
        
    end
end