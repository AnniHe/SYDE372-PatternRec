classdef Gaussian
    %GAUSSIAN Summary of this class goes here
    %   Detailed explanation goes here

    methods (Static)
        
        function p2d = plot2D(al,bl,cl)
            
           gridSize = 2;

            %Create a grid
            xGridDim = max([max(al(:,1)),max(bl(:,1)),max(cl(:,1))]) - min([min(al(:,1)),min(bl(:,1)),min(cl(:,1))]);
            yGridDim = max([max(al(:,2)),max(bl(:,2)),max(cl(:,2))]) - min([min(al(:,2)),min(bl(:,2)),min(cl(:,2))]);
            grid = zeros(xGridDim,yGridDim);
            
            pdfA = probability(xGridDim, yGridDim, al, grid);
            pdfB = probability(xGridDim, yGridDim, bl, grid);
            pdfC = probability(xGridDim, yGridDim, cl, grid);
            
            %Classify each square on grid using ML
            for i = 1:xGridDim-1
                for j = 1:yGridDim-1
                    grid(i,j) = ML.classify(i,j, pdfA, pdfB, pdfC);
                end
            end
            
            %Plot points and ML boundary
            figure();
            hold on;
            contour(grid, 'LineWidth',2);
            s1 = scatter(al(:,1),al(:,2));
            s2 = scatter(bl(:,1),bl(:,2));
            s3 = scatter(cl(:,1),cl(:,2)); 
            Plot.applyCase(2, s1, s2, s3);
            hold off;
        end
        
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
   
        function result = evalGuassian(x, mean, sigma, n)
            if n == 2
                sigmaMagnitude = det(sigma)^(.5);
                result = (1/(sigmaMagnitude*sqrt(2*pi)))*exp(-1/2*(inv(sigma)*(x-mean).^2));
            else
                result = (1/(sigma*sqrt(2*pi)))*exp(-1/(2*(sigma^2))*((x-mean).^2));
            end
        end
        
    end
    
end

function pdf = probability(x, y, data, grid)
            pdf = grid;
            for i = 1:x
                for j = 1:y
                    pdf(i,j) =  1/(sqrt(2*pi)*det(cov(data))) * exp(-0.5*([i j]-mean(data))*inv(cov(data))*([i j]-mean(data))') * 0.5;
                end
            end
        end

