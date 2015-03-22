classdef Plotter
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function plotKMeans(data, p0, pn)
            figure;
            aplot(data);
            hold on;
            %Initial Prototypes
            plot(p0(:,1), p0(:,2), 'd', 'MarkerFaceColor','red', 'MarkerSize',6);
            %Convered Prototypes
            plot(pn(:,1), pn(:,2), 'd', 'MarkerFaceColor','blue', 'MarkerSize',6);
            legend('Inital Prototype', 'Converged Prototype',2);
        end
    end
    
end

