classdef Plot
    %PLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function plotNormPdf(mu, x, sigma, colour)
            figure();
            pdf = normpdf(x, mu, sigma);
            plot(x, pdf, colour);
%             hold off;
        end
    end
    
end

