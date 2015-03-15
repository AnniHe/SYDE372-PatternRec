classdef Plot
    %PLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function plotNormPdf(x, mu, sigma, colour)
%             figure();
            pdf = normpdf(x, mu, sigma);
            plot(x, pdf, colour);
%             hold off;
        end
        function plotExponential(x,lamda, colour)
            pdf = exppdf(x, lamda);
            plot(x, pdf, colour);   
        end
    end
    
end

