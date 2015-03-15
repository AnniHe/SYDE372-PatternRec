classdef Gauss
    
    methods (Static)
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

