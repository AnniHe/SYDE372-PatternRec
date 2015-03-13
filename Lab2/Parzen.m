classdef Parzen
    %   Non-parametric estimation using a gaussian parzen function
    %   sample can be 1 dimensional or a 2dimensional dataset    
    methods (Static)
        function p = plot(sample, sigma)

            if size(sample, 2) == 2
                n = 2;
                sigma = cov(sample);
                range = [min(sample) - 1 : .1 ; max(sample) +1];
                
            else
                n = 1;
            end

            N = length(sample); %number of samples
            parzen = [];
            K = 100 %optimize this parameter
            h = K/sqrt(N); %width of window, formula in notes

            % height = 1/(N) * 1/width * 1;
            range = [0:.1:10];

            for x = range
                sum = 0;
                %sum the sample constributions at each discrete step
                for i = 1:N
                    mean = sample(i,:);  
                    sum = sum + 1/h* evalGuassian(x, mean', sigma, n);
                end

                %keep track of the y (1D) z (2D) vals at each discrete step
                parzen = [parzen 1/N*sum];
            end
            figure();
            plot(range,parzen);
            
        end
    end
        
        
end
    
    
function result = evalGuassian(x, mean, sigma, n)
    if n == 2
        sigmaMagnitude = det(sigma)^(.5);
        result = (1/(sigmaMagnitude*sqrt(2*pi)))*exp(-1/2*(inv(sigma)*(x-mean).^2));
    else
        result = (1/(sigma*sqrt(2*pi)))*exp(-1/(2*(sigma^2))*((x-mean)^2));
    end
    
end
        

