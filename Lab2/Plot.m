classdef Plot

    methods (Static)
        function plotNormPdf(x, mu, sigma, colour)
            pdf = normpdf(x, mu, sigma);
            plot(x, pdf, colour, 'LineWidth',2);
        end
        function plotExponentialPdf(x, mu, colour, scalingfactor)
            pdf = exppdf(x,mu);
            plot(x, pdf, colour, 'LineWidth',2);
        end
        function plotUniformPdf(x, a, b, colour)
            pdf = unifpdf(x,a,b).*heaviside(x);
            plot(x, pdf, colour, 'LineWidth',2);
        end
        function applyCase(graph, varargin)
            switch graph
                case 1
                    legend('Estimated Density','True Density');
                    xlabel('x');
                    ylabel('PDF(X)');
                case 2
                    names = ['al'; 'bl'; 'cl'];
                    for i = length(varargin)
                        legend(varargin{i}, names(i));
                    end
                    xlabel('x');
                    ylabel('y');
            end

        end
    end
    
end

