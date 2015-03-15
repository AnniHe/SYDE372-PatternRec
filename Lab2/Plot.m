classdef Plot

    methods (Static)
        function plotNormPdf(x, mu, sigma, colour)
%             figure();
            pdf = normpdf(x, mu, sigma);
            plot(x, pdf, colour, 'LineWidth',2);
%             hold off;
        end
        function plotExponentialPdf(x, mu, colour, scalingfactor)
            pdf = exppdf(x,mu);
            plot(x, pdf, colour, 'LineWidth',2);
        end
        function plotUniformPdf(x, a, b, colour)
            pdf = unifpdf(x,a,b).*heaviside(x);
            pd2 = makedist('Uniform','lower',a,'upper',b);

            plot(x, pdf, colour, 'LineWidth',2);
            hold on
%             prev = [0:.1:min(x)];
%             after = [max(x):.1:10];
%             plot(prev, zeros(1, length(prev)));
%             plot(after, zeros(1, length(after)));
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

