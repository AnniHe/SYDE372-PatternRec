classdef OneDimData
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        mu_true
        sigma_true
        lambda_true
        mu
        sigma
        lambda
        range
        pdf
        pdf_true
    end
    
    methods
        function OneDimData = OneDimData(data, m, v, l)
            OneDimData.data = data;
            OneDimData.mu_true = m;
            OneDimData.sigma_true = v;
            OneDimData.lambda_true = l;
            OneDimData.mu = mean(data);
            OneDimData.sigma = var(data);
            OneDimData.lambda = inv(OneDimData.mu);
            OneDimData.range = [min(data):0.01:max(data)];
            OneDimData.getNormPdf();
        end
        function plotAny(OneD)
            figure;
            plot(OneD.range, OneD.pdf, 'b');
            hold on;
            plot(OneD.range, OneD.pdf_true, 'r');
            title('Gaussian Parametric Estimation: Data');
            legend('Sample Distribution', 'True Distribution');
        end
        function plotGaussian(OneD)
            figure();
            Plot.plotNormPdf(OneD.range, OneD.mu_true, OneD.sigma_true, 'r');
            hold on;
            Plot.plotNormPdf(OneD.range, OneD.mu, OneD.sigma, 'b');
            hold off;
        end
        function plotExponential(OneD)
            figure();
            Plot.plotExponential(OneD.range, OneD.lambda_true, 'r');
            hold on;
            Plot.plotExponential(OneD.range, OneD.lambda, 'b');
            hold off;
        end
        function plotUniform(OneD)
%             figure();
%             Plot.plotExponential(OneD.range, OneD.lambda_true, 'r');
%             hold on;
%             Plot.plotExponential(OneD.range, OneD.lambda, 'b');
%             hold off;
        end
        function OneD = getNormPdf(OneD)
            OneD.pdf = normpdf(OneD.range, OneD.mu, OneD.sigma);
            OneD.pdf_true = normpdf(OneD.range, OneD.mu_true, OneD.sigma_true);          
        end
        function OneD = getExpPdf(OneD)
            OneD.pdf = exppdf(OneD.range, OneD.mu);
            OneD.pdf_true = exppdf(OneD.range, OneD.lambda_true);          
        end
        function OneD = getUniPdf(OneD)
           OneD.pdf = unifpdf(OneD.range, min(OneD.data), max(OneD.data)); 
           OneD.pdf_true = unifpdf(OneD.range, min(OneD.data), max(OneD.data));
        end
 
    end
    
end

