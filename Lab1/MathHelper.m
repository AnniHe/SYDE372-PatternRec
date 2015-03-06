classdef MathHelper
    %MATHHELPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function gauss = gen_normal_distribution(N, mu, covar)
            rand = randn(N,2);
            R = chol(covar);
            gauss = repmat(transpose(mu), length(rand), 1) + rand*R;
        end

        function [V,D] = get_eigenvalues(matrix)
            [V,D] = eig(matrix);
        end
    end
    
end

