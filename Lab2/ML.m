classdef ML
    %ML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function class = classify(x, y, varargin)
            class = -1;
            pdfs = [];
            pdfA = cell2mat(varargin(1));
            pdfB = cell2mat(varargin(2));
            pdfC = cell2mat(varargin(3));
            probA = pdfA(x, y);
            probB = pdfB(x, y);
            probC = pdfC(x, y);
%             for i = 1:length(varargin)
%                 pdfs = [pdfs varargin{1}];
%             end
            [~, class] = max([probA, probB, probC])
            
        end
    end
    
end

