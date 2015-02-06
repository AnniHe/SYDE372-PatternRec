%TEMP - not finished

classdef ClassA < handle
    properties
      N_A = 200
      muA = [5; 10];
      covarA = [8 0; 0 4];
      colourA = 'm';
      data = MathHelper.gen_normal_distribution(Na, muA, covarA);
    end
    
    properties
      
    end

    methods
        function obj = ClassA()
           data = MathHelper.gen_normal_distribution(N_A, muA, covarA); 
        end
    end

end