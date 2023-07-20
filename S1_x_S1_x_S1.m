function S1_x_S1_x_S1
  
    First      = [-0.6 -0.6 -1.1 -1.1 -4.1 -4.1];
    Last       = [0.6 0.6 1.1 1.1 4.1 4.1];
    Division   = [10 10 10 10 10 10]; 
    FirstPoint = [0.5 0.0 1.0 0.0 2.0 0.0];
    FileName   = 'S1xS1xS1.txt';

    %MarchingSimplex(6, 3, First, Last, Division, @s1xs1xs1, FileName);
    %ContinuationSimplex(6, 3, First, Last, Division, FirstPoint, @s1xs1xs1, FileName);
    %GeneralizedMarchingHyperCube(6, 3, First, Last, Division, FileName, @s1xs1xs1); 
    GeneralizedContinuationHyperCube(6, 3, First, Last, Division, FirstPoint, @s1xs1xs1, FileName);
    
    
    return 

    % S1 x S1 x S1
    function [f] = s1xs1xs1(x) 
       f(1) = x(1)*x(1) + x(2)*x(2) - 0.25;
       f(2) = x(3)*x(3) + x(4)*x(4) - 1.0;
       f(3) = x(5)*x(5) + x(6)*x(6) - 4.0;
       return
    end
    
 end 
