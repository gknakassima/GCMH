function Torus_R10
  
    First      = [-1.1 -1.1 -1.1 -1.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1];
    Last       = [1.1 1.1 1.1 1.1 0.1 0.1 0.1 0.1 0.1 0.1];
    Division   = [10 10 10 10 1 1 1 1 1 1]; 
    FirstPoint = [1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
    FileName   = 'Torus_R10.pol';
  
    %MarchingSimplex(6, 3, First, Last, Division, @s1xs1, FileName);
    %ContinuationSimplex(6, 3, First, Last, Division, FirstPoint, @s1xs1, FileName);
    %GeneralizedMarchingHyperCube(6, 3, First, Last, Division, FileName, @s1xs1); 
    GeneralizedContinuationHyperCube(10, 8, First, Last, Division, FirstPoint, @s1xs1, FileName);
    
    
    return 

    % S1 x S1 x S1
    function [f] = s1xs1(x) 
       f(1) = x(1)*x(1) + x(2)*x(2) - 1.0;
       f(2) = x(3)*x(3) + x(4)*x(4) - 1.0;
       f(3) = x(5);
       f(4) = x(6); 
       f(5) = x(7);
       f(6) = x(8);
       f(7) = x(9);
       f(8) = x(10);
       return
    end
    
 end 
