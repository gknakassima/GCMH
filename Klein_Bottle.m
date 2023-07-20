function Klein_Bottle
  
    First      = [-5 -5 -2 0 0];
    Last       = [5 5 2 2*pi 2*pi];
    Division   = [40 40 20 40 40]; 
    FirstPoint = [-3 0 0 pi pi];
    FileName   = 'Klein.pol';

    %MarchingSimplex(4, 2, First, Last, Division, @FKlein, FileName);
    %ContinuationSimplex(4, 2, First, Last, Division, FirstPoint, @FKlein, FileName);
    %GeneralizedMarchingHyperCube(5, 3, First, Last, Division, FileName, @FKlein);
    GeneralizedContinuationHyperCube(5, 3, First, Last, Division, FirstPoint, @FKlein, FileName);
    
    
    % Klein
    function [f] = FKlein(x) 
       a    = 3;
       f(1) = x(1) - (a+cos(x(4)/2)*sin(x(5))-sin(x(4)/2)*sin(2*x(5)))*cos(x(4));
       f(2) = x(2) - (a+cos(x(4)/2)*sin(x(5))-sin(x(4)/2)*sin(2*x(5)))*sin(x(4));
       f(3) = x(3) - sin(x(4)/2)*sin(x(5))+cos(x(4)/2)*sin(2*x(5));
       return
    end

    

 end 
