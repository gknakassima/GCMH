function Complex_Cosine
  
    First      = [-pi/2 -pi/2 -pi/2 -pi/2];
    Last       = [pi/2 pi/2 pi/2 pi/2];
    Division   = [5 5 5 5]; 
    FirstPoint = [0.0 0.0 1.0 0.0];
    FileName   = 'Cosine.pol';

    %MarchingSimplex(4, 2, First, Last, Division, @zcosw, FileName);
    %ContinuationSimplex(4, 2, First, Last, Division, FirstPoint, @zcosw, FileName);
    %GeneralizedMarchingHyperCube(4, 2, First, Last, Division, FileName, @zcosw);
    %GeneralizedContinuationHyperCube(4, 2, First, Last, Division, FirstPoint, @zcosw, FileName);
    
    
    return 

    % z^2 = cos(w)
    function [f] = z2cosw(x) 
       f(1) = x(3)*x(3) - x(4)*x(4) - 0.5*exp(-x(2))*cos(x(1)) - 0.5*exp(x(2))*cos(-x(1));
       f(2) = 2.0*x(3)*x(4) - 0.5*exp(-x(2))*sin(x(1)) - 0.5*exp(x(2))*sin(-x(1));
       return
    end
    
    % z = cos(w)
    function [f] = zcosw(x) 
       f(1) = x(3) - 0.5*exp(-x(2))*cos(x(1)) - 0.5*exp(x(2))*cos(-x(1));
       f(2) = x(4) - 0.5*exp(-x(2))*sin(x(1)) - 0.5*exp(x(2))*sin(-x(1));
       return
    end

    % z^2 = w^5 + (t + is) cond inic z = (1,0) w = (0,-1) */
    function [f] = z2w5(x) 
       t    = 1.0;
       s    = 1.0;
       f(1) = x(1)^2-x(2)^2 - x(3)^5 - 10.0*x(3)^3*x(4)^2 - 5*x(3)*x(4)^4 - t;
       f(2) = 2*x(1)*x(2) - 5.0*x(3)^4*x(4) - 10.0*x(3)^2*x(4)^3 - x(4)^5 - s;
       return
    end

 end 
