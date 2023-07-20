function Ambiguity_Manifold

    ndiv       = 5;    
    First      = [-1   -1   -1   -1];
    Last       = [ 1    1    1    1];
    Division   = ndiv*ones(1,4); 
    FirstPoint = [0.1   0.1   0.1   0.1];
    FileName   = ['AmbiguityManifold_div' num2str(ndiv) '.pol'];
    
    
    %MarchingSimplex(6, 3, First, Last, Division, @amb_manifold, FileName);
    %ContinuationSimplex(6, 3, First, Last, Division, FirstPoint, @amb_manifold, FileName);
    %GeneralizedMarchingHyperCube(6, 3, First, Last, Division, FileName, @amb_manifold); ;
    GeneralizedContinuationHyperCube_times(4, 1, First, Last, Division, FirstPoint, @amb_manifold, FileName);

    return 

    function [f] = amb_manifold(x)
        f = x(1)*x(2)*x(3)*x(4) - 0.0001;
        return
    end

end 
