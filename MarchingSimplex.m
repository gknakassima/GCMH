function MarchingSimplex(n, k, First, Last, Division, Func, filename) 
    file  = fopen(filename,'w');
    fprintf(file,'%3d %3d\n',n,k);
    fprintf(file,'%3d ',Division(1:n));
    fprintf(file,'\n\n');

    Delta = (Last-First)./Division;                 % Hypercube sizes
    Pert  = random('Normal',0,1,2^n,n)*0.000001;    % Perturbation to vertices to avoid problems at vertices
    Base  = InitGrid(n, Division);                  % Initialize grid indices

    % Initialize grid for simplex vertices inside cubes
    for i = 1:n
        DivisionP(i) = i;
    end
    [BaseP] = InitGrid(n,DivisionP);

    % Initialize grid for simplex edges inside cubes
    for i = 1:k+1
        DivisionC(i) = n-k+i;
    end
    [BaseC] = InitGrid(k+1,DivisionC);

    Ncube    = 0;
    Nsimplex = 0;

    % Go through each hypercube in the grid
    for g = 0:Base(n+1)-1 

      %Ncube  = Ncube + 1; 
      
      % Find hypercube indices (ix1,...,ixn) in the grid
      [Grid] = GridCoords(n,g,Base);
      % Calculates (x1,...,xn) of perturbed vertices and respective value of F
      [Vert, FVert] = GenVert(n, Grid, Pert, First, Delta, Func);

      % Loop through simplices in the hypercube g
      for s = 0:BaseP(n+1)-1

        %Nsimplex = Nsimplex + 1;  
        %fprintf('Simplexo: Ncube = %d  Nsimplex = %d  g = %d  s = %d\n',Ncube,Nsimplex,g,s);
        fprintf('Simplexo: g = %d  s = %d\n',g,s);

        % Indices of local hypercube vertices in the simplex (in permutation notation from the book)
        [f] = GridCoords(n,s,BaseP);
        f = f+1;
        % Map permutation notation back to usual number list
        [P] = Map_Perm(n,f);
        % Local hypercube indices of vertices in the simplex
        [Simp] = GenLabelSimplex(n,P);

        VertexManifold = [];
        FaceVertex     = [];
        NumberVertex   = 0;
        
        % Loop through all possible simplex (k+1)-faces in the hypercube
        for j = 0:BaseC(k+2)-1
            
            % Make a guess at possible vertex combination that may make a (k+1)-face
            [C]   = GridCoords(k+1,j,BaseC);
            C     = C+1;
            % Check if vertex combination is in lexicographic order; if not, discard the combination
            [lex] = Lexico(k+1,C);
            
            % If it is in lexicographic order, check for manifold
            if lex == 1
                % Map hypercube vertices to (k+1)-face
                Face = zeros(1,k+1);
                for i = 1:k+1
                    Face(i) = Simp(C(i));
                end
                % Check if (k+1)-face is transversal to manifold
                % Vertex = approx. coordinates of intersection point between face and manifold
                [Vertex, trans] = GetVertexManifold(n,k,Face,Vert,FVert);
                % If (k+1)-face is transversal, add to list
                if trans == 1
                    VertexManifold = [VertexManifold; Vertex];
                    FaceVertex     = [FaceVertex; Face];
                    NumberVertex = NumberVertex + 1;
                end
            end

        end

        if NumberVertex > 0

           [nSkel, Skel, AdjSkel] = Skeleton(n,Simp,k,NumberVertex,FaceVertex);  

           fprintf(file,'%3d ',g);
           fprintf(file,'%3d ',Grid);
           fprintf(file,'\n');
           fprintf(file,'  1\n');

           fprintf(file,'%3d \n',nSkel(1));
           for i = 1:nSkel(1)
              fprintf(file,'%3d ',Skel{1}(i,:));
              fprintf(file,'%15.8f ',VertexManifold(i,:));
              fprintf(file,'\n');
           end
           for j = 1:n-k
              fprintf(file,'%3d\n',nSkel(j+1));
              for i = 1:nSkel(j+1)
                 fprintf(file,'%3d ',AdjSkel{j}{i}(:));
                 fprintf(file,'\n');
              end
           end
           fprintf(file,'\n');

         end

      end

    end 
    
    % Print end of file
    fprintf(file,'-1\n');
    fclose(file);

    return
end

%% InitGrid: Initializes vector of sequential indices
function [Base] = InitGrid(n, Division) 
   Base(1) = 1;
   for i = 2:n+1
      Base(i) = Base(i-1)*Division(i-1);
   end
   return
end

%% GridCoords: Find n-dimensional indices (i1,...,in)
function [Grid] = GridCoords(n,i,Base) 
   copy = i;
   for j = n:-1:2
      aux     = mod(copy,Base(j));
      Grid(j) = (copy-aux)/Base(j);
      copy    = aux;
   end
   Grid(1) = copy;
   return
end  

%% GenVert: Calculate perturbed vertices' coordinates and value of F
function [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func) 
   Vert  = [];
   FVert = [];
   for i = 0:2^n-1 
      [Coords]    = HyperCubeCoords(n,i);
      [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
      [VHC]       = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
      Vert        = [Vert; VHC];
      [FVHC]      = Func(n,VHC); 
      FVert       = [FVert; FVHC];
   end
   return
end

%% HyperCubeCoords: Find n-dimensional indices of a vertex in a hypercube
function [Coords] = HyperCubeCoords(n,i) 
   copy = i;
   for j = 1:n-1
      Coords(j) = mod(copy,2);
      copy      = (copy-Coords(j))/2;
   end
   Coords(n) = copy;
   return
end

%% HyperCubePert: Apply perturbation to hypercube vertices
function [CoordPert] = HyperCubePert(n,Grid,Coords,Pert)
   pot = 1;
   label = 0;
   for i = 1:n
      p     = abs(Coords(i) - mod(Grid(i),2));
      label = label + pot*p;
      pot   = 2*pot;    
   end
   CoordPert = Pert(label+1,:);
   return
end

%% HyperCube: Find cartesian coordinates of the vertices of the hypercube
function [VHC] = HyperCube(n,First,Delta,I,Coords,Pert) 
   for i = 1:n
      VHC(i) = First(i) + (I(i)+Coords(i))*Delta(i) + Pert(i);
   end
   return
end

%% GetVertexManifold: Find intersection between face and manifold
function [Vertex, trans] = GetVertexManifold(n,k,Face,Vert,FVert)
    % Creates linear system to find approximate intersection point
    for i = 1:k+1
      A(1,i) = 1;
      A(2:k+1,i) = FVert(Face(i)+1,:);
   end
   b = [1; zeros(k,1)];
   lamb = A\b;
   trans = 0;
   Vertex = zeros(1,n);
   % If all lambda > 0, find intersection point as linear combination of vertices
   if lamb >= 0
       for i = 1:k+1
          Vertex = Vertex + lamb(i)*Vert(Face(i)+1,:);
       end
       trans = 1;
   end   
   return
end

%% GenVertHyperCube
function [V FV] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label,Func)     
   [Coords]    = HyperCubeCoords(n,Label);
   [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
   [V]         = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
   [FV]        = Func(n,V); 
   return
end

%% Map_Perm: Turns permutation notation from the book back to usual notation
function [F] = Map_Perm(n,f)
    for i = 1:n
        F(i) = 0;
    end
    for i = 1:n
        if F(f(i)) == 0
            F(f(i)) = i;
        else
            for j = n:-1:f(i)+1
                F(j) = F(j-1);
            end
            F(f(i)) = i;
        end
    end
    return
end

%% Lexico: Checks if permutation is in lexicoghaphic (ascendent) order
function [lex] = Lexico(k,f)
    for i = 1:k-1
        if f(i) >= f(i+1)
            lex = 0;
            return
        end
    end
    lex = 1;
    return
end

%% GenLabelSimplex
function [Simp] = GenLabelSimplex(n,P)
    Simp(1) = 0;
    for i = 1:n
        Simp(i+1) = Simp(i)+2^(P(i)-1);
    end
end

%% Skeleton: Generates the combinatorial skeleton
function [nSkel,Skel,AdjSkel] = Skeleton(n, Simp, k, nVertex, FaceVertex)   
    nSkel(1) = nVertex;
    Skel{1}  = FaceVertex;
    % Loop through dimensions
    for s = k+1:n
        
        % Generates labels for the s-faces of the approximation
        for i = 1:s+1
            DivisionC(i) = n-s+i;
        end
        [BaseC]      = InitGrid(s+1,DivisionC);
        FacesSimp    = [];
        FacesAdj     = [];
        nSkel(s-k+1) = 0;
        
        % Loops through all possible vertex combinations
        for j = 0:BaseC(s+2)-1
            [C]   = GridCoords(s+1,j,BaseC);
            C     = C+1;
            [lex] = Lexico(s+1,C);
            
            if lex == 1
                % Generates face label using vertices
                for i = 1:s+1
                    Face(i) = Simp(C(i));
                end
                cont = 0;
                Adj  = [];
                
                % Loops through all faces in combinatorial skeleton
                for i = 1:nSkel(s-k)
                    % Checks if s-face contains (s-1)-face in the skeleton
                    [np] = Include(s,Skel{s-k}(i,:),s+1,Face);
                    % If positive, update adjacencies
                    if np == s
                        cont      = cont+1;
                        Adj(cont) = i;
                    end
                end
                
                % If current s-face contains any (s-1)-face
                if cont > 0
                    % Update number of (s-k)-faces in the skeleton
                    nSkel(s-k+1)           = nSkel(s-k+1)+1;
                    % Include vertex labels in the s-face
                    FacesSimp              = [FacesSimp; Face];
                    % Includes (s-1)-faces adjacent to s-face
                    FacesAdj{nSkel(s-k+1)} = Adj;                    
                end
            end
        end
        Skel{s-k+1}  = FacesSimp;
        AdjSkel{s-k} = FacesAdj;
    end
    return
end

%% Include: Verifies whether (m-1)-face y contains the (n-1)-face x
function [k] = Include(n,x,m,y)
   k = 0;
   for i = 1:n
       j = i;
       while (j < m) && (y(j) < x(i))
           j = j+1;
       end
       if y(j) == x(i)
           k = k+1;
       end
   end
   return
end
