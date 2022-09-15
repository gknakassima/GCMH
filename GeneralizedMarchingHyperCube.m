function GeneralizedMarchingHyperCube(n, k, First, Last, Division, filename, Func)
    
    file = fopen(filename,'w');
    file2 = fopen('FullCoef.txt','w');

    fprintf(file,'%3d %3d\n',n,k);
    fprintf(file,'%3d ',Division(1:n));
    fprintf(file,'\n\n');

    Delta = (Last-First)./Division;              % Hypercube sizes
    Pert  = random('Normal',0,1,2^n,n)*1e-8; % Perturbation to vertices to avoid problems
    Base  = InitGrid(n, Division);
    
    fprintf('Cubes processed: \n');

    % Loops through all hypercubes
    for g = 0:Base(n+1)-1
        
        % Find hypercube indices (ix1,...,ixn) in the grid
        [Grid] = GridCoords(n,g,Base);
        
        fprintf('g = %3d\n',g);
        
        % Generates vertex coordinates and function values on them
        [Vert,FVert] = GenVert(n, Grid, Pert, First, Delta, Func);
    
        % Split hypercube faces into simplices
        [CubeFace,nCubeFace,FaceVert,CubeSimp,nCubeSimp] = SplitFaceSimplex(n, k);
        
        VertexManifold = [];
        FaceCubeVertex = [];
        SimpCubeVertex = [];
        CubeFaceLabel = [];
        NumberVertex   = 0;
        FullCoefList = [];
    
        % Loop through hypercube faces 
        for f = 1:nCubeFace
    
            % Loop through all hypercube face simplices
            for s = 1:nCubeSimp

                % Get hypercube vertices for the current simplex
                FaceIndex = CubeFace(f,:);
                FaceVertex = FaceVert(f,:);
                SimpFace = CubeSimp{f}(s,:);
    
                % Check if k-face is transversal to manifold
                % Vertex = approx. coordinates of intersection point between face and manifold
                [Vertex, trans] = GetVertexManifold(n,k,SimpFace,Vert,FVert);
    
                % If k-face is transversal, add to list
                if trans == 1
                    % Apply Newton's Method for checking
%                     [FullCoefVertex] = PseudoNewton(40,Vertex);
%                     [FullCoefVertex] = NewtonMethod(40,Vertex,FaceIndex);
%                     fprintf(file2,'%15.8f ',FullCoefVertex);
%                     fprintf(file2,'\n');
%                     [FullCoefList] = [FullCoefList; FullCoefVertex];

                    VertexManifold = [VertexManifold; Vertex];
                    FaceCubeVertex = [FaceCubeVertex; FaceVertex];
                    SimpCubeVertex = [SimpCubeVertex; SimpFace];
                    CubeFaceLabel  = [CubeFaceLabel; FaceIndex];
                    NumberVertex   = NumberVertex + 1;
                end

            end
        end
        
        % Connect manifold vertices into manifold edges
        [VertexEdgeManifold1,FaceEdgeManifold,NumberEdge,NumberFace] = GetEdgeManifold(n,k,NumberVertex,CubeFaceLabel);

        % Loop through (k+1)-faces with manifold vertices for ambiguities
        VertexEdgeManifold = [];
        NumberEdgeManifold = 0;
        for j = 1:NumberFace
            % If a (k+1)-face has ambiguities
            if NumberEdge(j) > 1
                [VertexEdgeManifold2,NumberEdgeManifold2] = SolveAmbiguity(n,k,Vert,FVert,FaceCubeVertex,NumberVertex,FaceEdgeManifold(j,:));
                [VertexEdgeManifold] = [VertexEdgeManifold; VertexEdgeManifold2];
                NumberEdgeManifold = NumberEdgeManifold + NumberEdgeManifold2;
            else
                [VertexEdgeManifold] = [VertexEdgeManifold; VertexEdgeManifold1(j,1:2)];
                NumberEdgeManifold = NumberEdgeManifold + 1;
            end
        end
        
        % Constructs the combinatorial skeleton
        if NumberVertex > 0
            fprintf(file,'%3d ',g);
            fprintf(file,'%3d ',Grid);
            fprintf(file,'\n');
            SkeletonHyperCube(file,n,k,NumberVertex,NumberEdgeManifold,VertexManifold,CubeFaceLabel,FaceCubeVertex,VertexEdgeManifold);
            fprintf(file,'\n');
        end

    end
    fprintf(file,'-1\n');
    fclose(file);
    fclose(file2);

end


%% InitGrid: Initialize base vector for sequential indexing
function [Base] = InitGrid(n, Division) 
  
    Base(1) = 1;
    for i = 2:n+1
        Base(i) = Base(i-1)*Division(i-1);
    end
    return
end

%% GenVert: Calculate perturbed vertices' coordinates and value of F
function [Vert,FVert] = GenVert(n, Grid, Pert, First, Delta, Func) 
  
    Vert  = [];
    FVert = [];
    for i = 0:2^n-1 
        [Coords]    = HyperCubeCoords(n,i);
        [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
        [VHC]       = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
        Vert        = [Vert; VHC];
        [FVHC]      = Func(VHC); 
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

%% HyperCube: Calculates Cartesian coordinates of vertices of hypercube
function [VHC] = HyperCube(n,First,Delta,I,Coords,Pert) 
    for i = 1:n
        VHC(i) = First(i) + (I(i)+Coords(i))*Delta(i) + Pert(i);
    end
    return
end

%% GenFace: Generate labels of k-faces
function [Face,m] = GenFace(n,k)
    if k == 0
        % Generate labels for vertices more efficiently
        Face = zeros(2^n,n);
        f1 = 0:n-1;
        f2 = n:2*n-1;
        for g = 0:2^n-1
            F = flip(dec2bin(g,n) - '0');
            Face(g+1,:) = sort(~F.*f1 + F.*f2);
        end
        m = 2^n;
    else
        Face = []; 
        Division(1:n-k) = 2*n;
        [Base]  = InitGrid(n-k, Division);
        m = 0;
        for g = 0:Base(n-k+1)-1
            [F]      = GridCoords(n-k,g,Base);
            [isface] = IsFace(n,n-k,F);
            if isface == 1
                Face = [Face; F];
                m    = m+1;
            end
        end
    end
    return;
end

%% SplitFaceSimplex: Splits k-face into simplices
function [Face,NumberFace,VertFace,Simp,NumberSimp] = SplitFaceSimplex(n, k)
    
    [Face,NumberFace] = GenFace(n,k);        % Generate face labels
    [VertCube,NumberVert] = GenFace(n,0);    % Generate vertex labels
    
    VertFace = [];
    % Loops through all hypercube k-faces
    for IndFace = 1:NumberFace
        VertFace1 = [];
       
        % Loops through all hypercube vertices
        for IndVert = 1:NumberVert
           % Check whether vertex belongs to face, if so, add to list of vertices in that face
           VertInFace = Include(n-k,Face(IndFace,:),n,VertCube(IndVert,:));
           if VertInFace == n-k
               VertFace1 = [VertFace1 IndVert-1];
           end
        end
        
        % Split k-face into simplices using the vertex list
        % If vertices are in lexicographic order, the k-faces are compatible
        [Simp{IndFace},NumberSimp] = SplitHypercubeSimplex(k,VertFace1);
        VertFace = [VertFace; VertFace1];
    end
end

%% SplitHypercubeSimplex: Splits a hypercube into simplices
function [Simp,nSimp] = SplitHypercubeSimplex(n,VertList,varargin)
  
    % Initialize grid for simplex vertices inside cubes
    for i = 1:n
        DivisionP(i) = i;
    end
    [BaseP] = InitGrid(n,DivisionP);
    nSimp = 0;
    
    % Split hypercube in simplices
    for s = 0:BaseP(n+1)-1

        % Directions traveled in grid from vertex 0 to make simplex 
        % (which is a permutation of coordinate axes) in book notation
        [f] = GridCoords(n,s,BaseP);
        f = f+1;
        % Directions traveled in grid from vertex 0 to make simplex 
        % (which is a permutation of coordinate axes) in usual notation
        [P] = Map_Perm(n,f);

        % Local hypercube indices of vertices in the simplex
        % Every simplex vertex list starts at 0 and ends at 2^n-1 due to how they are generated
        Simp(s+1,:) = GenLabelSimplex(n,P);
        
        % If vertex list is provided (for hypercube face splitting) use the list
        % List must be in lexicographic order
        if nargin > 1
            Simp(s+1,:) = VertList(Simp(s+1,:)+1);
        end
        nSimp = nSimp+1;
    end
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
%     
%     copy2(n) = 0;
%     copy2 = i*ones(size(copy2));
%     aux2 = mod(copy2,Base(2:end));
%     Grid2 = (copy-aux)./Base(2:end);
%     Grid2(1) = i;

    return
end  

%% IsFace: Check if a face is shared by adjacent (n-1)-faces
function [isface] = IsFace(dim,n,F)
  
    for i = 1:n-1
        for j = i+1:n
            % Checks if indices of F are in lexicographic order
            if F(i) >= F(j)
                isface = 0;
                return;        
            end
            % Checks if (n-1)-faces are not opposite to each other
            if mod(F(i),dim) == mod(F(j),dim)
                isface = 0;
                return;        
            end
        end
    end
    isface = 1;
    return;
end

%% Include: Verifies whether (m-1)-face labels y contains the (n-1)-face labels x
function [k] = Include(n,x,m,y)
  
    k = 0;  % Counter for matching labels; if k = n then y contains x
    for i = 1:n
        j = i;
        % Checks whether the face labels of x are in the labels of y (assuming both are in lexicographic order)
        while (j < m) && (y(j) < x(i))           
            j = j+1;
        end
        if y(j) == x(i)
            k = k+1;
        end
    end
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

%% GenLabelSimplex: Find hypercube vertices defining a simplex
function [Simp] = GenLabelSimplex(n,P)
  
    Simp(1) = 0;
    for i = 1:n
        Simp(i+1) = Simp(i)+2^(P(i)-1);
    end
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

%% GetEdgeManifold: Find (k+1)-faces that contain manifold vertices and edges
function [VertexEdgeManifold,FaceEdgeManifold,NumberEdgeManifold,NumberFace] = GetEdgeManifold(n,k,NumberVertex,LabelVertexManifold)
    FaceEdgeManifold   = [];
    NumberEdgeManifold = [];
    VertexEdgeManifold = [];
    NumberFace         = 0;
    
    % Initialize sequential numbering for hypercube (k+1)-faces
    if n-k > 1
        for i = 1:n-k-1
            Division(i) = 2*n;
        end
    else
        Division = 1;
    end
    [Base]  = InitGrid(n-k-1, Division);
    nVertexMax = 2*(k+1)*factorial(k);   % Maximum possible number of vertices in a (k+1)-face
    
    % Loops through all possible (k+1)-face labels 
    for g = 0:Base(n-k)-1
        % Convert sequential numbering to (k+1)-face labels
        [F]      = GridCoords(n-k-1,g,Base);
        % Check if (k+1)-face label corresponds to actual face
        [isface] = IsFace(n,n-k-1,F);
        if isface == 1
            % Checks if face has a manifold vertex
            [LabelVertex,NFace] = InFace(n,k,NumberVertex,LabelVertexManifold,F);
            if NFace > 0
                nVertex = length(LabelVertex);
                LabelVertex = [LabelVertex -ones(1,nVertexMax-nVertex)];
                NumberEdgeManifold = [NumberEdgeManifold; NFace/2];     % Update number of manifold edges in (k+1)-face
                VertexEdgeManifold = [VertexEdgeManifold; LabelVertex]; % Update list of indices (LabelVertexManifold indices)
                FaceEdgeManifold   = [FaceEdgeManifold; F];             % Update face labels that contain edges
                NumberFace         = NumberFace+1;                      % Update counter of k-faces
            end
        end
    end
    return
end

%% InFace: Checks which manifold vertices are in hypercube face 
function [LabelVertex,NumberFace] = InFace(n,k,NumberVertex,LabelVertexManifold,F)
    LabelVertex = [];
    NumberFace  = 0;
    for i = 1:NumberVertex
        % Checks if k-face labels in LabelVertexManifold are in (k+1)-face F
        [j] = Include(n-k-1,F,n-k,LabelVertexManifold(i,:));
        % If egde is in F
        if j == n-k-1
            LabelVertex = [LabelVertex i];  % Include which vertex inside LabelVertexManifold is in 2-face
            NumberFace = NumberFace+1;      % Update number of edges in 2-face F
        end
    end
    return
end

%% SolveAmbiguity: Solves manifold edge ambiguities when there are more than 2 vertices in (k+1)-face
function [VertexEdgeManifold,NumberEdgeManifold] = SolveAmbiguity(n,k,Vert,FVert,FaceVertManifold,NVertManifold,FaceLabel)

    % Find (k+1)-face vertices from label (same as SplitFaceSimplex)
    [FaceVert] = LabelToVert(n,k+1,FaceLabel);
    
    % Split (k+1)-face into simplices
    [Simp,NSimp] = SplitHypercubeSimplex(k+1, FaceVert);

    VertexManifoldInternal = [];
    SimpFacesInternal      = [];
    NVertexTotal = 0;

    % Loop through each simplex
    for indSimp = 1:NSimp
        % Find manifold vertices in each k-face of simplex
        [SimpVertex,SimpFaces,NVertexSimp] = VertexSimplexFaces(n,k,Vert,FVert,Simp(indSimp,:));
        % Add manifold vertices and simplex k-faces to list for sorting later
        if NVertexSimp > 0
            VertexManifoldInternal = [VertexManifoldInternal; SimpVertex];
            SimpFacesInternal = [SimpFacesInternal; SimpFaces];
            NVertexTotal = NVertexTotal + NVertexSimp;
        end
    end

    % Find k-faces in (k+1)-face
    [KFaceLabels,NKFace] = GenFace(n,k);
    KFaceVert = [];
    NKFVert = 0;
    for indKFace = 1:NKFace
       iskface = Include(n-k-1,FaceLabel,n-k,KFaceLabels(indKFace,:));
       if iskface == n-k-1
           KVert = LabelToVert(n,k,KFaceLabels(indKFace,:));
           KFaceVert = [KFaceVert; KVert];
           NKFVert = NKFVert + 1;
       end
    end
    
    % Connect manifold vertices along the inside of the (k+1)-face and find the endpoints to solve ambiguity
    [VertexEdgeManifold,NumberEdgeManifold] = ConnectEndpoints(k,FaceVertManifold,NVertManifold,Simp,SimpFacesInternal,KFaceVert,NSimp,NKFVert,NVertexTotal);    

end

%% LabelToVert: Find hypercube vertices that composes a face label
function [FaceVertex] = LabelToVert(n,k,FaceLabel)
    
    VertCube = GenFace(n,0);
    FaceVertex = [];
    NCVert = 2^n;
    for indVert = 1:NCVert
       VertInFace = Include(n-k,FaceLabel,n,VertCube(indVert,:));
       if VertInFace == n-k
           FaceVertex = [FaceVertex indVert-1];
       end
    end

end

%% VertexSimplexFaces: Find vertices in simplex faces for solving ambiguity 
% This essentially works as a Marching Simplex iteration for a single simplex
function [VertexList,FaceList,NVertex] = VertexSimplexFaces(n,k,Vert,FVert,Simp)
    
    for i = 1:k+1
        DivisionC(i) = i+1;
    end
    [BaseC] = InitGrid(k+1,DivisionC);

    VertexList = [];
    FaceList   = [];
    NVertex    = 0;

    % Loop through all possible simplex faces in the hypercube
    for j = 0:BaseC(k+2)-1
        
        % Make a guess at possible vertex combination that makes a simplex face
        [C] = GridCoords(k+1,j,BaseC);
        C   = C+1;
        % Check if vertex combination is in lexicographic order; if not, discard the combination
        [lex] = Lexico(k+1,C);
        
        % If it is in lexicographic order, check for manifold
        if lex == 1
            % Map hypercube vertices to simplex face
            for i = 1:k+1
                Face(i) = Simp(C(i));
            end
            % Check if (k+1)-face is transversal to manifold
            [Vertex, trans] = GetVertexManifold(n,k,Face,Vert,FVert);
            % If (k+1)-face is transversal, add to list
            if trans == 1
                VertexList = [VertexList; Vertex];
                FaceList   = [FaceList; Face];
                NVertex    = NVertex + 1;
            end
        end

    end

end

%% ConnectEndpoints: Connect manifold edges along a (k+1)-face and output endpoints
function [VertexManifoldEdge,NumberManifoldEdge] = ConnectEndpoints(k,FaceVertManifold,NVertManifold,Simp,SimpFacesVert,KFaceVertex,NCSimp,NCFace,NSFace)
    
% Inputs:
% k: Dimension of face
% FaceVertManifold: 
% NVertManifold
% Simp
% SimpFacesVert
% KFaceVertex
% NCSimp
% NCFace
% NSFace


    VertexManifoldEdge = [];        % List of connected faces
    SimpFacesLeft = SimpFacesVert;  % Remaining simplex faces that contain manifold vertices
    SimpLeft = Simp;                % Remaining simplices that contain manifold vertices
    StartAgain = 1;                 % Flag for initializing from cube face
    NFLeft = NSFace;                % Number of remaining unverified faces 
    NSLeft = NCSimp;                % Number of remaining unverified simplices
    CubeFaceFlag = zeros(NVertManifold,1); % Flags for already verified cube faces
    NumberManifoldEdge = 0;         % Number of manifold edges

    % While list is not exhausted
    while NFLeft > 0
        
        % Get initial vertex that start with face
        if StartAgain == 1
            % Find first simplex face and cube face match
            [~,~,indFace] = LabelMatch(2^k,k+1,KFaceVertex,SimpFacesLeft,NCFace,NFLeft);

            % If there are vertices left not connected to a cube face, throw warning for refinement and return
            if indFace < 0
                warning('Internal loop detected; refinement recommended')
                return
            end

            % Store simplex face that is in a cube face
            ConnSimpFaces = SimpFacesLeft(indFace,:);

            % Store manifold vertex label to connect as edge
            [~,VertexLabel1,~] = LabelMatch(2^k,k+1,FaceVertManifold,SimpFacesLeft(indFace,:),NVertManifold,1,CubeFaceFlag); 

            % Remove simplex face already matched from verification
            SimpFacesLeft(indFace,:) = [];
            StartAgain = 0;
            NFaceConn = 1;
            NFLeft = NFLeft - 1;

            % Flag cube face as already verified
            CubeFaceFlag(VertexLabel1) = 1;
        end

        % Find next simplex that share the face from previous/starting iteration
        [~,indSimp,~] = LabelMatch(k+2,k+1,SimpLeft,ConnSimpFaces(NFaceConn,:),NSLeft,1);

        % Find 2nd face from the simplex that contains a manifold vertex
        [~,~,indFace] = LabelMatch(k+2,k+1,SimpLeft(indSimp,:),SimpFacesLeft,1,NFLeft);
        
        % Connect faces
        ConnSimpFaces = [ConnSimpFaces; SimpFacesLeft(indFace,:)];
        NFaceConn = NFaceConn + 1;

        % Remove new face and verfied simplex from list of unverified faces
        SimpFacesLeft(indFace,:) = [];
        NFLeft = NFLeft - 1;

        % Check if new face is in a cube face
        [FaceMatch,~,~] = LabelMatch(2^k,k+1,KFaceVertex,ConnSimpFaces(NFaceConn,:),NCFace,1);
        % If so, connect manifold vertices as single edge and start again
        if FaceMatch == 1
            [~,VertexLabel2,~] = LabelMatch(2^k,k+1,FaceVertManifold,ConnSimpFaces(NFaceConn,:),NVertManifold,1,CubeFaceFlag); 
            VertexManifoldEdge = [VertexManifoldEdge; VertexLabel1 VertexLabel2];
            NumberManifoldEdge = NumberManifoldEdge + 1;
            StartAgain = 1;
            CubeFaceFlag(VertexLabel2) = 1;
        % If not, remove duplicate of new face
        else
            [~,~,indFace] = LabelMatch(k+2,k+1,SimpLeft(indSimp,:),SimpFacesLeft,1,NFLeft);
            SimpFacesLeft(indFace,:) = [];
            NFLeft = NFLeft - 1;
        end

        % Remove verfied simplex from list of unverified simplices
        SimpLeft(indSimp,:) = [];
        NSLeft = NSLeft - 1;
        
    end

end

%% LabelMatch: Finds first match between lists of labels
function [FaceMatch,indLong,indShort] = LabelMatch(n,k,LongLabelList,ShortLabelList,NLongLabel,NShortLabel,LongLabelFlags,varargin)
    if nargin < 7
        LongLabelFlags = zeros(size(LongLabelList,1));
    end
    FaceMatch = 0;
    for indShort = 1:NShortLabel
        for indLong = 1:NLongLabel
            j = Include(k,ShortLabelList(indShort,:),n,LongLabelList(indLong,:));
            if indLong > length(LongLabelFlags)
                disp('Here be error')
            end
            if j == k && LongLabelFlags(indLong) == 0
                FaceMatch = 1;
                return
            end
        end
    end
    indShort = -1;
    indLong = -1;
    return
end

%% SkeletonHyperCube: Constructs the combinatorial skeleton and print to file
function SkeletonHyperCube(file, n, k, nVertex, nEdge, Vertex, FaceVertex, LabelVertex, Edges)
    % Make connectivity matrix between manifold vertices
    A     = zeros(nVertex);
    for indEdge = 1:nEdge
        A(Edges(indEdge,1),Edges(indEdge,2)) = 1;
        A(Edges(indEdge,2),Edges(indEdge,1)) = 1;
    end
    [nConn,nVertConn,VertConn] = Split_Edges(nVertex,A);
    fprintf(file,'%3d\n',nConn);
    for indConn = 1:nConn
        nVertex = nVertConn(indConn);
        fprintf(file,'%3d\n',nVertex);
        for j = 1:nVertex
            FaceVertexComp(j,:)  = FaceVertex(VertConn(indConn,j),:);
            LabelVertexComp(j,:) = LabelVertex(VertConn(indConn,j),:);
            VertexComp(j,:)      = Vertex(VertConn(indConn,j),:);
            fprintf(file,'%3d ',LabelVertexComp(j,:));
            fprintf(file,'%15.8f ',VertexComp(j,:));
            fprintf(file,'\n');
        end
        [nSkel,Skel,AdjSkel] = SkeletonComp(n, k, nVertex, FaceVertexComp);
        for j = 2:n-k+1
            fprintf(file,'%3d\n',nSkel(j));
            for indEdge = 1:nSkel(j)
                fprintf(file,'%3d ',AdjSkel{j}{indEdge}(:));
                fprintf(file,'\n');
            end
        end
    end
    return
end

%% Split_Edges: Find connected components 
function [nc,nvc,vc] = Split_Edges(nv,A)
    nc   = 0;
    % Find first non-zero row of connectivity matrix A
    [k]  = FindNextComp(nv,A);
    % Loops until all non-zero rows of A are exhausted
    while k ~= 0
        nc             = nc+1;  
        nvc(nc)        = 1;      
        vc(nc,nvc(nc)) = k;
        av             = 1;
        while av <= nvc(nc)
            k = vc(nc,av);
            for i = 1:nv
                if A(k,i) ~= 0
                    [value] = InComp(nvc(nc),i,vc(nc,:));
                    if value == 0
                        nvc(nc)        = nvc(nc)+1;
                        vc(nc,nvc(nc)) = i;
                    end
                    % Set connectivity to zero since we already looped through this one
                    A(k,i) = 0;
                    A(i,k) = 0;
                end
            end
            av = av+1;
        end
        % Find next non-zero row
        [k]  = FindNextComp(nv,A);
    end
    return
end

%% FindNextComp: Find first non-zero row of matrix
function [k] = FindNextComp(n, A)
    for i = 1:n-1
        for j = i+1:n
            if A(i,j) ~= 0
                k = i;
                return
            end
        end
    end
    k = 0;
    return
end

%% InComp: Checks if y is in n-by-1 array x
function [value] = InComp(n, y, x)
    value = 0;
    for i = 1:n
        if x(i) == y
            value = 1;
            return
        end
    end
    return
end

%% SkeletonComp
function [nSkel,Skel,AdjSkel] = SkeletonComp(n, j, nVertex, FaceVertex) 
    nSkel(1) = nVertex;
    Skel{1}  = FaceVertex;
    % Loop through dimensions
    for k = j+1:n-1
        
        % Initialize sequential numbering for k-faces
        for i = 1:n-k
            Division(i) = 2*n;
        end
        [Base]    = InitGrid(n-k, Division);
        FacesSimp = [];
        FacesAdj  = [];
        nSkel(k-j+1)  = 0;
        
        % Loops through all possible k-face labels (given by which (n-1)-faces the k-face belongs to)
        for g = 0:Base(n-k+1)-1
            % Convert sequential numbering to face label
            [Face]    = GridCoords(n-k,g,Base);
            % Check if k-face label corresponds to actual face
            [isface]  = IsFace(n,n-k,Face);
            % If positive, find all (k-1)-faces in skeleton that are part of this k-face
            if isface == 1
                cont = 0;
                Adj  = [];
                % Loops through all (k-1)-faces in skeleton
                for i = 1:nSkel(k-j)
                    % Checks if k-face contains (k-1)-face in the skeleton
                    [np] = Include(n-k,Face,n-k+1,Skel{k-j}(i,:));
                    % If positive, update adjacencies
                    if np == n-k
                        cont      = cont+1;
                        Adj(cont) = i;
                    end
                end
                
                if cont > 0 
                    % Update number of k-faces in the skeleton
                    nSkel(k-j+1)           = nSkel(k-j+1)+1;
                    % Include k-face labels in the face structure
                    FacesSimp          = [FacesSimp; Face];
                    % Include (k-1)-faces adjacent to k-face
                    FacesAdj{nSkel(k-j+1)} = Adj;
                end    
            end
        end
        Skel{k-j+1}    = FacesSimp;
        AdjSkel{k-j+1} = FacesAdj;
    end 
    nSkel(n-j+1)      = 1;
    AdjSkel{n-j+1}{1} = [1:nSkel(n-j)];
    return
end
