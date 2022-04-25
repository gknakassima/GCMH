function ContinuationSimplex(n, k, First, Last, Division, FirstPoint, Func, filename)
   file  = fopen(filename,'w');
   fprintf(file,'%3d %3d\n',n,k);
   fprintf(file,'%3d ',Division(1:n));
   fprintf(file,'\n\n');
   
   Delta = (Last-First)./Division;
   Pert  = random('Normal',0,1,2^n,n)*0.000001;
   Base  = InitGrid(n, Division);
   for i = 1:n
        DivisionP(i) = i;
   end
   [BaseP] = InitGrid(n,DivisionP);
   for i = 1:k+1
        DivisionC(i) = n-k+i;
   end
   [BaseC] = InitGrid(k+1,DivisionC);
    
   NSimpTrans = 0;
   SimpTrans  = [];
   NSimpProc  = 0;
   SimpProc   = [];
   
   [g s value] = GetFirstSimplex(n,k,First,Last,Division,FirstPoint,Func);
   
   if value < 0
       fprintf('Nenhum simplexo transversal\n');
       return
   end
   
   [NSimpTrans SimpTrans] = InsertList(NSimpTrans,SimpTrans,g,s);
   gold                   = -1;
   
   Nsimplex = 0;
   Ncube    = 1;
   while NSimpTrans > 0
        [g s NSimpTrans SimpTrans] = GetAndRemoveList(NSimpTrans,g,SimpTrans);
        [NSimpProc SimpProc]       = InsertList(NSimpProc,SimpProc,g,s);
        Nsimplex = Nsimplex + 1;
        fprintf('Simplex: Ncube = %d  Nsimplex = %d  g = %d  s = %d  N = %d\n',Ncube,Nsimplex,g,s,NSimpTrans);
        if g ~= gold
            Ncube        = Ncube + 1;
            [Grid]       = GridCoords(n,g,Base);
            [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func);
        end
        [P]              = Get_Perm(n,BaseP,s);
        [Simp]           = GenLabelSimplex(n,P); 
        NumberVertex     = 0;
        VertexManifold   = [];
        FaceVertex       = [];
        for j = 0:BaseC(k+2)-1
            [C]   = GridCoords(k+1,j,BaseC);
            C     = C+1;
            [lex] = Lexico(k+1,C);
            if lex == 1
                for i = 1:k+1
                    Face(i) = Simp(C(i));
                end 
                [Vertex trans] = GetVertexManifold(n,k,Face,Vert,FVert);
                if trans == 1
                    VertexManifold = [VertexManifold; Vertex];
                    FaceVertex     = [FaceVertex; Face];
                    NumberVertex   = NumberVertex + 1;
                end
            end 
        end     
        if NumberVertex > 0
           [nSkel Skel AdjSkel] = Skeleton(n,Simp,k,NumberVertex,FaceVertex);  
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
         
           for i = 1:nSkel(n-k)
               [r]           = GetPivotSimplex(n+1,Skel{n-k+1}(1,:),Skel{n-k}(i,:));
               [value g1 s1] = Pivoting(n,r-1,Grid,P,Division,Base,BaseP);
               if value == 1
                   [value] = InList(NSimpProc, SimpProc, g1, s1);
                   if value == 0
                       [NSimpTrans SimpTrans] = InsertList(NSimpTrans,SimpTrans,g1,s1);
                   end
               end
           end
        end
        gold = g;
   end
   fprintf(file,'-1\n');
   
   fclose(file);
   return
end


function [value] = InList(NList, List, g, s)
    value = 0;
    if NList > 0
        for i = 1:NList
            if List(i,1) >= g
                if (List(i,1) == g) && (List(i,2) >= s)
                    if List(i,2) > s
                        return
                    else
                        value = 1;
                        return
                    end
                else
                    if List(i,1) > g
                        return
                    end
                end
            end
        end
    end
    return
end


function [g s NList List] = GetAndRemoveList(NList, g, List)
    if NList == 0
        g  = -1;
        s  = -1;
        return
    end
    for i = 1:NList
        if List(i,1) == g
            s     = List(i,2);
            List  = [List(1:i-1,:); List(i+1:NList,:)];
            NList = NList - 1;
            return
        end
    end
    g     = List(1,1);
    s     = List(1,2);
    List  = List(2:NList,:);
    NList = NList - 1;
    return
end


function [NList List] = InsertList(NList, List, g, s)
    if NList > 0
        for i = 1:NList
            if List(i,1) >= g
                if (List(i,1) == g) && (List(i,2) >= s)
                    if List(i,2) > s
                        List  = [List(1:i-1,:); g s; List(i:NList,:)];
                        NList = NList + 1;
                        return
                    else
                        return
                    end
                else
                    if List(i,1) > g
                        List  = [List(1:i-1,:); g s; List(i:NList,:)];
                        NList = NList + 1;
                        return
                    end
                end
            end
        end
    end
    List  = [List; g s];
    NList = NList + 1;
    return
end

    

function [value g s] = Pivoting(n, i, Grid, P, Division, Base, BaseP)
    if i == 0
        if Grid(P(1)) < Division(P(1))-1
            G       = Grid;
            G(P(1)) = G(P(1)) + 1;
            g       = Get_Label_Prod(n,Base,G);
            Q       = [P(2:n) P(1)];
            s       = Get_Label_Perm(n,BaseP,Q);
            value   = 1;
            return
        else
            g       = 0;
            s       = 0;
            value   = 0;
            return
        end
    end
    if i == n
        if Grid(P(n)) > 0
            G       = Grid;
            G(P(n)) = G(P(n)) - 1;
            g       = Get_Label_Prod(n,Base,G);
            Q       = [P(n) P(1:n-1)];
            s       = Get_Label_Perm(n,BaseP,Q);
            value   = 1;
            return
        else
            g       = 0;
            s       = 0;
            value   = 0;
            return
        end
    end
    g     = Get_Label_Prod(n,Base,Grid);
    Q     = [P(1:i-1) P(i+1) P(i) P(i+2:n)];
    s     = Get_Label_Perm(n,BaseP,Q);
    value = 1;
    return
end
          

function [g s value] = GetFirstSimplex(n, k, First, Last, Division, FirstPoint, Func)
   tol  = 0.000001;
   [F]  = Func(n,FirstPoint);
   if norm(F(1:k),2) > tol
       g     = -1;
       s     = -1;
       value = -1;
       fprintf('Ponto nao pertence a variedade\n');
       return
   end
   for i = 1:n
        if (FirstPoint(i) < First(i)) || (FirstPoint(i) > Last(i))
            g     = -1;
            s     = -1;
            value = -2;
            fprintf('Ponto fora do dominio\n');
            return
        end
   end 
   Delta = (Last-First)./Division;
   Pert  = random('Normal',0,1,2^n,n)*tol;
   Base  = InitGrid(n, Division);
   for i = 1:n
        DivisionP(i) = i;
   end
   [BaseP] = InitGrid(n,DivisionP);    
   Grid    = fix((FirstPoint-First)./Delta);
   g    = Get_Label_Prod(n,Base,Grid);
   s    = 0;
   gold = -1; 
   while g > 0
        fprintf('Simplexo: g = %d  s = %d\n',g,s);
        if g ~= gold
            [Grid]       = GridCoords(n,g,Base);
            [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func);
        end
        [P]     = Get_Perm(n,BaseP,s);
        [Simp]  = GenLabelSimplex(n,P); 
        [trans] = GetBestSimplex(n,Simp,Vert,FirstPoint);
        if trans < 0
            value = 1;
            return
        end 
        gold = g;
        [value g1 s1] = Pivoting(n,trans,Grid,P,Division,Base,BaseP);
        if value == 0
            g  = -1;
        else
            g  = g1;
            s  = s1;
        end 
   end
   g     = -1;
   s     = -1;
   value = -2;
   fprintf('Ponto fora do dominio\n');
   return
end


function [trans] = GetBestSimplex(n,Simp,Vert,FirstPoint)        
   for i = 1:n+1
      A(1,i) = 1;
      A(2:n+1,i) = Vert(Simp(i)+1,:);
   end
   b = [1; FirstPoint(1:n)'];
   lamb = A\b;
   trans = 0;
   if lamb >= 0.0    
       trans = -1;
   else
       min  = lamb(1);
       imin = 1;
       for i = 1:n+1
           if lamb(i) < min
               min  = lamb(i);
               imin = i;
           end
       end
       trans = imin-1;
   end   
   return
end


function [Base] = InitGrid(n, Division) 
   Base(1) = 1;
   for i = 2:n+1
      Base(i) = Base(i-1)*Division(i-1);
   end
   return
end

    
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


function [Coords] = HyperCubeCoords(n,i) 
   copy = i;
   for j = 1:n-1
      Coords(j) = mod(copy,2);
      copy      = (copy-Coords(j))/2;
   end
   Coords(n) = copy;
   return
end


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


function [VHC] = HyperCube(n,First,Delta,I,Coords,Pert) 
   for i = 1:n
      VHC(i) = First(i) + (I(i)+Coords(i))*Delta(i) + Pert(i);
   end
   return
end


function [Vertex trans] = GetVertexManifold(n,k,Face,Vert,FVert)        
   for i = 1:k+1
      A(1,i) = 1;
      A(2:k+1,i) = FVert(Face(i)+1,:);
   end
   b = [1; zeros(k,1)];
   lamb = A\b;
   trans = 0;
   Vertex = zeros(1,n);
   if lamb >= 0      
       for i = 1:k+1
          Vertex = Vertex + lamb(i)*Vert(Face(i)+1,:);
       end
       trans = 1;
   end   
   return
end


function [V FV] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label,Func)     
   [Coords]    = HyperCubeCoords(n,Label);
   [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
   [V]         = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
   [FV]        = Func(n,V); 
   return
end


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

function [Simp] = GenLabelSimplex(n,P)
    Simp(1) = 0;
    for i = 1:n
        Simp(i+1) = Simp(i)+2^(P(i)-1);
    end
end


function [nSkel Skel AdjSkel] = Skeleton(n, Simp, k, nVertex, FaceVertex)   
    nSkel(1) = nVertex;
    Skel{1}  = FaceVertex;
    for s = k+1:n
        for i = 1:s+1
            DivisionC(i) = n-s+i;
        end
        [BaseC]      = InitGrid(s+1,DivisionC);
        FacesSimp    = [];
        FacesAdj     = [];
        nSkel(s-k+1) = 0;
        for j = 0:BaseC(s+2)-1   
          [C]   = GridCoords(s+1,j,BaseC);
          C     = C+1;
          [lex] = Lexico(s+1,C);
          if lex == 1
              for i = 1:s+1
                  Face(i) = Simp(C(i));
              end 
              cont = 0;
              Adj  = [];
              for i = 1:nSkel(s-k)
                  [np] = Include(s,Skel{s-k}(i,:),s+1,Face);
                  if np == s
                      cont      = cont+1;
                      Adj(cont) = i;
                  end
              end
              if cont > 0              
                 nSkel(s-k+1)           = nSkel(s-k+1)+1;
                 FacesSimp              = [FacesSimp; Face];
                 FacesAdj{nSkel(s-k+1)} = Adj;
              end    
           end
        end
        Skel{s-k+1}  = FacesSimp;
        AdjSkel{s-k} = FacesAdj;
    end
    return
end
              

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

function [i] = GetPivotSimplex(n, Simp, Face)
    for i = 1:n
        vertex = Simp(i);
        [k]    = Include(1,vertex,n-1,Face);
        if k == 0
            return
        end
    end
    i = 0;
    return
end

function [f] = Get_Prod(n, Base, labelp)
    [f] = GridCoords(n,labelp,Base);
    return
end

function [p] = Get_Perm(n, Base, labelp)
    labelp = Base(n+1)-1-labelp;
    [f]    = GridCoords(n,labelp,Base);
    f      = f + 1;
    [p]    = Map_Perm(n,f);
    return
end

function [label] = Get_Label_Perm(n, Base, p)
    [f]     = Map_Inv_Perm(n,p);
    [label] = Get_Label_Prod(n,Base,f);
    label   = Base(n+1)-1 - label;
    return
end

function [label] = Get_Label_Prod(n, Base, f)
    label = 0;
    for i = 1:n
        label = label + f(i)*Base(i);
    end
    return
end

function [f] = Map_Inv_Perm(n,p)
    q = p;
    for i = n:-1:2
        [j]  = Find(i,i,q);
        f(i) = j;
        q    = [q(1:j-1) q(j+1:i)];
    end
    f(1) = 1;
    f    = f - 1;
    return
end

function [j] = Find(n,i,x)
    j = 0;
    for k = 1:n
        if x(k) == i
            j = k;
            return
        end
    end
    return
end


