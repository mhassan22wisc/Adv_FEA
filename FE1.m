classdef FE1
   properties (GetAccess=public, SetAccess=public)
       params
       dirichletX0          {mustBeNumeric}
       dirichletX1          {mustBeNumeric}
       dirichletY0          {mustBeNumeric}
       dirichletY1          {mustBeNumeric}
       dirichletZ0          {mustBeNumeric}
       dirichletZ1          {mustBeNumeric}
       NodeMatrix
       NumberNodes          {mustBeNumeric}
       ElementMatrix 
       NumberElement        {mustBeNumeric}
       disp
       currentIncrement     {mustBeNumeric}
       currentIteration     {mustBeNumeric}
       problem              {mustBeNumeric}
       KGlobal
       FGlobal
   end
   
   methods( Static ) 
       %%%Function to read the Mesh
       function [out1, out2]                = readMesh( str1 , str2 )
            opts                            = delimitedTextImportOptions("NumVariables", 4);
            opts.DataLines                  = [1, Inf];
            opts.Delimiter                  = ",";
            opts.VariableNames              = ["VarName1", "e01", "e01_1", "e01_2"];
            opts.VariableTypes              = ["double", "double", "double", "double"];
            opts.ExtraColumnsRule           = "ignore";
            opts.EmptyLineRule              = "read";
            out1                            = table2array(readtable(str1, opts));
            out1(:,1)                       = [];
            clear opts
            opts                            = delimitedTextImportOptions("NumVariables", 9);
            opts.DataLines                  = [1, Inf];
            opts.Delimiter                  = ",";
            opts.VariableNames              = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
            opts.VariableTypes              = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];
            opts.ExtraColumnsRule           = "ignore";
            opts.EmptyLineRule              = "read";
            out2                            = table2array(readtable(str2, opts));
            out2(:,1)                       = [];
       end
       
       %%%%%%%function to store boundary condition for respective DOF's
       function [ SurfNode ]                = boundaryCondition( obj1 )                                         
           if( isempty( obj1.dirichletX0 )  == 0 )
               SurfNode.x0                  = [];
               SurfNodeIndex.x0             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 1 ) - obj1.params.boundaryX0 ) < 1.0e-3 )
                       if (obj1.problem     == 1)
                           SurfNode.x0      = [ SurfNode.x0 ; obj1.params.dim * ( index1 - 1 ) + 1];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                           %Apply Dirichlet boundary condition to respective nodes.
                           if ( abs( obj1.NodeMatrix( index1 , 2 ) - obj1.params.boundaryY0 ) < 1.0e-3 & abs( obj1.NodeMatrix( index1 , 3 ) - obj1.params.boundaryZ0) < 1.0e-3)
                               SurfNode.x0  = [ SurfNode.x0 ; obj1.params.dim * ( index1 - 1 ) + 2; obj1.params.dim * ( index1 - 1 ) + 3]; %u2=u3=0 at x=0
                           end
                           if ( abs( obj1.NodeMatrix( index1 , 2 ) - obj1.params.boundaryY0 ) < 1.0e-3 & abs( obj1.NodeMatrix( index1 , 3 ) - obj1.params.boundaryZ1) < 1.0e-3)
                               SurfNode.x0  = [ SurfNode.x0 ; obj1.params.dim * ( index1 - 1 ) + 2]; % u2 = 0 at x = 3 e3
                           end
                       end
                       if (obj1.problem     == 2)
                           SurfNode.x0      = [ SurfNode.x0 ; obj1.params.dim * ( index1 - 1 ) + 1; obj1.params.dim * ( index1 - 1 ) + 2; obj1.params.dim * ( index1 - 1 ) + 3 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       end
                       SurfNodeIndex.x0     = [ SurfNodeIndex.x0 ; obj1.params.dim * ( index1 - 1 ) + 1 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
           if( isempty( obj1.dirichletX1 )  == 0)
               SurfNode.x1                  = [];
               SurfNodeIndex.x1             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 1 ) - obj1.params.boundaryX1 ) < 1.0e-3 )
                       SurfNode.x1          = [ SurfNode.x1 ; obj1.params.dim * ( index1 - 1 ) + 1 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex.x1     = [ SurfNodeIndex.x1 ; obj1.params.dim * ( index1 - 1 ) + 1 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
           if ( isempty( obj1.dirichletY0 ) == 0 )
               SurfNode.y0                  = [];
               SurfNodeIndex.y0             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 2 ) - obj1.params.boundaryY0 ) < 1.0e-3)
                       SurfNode.y0          = [ SurfNode.y0 ; obj1.params.dim * ( index1 - 1 ) + 2 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex.y0     = [ SurfNodeIndex.y0 ; obj1.params.dim * ( index1 - 1 ) + 2 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
           if( isempty( obj1.dirichletY1 )  == 0)
               SurfNode.y1                  = [];
               SurfNodeIndex.y1             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 2 ) - obj1.params.boundaryY1 ) < 1.0e-3 )
                       SurfNode.y1          = [ SurfNode.y1 ; obj1.params.dim * ( index1 - 1 ) + 2 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex.y1     = [ SurfNodeIndex.y1 ; obj1.params.dim * ( index1 - 1 ) + 2 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
           if( isempty( obj1.dirichletZ0 )  == 0 )
               SurfNode.z0                  = [];
               SurfNodeIndex.z0             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 3 ) - obj1.params.boundaryZ0 ) < 1.0e-3)
                       SurfNode.z0          = [ SurfNode.z0; obj1.params.dim * ( index1 - 1 ) + 3 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex.z0     = [ SurfNodeIndex.z0 ; obj1.params.dim * ( index1 - 1 ) + 3 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
           if( isempty( obj1.dirichletZ1 )  == 0 )
               SurfNode.z1                  = [];
               SurfNodeIndex.z1             = [];
               for index1                   = 1 : obj1.NumberNodes
                   if ( abs( obj1.NodeMatrix( index1 , 3 ) - obj1.params.boundaryZ1 ) < 1.0e-3 )
                       SurfNode.z1          = [ SurfNode.z1 ; obj1.params.dim * ( index1 - 1 ) + 3 ];%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex.z1     = [ SurfNodeIndex.z1 ; obj1.params.dim * ( index1 - 1 ) + 3 ]; %DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   end
               end 
           end
           
       end
       
       %%Function for quadrature rule
       function [ pts, wts ]                = obtainQuadRule( QuadRule )
           if ( QuadRule                    == 2 )
               val1                         = sqrt(1/3);
               pts                          = [val1 val1 val1; -val1 val1 val1; -val1 -val1 val1; -val1 val1 -val1; -val1 -val1 -val1; 
                                               val1 -val1 val1; val1 -val1 -val1; val1 val1 -val1];
               wts                          = [0.125; 0.125; 0.125; 0.125; 0.125; 0.125; 0.125; 0.125]*8.0;
           end
           if (QuadRule ==3)
               val= sqrt(3/5)
               pts =[ -val val val; 0.0 val val; val val val;   -val val 0; 0.0 val 0; val val 0;   -val val -val; 0.0 val -val; val val -val;
                      -val 0 val; 0.0 0.0 val; val 0.0 val;   -val 0.0 0; 0.0 0.0 0; val 0.0 0.0;   -val 0.0 -val; 0.0 0.0 -val; val 0.0 -val;
                      -val -val val; 0.0 -val val; val -val val;   -val -val 0; 0.0 -val 0; val -val 0;   -val -val -val; 0.0 -val -val; val -val -val];
                  wts =[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]*8.0/27.0;
           end
       end
       
       %Function for matrix assembly
       function [ KGLOBAL , FGLOBAL ]       = AssembleMatrix( Length , Klocal, Flocal, Node, DIMS, KGLOBAL, FGLOBAL)
           for index1 = 1 : Length
               for index2 = 1 : Length
                   for DIMS1 = 1 : DIMS
                       for DIMS2 = 1 : DIMS
                           KGLOBAL(DIMS*(Node(index1)-1)+DIMS1,DIMS*(Node(index2)-1)+DIMS2)=KGLOBAL(DIMS*(Node(index1)-1)+DIMS1,DIMS*(Node(index2)-1)+DIMS2)+Klocal(DIMS*(index1-1)+DIMS1,DIMS*(index2-1)+DIMS2);
                       end
                   end
               end
               for DIMS1 = 1 : DIMS
                   FGLOBAL( DIMS *( Node( index1 ) - 1 ) + DIMS1 ) = FGLOBAL( DIMS * ( Node( index1 ) -1 ) + DIMS1 ) + Flocal( DIMS * (index1-1) + DIMS1);
               end
           end
       end
       
       %Function to obtain shape function at every quadrature point
       function [N, Nfirst, Jacobian] = GetShapeFunction(pts, WeightIndex, CoordVec)

            x           = pts( WeightIndex , 1 ); y = pts( WeightIndex , 2 ); z = pts( WeightIndex , 3 ); 
            N           = [0.125*(1-x)*(1-y)*(1-z), 0.125*(1+x)*(1-y)*(1-z), 0.125*(1+x)*(1+y)*(1-z), 0.125*(1-x)*(1+y)*(1-z), 
                          0.125*(1-x)*(1-y)*(1+z), 0.125*(1+x)*(1-y)*(1+z), 0.125*(1+x)*(1+y)*(1+z), 0.125*(1-x)*(1+y)*(1+z) ]; % 8-noded linear Hex element.
            dN          = [-0.125*(1-y)*(1-z), -0.125*(1-x)*(1-z), -0.125*(1-y)*(1-x); 
                            0.125*(1-y)*(1-z), -0.125*(1+x)*(1-z), -0.125*(1-y)*(1+x);
                            0.125*(1+y)*(1-z), 0.125*(1+x)*(1-z), -0.125*(1+x)*(1+y);
                            -0.125*(1+y)*(1-z), 0.125*(1-x)*(1-z), -0.125*(1+y)*(1-x);
                            -0.125*(1-y)*(1+z), -0.125*(1-y)*(1+z), 0.125*(1-x)*(1-y);
                            0.125*(1-y)*(1+z), -0.125*(1+x)*(1+z), 0.125*(1+x)*(1-y);
                            0.125*(1+y)*(1+z), 0.125*(1+x)*(1+z), 0.125*(1+x)*(1+y);
                            -0.125*(1+y)*(1+z), 0.125*(1-x)*(1+z), 0.125*(1-x)*(1+y)]';
            Jacobian    = CoordVec * dN';
            Nfirst      = [ inv( Jacobian' )] * dN;
       end
       
       %Function to apply boundary conditions on the assembled matrix
       function[ KGlobal1 , FGlobal1 ] = ApplyBoundaryCondition( obj1 , SurfNode)
            if( isempty( obj1.dirichletX0 ) == 0 )
            %%%%%%%%%FIX B.C. AT X==0
                    for index1 = 1 : length( SurfNode.x0 )
                        if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal-obj1.dirichletX0*obj1.KGlobal(:,SurfNode.x0(index1));
                        end
                        obj1.KGlobal( SurfNode.x0( index1 ) , : )                   = 0.0;   
                        obj1.KGlobal( : , SurfNode.x0( index1 ) )                   = 0.0;
                        obj1.KGlobal( SurfNode.x0( index1 ) , SurfNode.x0( index1 ))= 1.0;
                        if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.x0( index1 ) , 1 )               = obj1.dirichletX0;
                        else 
                            obj1.FGlobal( SurfNode.x0( index1 ) , 1 )               = 0.0;
                        end
                    end        
            end
            if( isempty( obj1.dirichletX1 ) == 0 )
            %%%%%%%%%FIX B.C. AT X==1
                    for index1 = 1 : length( SurfNode.x1 )
                         if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal - obj1.dirichletX1 * obj1.KGlobal( : , SurfNode.x1( index1 ));
                         end 
                         obj1.KGlobal( SurfNode.x1( index1 ) , : )                  = 0.0;   
                         obj1.KGlobal( : , SurfNode.x1( index1 ))                   = 0.0;
                         obj1.KGlobal( SurfNode.x1( index1 ), SurfNode.x1( index1 ))= 1.0;
                         if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.x1( index1 ) , 1 )               = obj1.dirichletX1;
                         else 
                            obj1.FGlobal( SurfNode.x1 ( index1 ) , 1 )              = 0.0;
                         end
                    end        
            end

            if( isempty( obj1.dirichletY0 ) == 0 )
            %%%%%%%%%FIX B.C. AT y==0
                    for index1 = 1 : length( SurfNode.y0 )
                        if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal - obj1.dirichletY0 * obj1.KGlobal( : , SurfNode.y0( index1 ));
                         end
                        obj1.KGlobal( SurfNode.y0( index1 ) , : )                   = 0.0;   
                        obj1.KGlobal( : , SurfNode.y0( index1 ))                    = 0.0;
                        obj1.KGlobal( SurfNode.y0( index1 ) , SurfNode.y0( index1 ))= 1.0;
                        if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.y0( index1 ) , 1 )               = obj1.dirichletY0;
                        else 
                            obj1.FGlobal( SurfNode.y0 ( index1 ) , 1 ) = 0.0;
                        end
                    end        
            end
            % % 

            if( isempty( obj1.dirichletY1 ) == 0 )
            %%%%%%%%%FIX B.C. AT y==1
                    for index1 = 1 : length( SurfNode.y1 )
                        if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal - obj1.dirichletY1 * obj1.KGlobal( : , SurfNode.y1( index1 ));
                         end
                        obj1.KGlobal( SurfNode.y1( index1 ) , : )                   = 0.0;   
                        obj1.KGlobal( : , SurfNode.y1( index1 ))                    = 0.0;
                        obj1.KGlobal( SurfNode.y1( index1 ) , SurfNode.y1( index1 ))= 1.0;
                        if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.y1( index1 ) , 1 )               = obj1.dirichletY1;
                        else 
                            obj1.FGlobal( SurfNode.y1( index1 ) , 1 )               = 0.0;
                        end
                    end        
            end


            if( isempty( obj1.dirichletZ0 ) == 0 )
                %%%%%%%%%FIX B.C. AT z==0
                    for index1 = 1 : length( SurfNode.z0 )
                        if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal-obj1.dirichletZ0 * obj1.KGlobal( : , SurfNode.z0( index1 ));
                         end
                        obj1.KGlobal( SurfNode.z0( index1 ) , : )                   = 0.0;   
                        obj1.KGlobal( : , SurfNode.z0( index1 ))                    = 0.0;
                        obj1.KGlobal( SurfNode.z0( index1 ) , SurfNode.z0( index1 ))= 1.0;
                        if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.z0( index1 ) , 1 )               = obj1.dirichletZ0;
                        else 
                            obj1.FGlobal( SurfNode.z0( index1 ) , 1 )               = 0.0;
                        end
                    end        
            end

            if( isempty( obj1.dirichletZ1 ) == 0 )
            %%%%%%%%%FIX B.C. AT z==1
                    for index1 = 1 : length( SurfNode.z1 )
                        if obj1.currentIteration == 0
                            obj1.FGlobal                                            = obj1.FGlobal-obj1.dirichletZ1 * obj1.KGlobal( : , SurfNode.z1( index1 ));
                         end
                        obj1.KGlobal( SurfNode.z1( index1 ) , : )                   = 0.0;   
                        obj1.KGlobal( : , SurfNode.z1( index1 ))                    = 0.0;
                        obj1.KGlobal( SurfNode.z1( index1 ) , SurfNode.z1( index1 ))= 1.0;
                        if obj1.currentIteration == 0
                            obj1.FGlobal( SurfNode.z1( index1 ) , 1)                = obj1.dirichletZ1;
                        else 
                            obj1.FGlobal( SurfNode.z1( index1 ) , 1 )               = 0.0;
                        end
                    end        
            end
            KGlobal1 = obj1.KGlobal;
            FGlobal1 = obj1.FGlobal;

            end

       
   end
end