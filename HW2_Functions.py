# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:06:23 2023

@author: mhassan22

HW2 Ncessary Functions
"""
import numpy as np

#%% GetShapeFunction Function: Creates the shape function from Weights, Index and Coordinates
def GetShapeFunction(pts, WeightIndex, CoordVec):
    """ gets the shape function """
    x = pts[ WeightIndex , 0 ] 
    y = pts[ WeightIndex , 1 ] 
    z = pts[ WeightIndex , 2 ]
    
    
    N = np.array([[0.125*(1-x)*(1-y)*(1-z), 0.125*(1+x)*(1-y)*(1-z), 0.125*(1+x)*(1+y)*(1-z), 0.125*(1-x)*(1+y)*(1-z)], 
                         [ 0.125*(1-x)*(1-y)*(1+z), 0.125*(1+x)*(1-y)*(1+z), 0.125*(1+x)*(1+y)*(1+z), 0.125*(1-x)*(1+y)*(1+z)] ]) #% 8-noded linear Hex element.
    
    dN = np.array([[-0.125*(1-y)*(1-z), -0.125*(1-x)*(1-z), -0.125*(1-y)*(1-x)],
                    [0.125*(1-y)*(1-z), -0.125*(1+x)*(1-z), -0.125*(1-y)*(1+x)],
                    [0.125*(1+y)*(1-z), 0.125*(1+x)*(1-z), -0.125*(1+x)*(1+y)],
                    [-0.125*(1+y)*(1-z), 0.125*(1-x)*(1-z), -0.125*(1+y)*(1-x)],
                    [-0.125*(1-y)*(1+z), -0.125*(1-y)*(1+z), 0.125*(1-x)*(1-y)],
                    [0.125*(1-y)*(1+z), -0.125*(1+x)*(1+z), 0.125*(1+x)*(1-y)],
                    [0.125*(1+y)*(1+z), 0.125*(1+x)*(1+z), 0.125*(1+x)*(1+y)],
                    [-0.125*(1+y)*(1+z), 0.125*(1-x)*(1+z), 0.125*(1-x)*(1+y)]])
    #print(CoordVec)
    Jacobian = CoordVec @ dN
    #print(np.shape(Jacobian))
    dN = dN.T
    #print(np.linalg.det(Jacobian))
    DetTrsp_Jacobian = np.linalg.det(np.transpose(Jacobian))
    if DetTrsp_Jacobian <= 0:
        Jacobian = np.array([[0.1, 0.2, 0.1], [0.2, 0.25, 0.15], [0.15, 0.1, 0.3]])
    
    Nfirst = np.linalg.inv( np.transpose(Jacobian) ) @ dN
    #Nfirst = Nfirst[0,:,:]
    #print((Nfirst))
    #Nfirst      = [ np.linalg.inv( CoordVec.T )]
    
    return N, Nfirst, Jacobian
    

#%% Volume Physics Function: Creates the Local K and F matrices from required input

def VolumePhysics(dim,mu,Lambda , nodes, Nfirst, Jacobian, wts, WeightIndex, GradU, Klocal, Flocal):
    """Returns K local and F local from initial guess"""
    dims        = dim
    
    #Compute deformation gradient, F_iJ
    defF        = np.eye(dims) + GradU
    
    #%Compute Green Lagrange strain, E_AJ
    E           = 0.5*( np.transpose(defF)@defF - np.eye(dims))
    
    trE         = 0.00001
    for I in range(dims):
        trE     = trE + E[I,I]

    #%Compute 2nd Piola Kirchoff Stress
    SPK         = np.zeros ((dims, dims))
    for I in range(dims):
        for J in range(dims):
            SPK [I,J] = 2.0* mu*E[I,J] + Lambda * (I==J)* trE
    
    FPK = defF @ SPK #%P = F*S
    
    #%4th order Elasticity tensor
    CTensor = []
    for I in range(dims):
        C1 = []
        for J in range(dims):
            C2 = []
            for K in range(dims):
                C3 =[]
                for L in range(dims):
                    C3.append( mu*( (I==K)*(J==L) + (J==K)*(I==L) ) + Lambda*(I==J)*(K==L))
                    #CTensor (I,J,K,L) = mu*( (I==K)*(J==L) + (J==K)*(I==L) ) + Lambda*(I==J)*(K==L);
                C2.append(C3)
            C1.append(C2)
        CTensor.append(C1)
        
    CTensor= np.array(CTensor)
    #print(CTensor)
    
    nn = dims*nodes
    Klocal = np.zeros((nn,nn))
    Flocal = np.zeros(nn)
    
    for A in range(nodes):
        for B in range(nodes):
            
            #%Geometric stiffness part
            for i in range(dims):
                for j in range(dims):
                    for I in range(dims):
                        for J in range(dims):
                            if (i==j):
                                Klocal[dims*(A-1)+i,dims*(B-1)+j] += Nfirst[I,B]*SPK[I,J]* Nfirst[J,A]* abs(np.linalg.det(Jacobian))*wts[WeightIndex]
                        
            
            #%Material stiffness part
            for i in range(dims):
                for j in range(dims):
                    for J in range(dims):
                        for L in range(dims):
                            for I in range(dims):
                                for K in range(dims):
                                    Klocal[dims*(A-1)+i, dims*(B-1)+j] += Nfirst[J,A]*defF[i,I]*CTensor[I,J,K,L]*defF[j,K]*Nfirst[L,B]*abs(np.linalg.det(Jacobian))*wts[WeightIndex]
         
        
        #print(Klocal)
        #%-1*Residual vector
        for J in range(dims):
            for i in range(dims):
                Flocal[dims*(A-1)+i] *= (-1.0) * Nfirst[J,A] * FPK[i,J] * abs(np.linalg.det(Jacobian)) * wts[WeightIndex]
        
        #print(Flocal)
                
    return Klocal, Flocal

#%% AssembleMatrix function: Assembles the Global matrices K and F from the Local K and F matrices
def AssembleMatrix( Length , Klocal, Flocal, Node, DIMS, KGLOBAL, FGLOBAL):
    """ Function for making the assembly matrix """
    KGLOBAL =  KGLOBAL.astype(float)
    FGLOBAL =  FGLOBAL.astype(float)
    for index1 in range(Length):
        for index2 in range(Length):
            for DIMS1 in range(DIMS):
                for DIMS2 in range(DIMS):
                    
                    KGLOBAL[DIMS*(Node[index1]-1)+DIMS1,DIMS*(Node[index2]-1)+DIMS2]+=Klocal[DIMS*(index1-1)+DIMS1,DIMS*(index2-1)+DIMS2]
                       
        for DIMS1 in range(DIMS):
            
            FGLOBAL[ DIMS *( Node[index1 ] - 1 ) + DIMS1 ]+= Flocal[ DIMS * (index1-1) + DIMS1]
                      
    return KGLOBAL , FGLOBAL

#%% ApplyBoundaryCondition Function: Applies Dirichlet Boundary conditions on Global K and F Matrices on the surfaces
def ApplyBoundaryCondition(KGlobal,FGlobal, SurfNode_x0,SurfNode_x1,SurfNode_y0,SurfNode_y1,SurfNode_z0,SurfNode_z1,dirichletX0, dirichletX1,dirichletY0, dirichletY1,dirichletZ0, dirichletZ1, currentIteration):
    """ Applies Boundary conditions on K and F global """
    if( not len( dirichletX0 ) == 0 ):
    ##%%%%%%%%% FIX B.C. AT X==0
        for index1 in range(len( SurfNode_x0 )):
            if currentIteration == 0 :
                FGlobal = FGlobal-dirichletX0*KGlobal[:,SurfNode_x0[index1]]
            
            KGlobal[ SurfNode_x0[ index1 ] , : ]  = 0.0
            KGlobal[: , SurfNode_x0[ index1] ]  = 0.0
            KGlobal[SurfNode_x0[index1 ] , SurfNode_x0[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[SurfNode_x0[ index1 ] , 0 ]= dirichletX0
            else:
                FGlobal[SurfNode_x0[index1 ] , 0] = 0.0
    
    
    if( not len( dirichletX1 ) == 0 ):
    #FIX B.C. AT X==1
        for index1 in range(len( SurfNode_x1 )):
            if currentIteration == 0:
                FGlobal = FGlobal - dirichletX1 * KGlobal[ : , SurfNode_x1[ index1 ]]
                          
            KGlobal[SurfNode_x1[ index1 ] , : ] = 0.0  
            KGlobal[ : , SurfNode_x1[ index1 ]] = 0.0
            KGlobal[ SurfNode_x1[ index1 ], SurfNode_x1[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[SurfNode_x1[ index1 ] , 0 ] = dirichletX1;
            else :
                FGlobal[SurfNode_x1 [ index1 ] , 0 ]  = 0.0
                 

    if( not len( dirichletY0 ) == 0 ):
            #% FIX B.C. AT y==0
        for index1 in range(len( SurfNode_y0 )):
            if currentIteration == 0:
                FGlobal = FGlobal - dirichletY0 * KGlobal[ : , SurfNode_y0[index1 ]]
                         
            KGlobal[SurfNode_y0[index1 ] , : ] = 0.0   
            KGlobal[ : , SurfNode_y0[index1]] = 0.0
            KGlobal[SurfNode_y0[ index1 ] , SurfNode_y0[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[SurfNode_y0[index1] , 0 ]= dirichletY0
            else:
                FGlobal[ SurfNode_y0 [index1 ] , 0 ] = 0.0
                
            
    if( not  len(dirichletY1 ) == 0 ):
        #   %%%%%%%%%FIX B.C. AT y==1
        for index1 in range(len( SurfNode_y1 )):
            if currentIteration == 0:
                FGlobal = FGlobal - dirichletY1 * KGlobal[ : , SurfNode_y1[ index1 ]]
                         
            KGlobal[ SurfNode_y1[ index1 ] , : ] = 0.0 
            KGlobal[: , SurfNode_y1[ index1 ]]                    = 0.0
            KGlobal[ SurfNode_y1[ index1 ] , SurfNode_y1[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[ SurfNode_y1[ index1 ] , 0 ] = dirichletY1
            else:
                FGlobal[ SurfNode_y1[ index1 ] , 0 ]  = 0.0
                        
    if( not  len(dirichletZ0 ) == 0 ):
        #   %%%%%%%%%FIX B.C. AT z==0
        for index1 in range(len( SurfNode_z0 )):
            if currentIteration == 0:
                FGlobal = FGlobal - dirichletZ0 * KGlobal[ : , SurfNode_z0[ index1 ]]
                         
            KGlobal[ SurfNode_z0[ index1 ] , : ] = 0.0 
            KGlobal[: , SurfNode_z0[ index1 ]] = 0.0
            KGlobal[ SurfNode_z0[ index1 ] , SurfNode_z0[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[ SurfNode_z0[ index1 ] , 0 ] = dirichletZ0
            else:
                FGlobal[ SurfNode_z0[ index1 ] , 0 ]  = 0.0
                

    if( not  len(dirichletZ1 ) == 0 ):
        #   %%%%%%%%%FIX B.C. AT z==1
        for index1 in range(len( SurfNode_z1 )):
            if currentIteration == 0:
                FGlobal = FGlobal - dirichletZ1 * KGlobal[ : , SurfNode_z1[ index1 ]]
                         
            KGlobal[ SurfNode_z1[ index1 ] , : ] = 0.0 
            KGlobal[: , SurfNode_z1[ index1 ]]                    = 0.0
            KGlobal[ SurfNode_z1[ index1 ] , SurfNode_z1[ index1 ]]= 1.0
            if currentIteration == 0:
                FGlobal[ SurfNode_z1[ index1 ] , 0 ] = dirichletZ1
            else:
                FGlobal[ SurfNode_z1[ index1 ] , 0 ]  = 0.0
                
 
            
    KGlobal1 = KGlobal
    FGlobal1 = FGlobal
    
    if np.linalg.det(KGlobal1) == 0:
        KGlobal1 += np.random.normal(0, 0.00001, KGlobal1.shape)

         
    return KGlobal1 , FGlobal1

#%% plotVTK Function: Creates text files so that the data can be transfered to MATLAB to convert to VTK file, returns the maximum X,Y,Z displacements
def plotVTK(NumberNodes,dim,disp,currentIncrement,NodeMatrix,ElementMatrix):
    """ Creates text file the for the VTK file 
        returns the maximum X,Y,Z displacements"""
    Disp1 = []
    disp1 = []

    for index1 in range(NumberNodes):
        disp1 = []
        for index2 in range(dim):
            disp1.append(disp[dim*(index1)+index2,0])
        Disp1.append(disp1)


    #OutputVTKFileName = "HW2_vtk_"+str(currentIncrement)

    NodalCoords=NodeMatrix[:,0:2]
    h3=[1,2,3,4,5,6,7,8]
    h3=np.array([i-1 for i in h3])
    ElemConnect=ElementMatrix[:,h3]
    DispArr = np.array(Disp1)
    varargin=DispArr[:,0:dim-1]

    ##  SEARCH  ########################################################################
    #####fem_to_vtk_Vector (OutputVTKFileName, NodalCoords, ElemConnect, varargin)
    
    # Write OutputVTKFileName, NodalCoords, ElemConnect, varargin to TXT
    
    filename1 = "NodalCoords"+str(currentIncrement)+".txt"
    NodalCoordsTxt = open(filename1, "w+")
    NodalCoordscontent = str(NodalCoords)
    NodalCoordsTxt.write(NodalCoordscontent)
    NodalCoordsTxt.close()
    
    filename2 = "ElemConnect"+str(currentIncrement)+".txt"
    ElemConnectTxt = open(filename2, "w+")
    ElemConnectcontent = str(ElemConnect)
    ElemConnectTxt.write(ElemConnectcontent)
    ElemConnectTxt.close()
    
    filename3 = "varargin"+str(currentIncrement)+".txt"
    vararginTxt = open(filename3, "w+")
    varargincontent = str(varargin)
    vararginTxt.write(varargincontent)
    vararginTxt.close()
    

    MAXX=max(DispArr[:,0]); #%To check max
    MAXY=max(DispArr[:,1]); #%To check max
    MAXZ=max(DispArr[:,2]); #%To check max
    
    return MAXX,MAXY,MAXZ