# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:08:25 2023

@author: mhassan22
"""
## HW 2, ME 964: Adv. FEM

#%% Import Libraries
import numpy as np
from scipy.sparse import csr_matrix
from scipy.linalg import solve
import HW2_Functions as HW2 ##User Defined Functions



###### parameters
Dim = 3
Lambda = 6.0e10
mu = 2.0e10           
boundaryX0 = -0.05  
boundaryX1 = 0.05
boundaryY0 = -0.015
boundaryY1 = 0.015   
boundaryZ0 = -0.015
boundaryZ1 = 0.015
numLoadStep =10

W=10 #%cm  % width
H=3 #%cm  % height
L=3 #%cm   % length

E=1000 #% Pa  % youngs modulus
nu=0.3  # % poisson's ratio

rho=1   #% kg/m^-3 % density
delta=np.identity(3)

#%% Create Mesh

# number of divisoin in X Y Z
nElx=20 
nEly=6
nElz=6

nNoEl=8     #%% 8 nodes in linear Hexagoan

nDof=3       #%% Degrees or Freedom


nNodes = (nElx+1)*(nEly+1)*(nElz+1) # total Number of nodes

nEl = nElx*nEly*nElz # Total numner of elements

# total number of elements directionwise
elemL=L/nElz
elemH=H/nEly
elemW=W/nElx

QuadRule = 2 #for Linear Hexahedral elements
#QuadRule = 3 #for Quadratic Hexahedral elements

#% Node coordinates/Matrix
NodesXYZ=np.zeros((nNodes,Dim))


#for i in range(1,nElz+1):
#    for j in range(1,nEly+1):
#        for k in range(1,nElx+1):
#            NodesXYZ[[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1)-1,0]]=(k-1)*elemW  #x-coordinate
#            NodesXYZ[[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1)-1,1]]=(j-1)*elemH  #y-coordinate
#            NodesXYZ[[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1)-1,2]]=(i-1)*elemL  #z-coordinate
            
#%% Global node coordinates
for i in range(nElz):
    for j in range(nEly):
        for k in range(nElx):
            NodesXYZ[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),0]=(k)*elemW #%x-coordinate
            NodesXYZ[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),1]=(j)*elemH # %y-coordinate
            NodesXYZ[k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),2]=(i)*elemL  # %z-coordinate

#%% Connectivity arrays/element matrix
connArr=np.zeros((nEl,nNoEl))

for i in range(1,nElz):
    for j in range(1,nEly):
        for k in range(1,nElx):

            l=k+(j-1)*nElx+(i-1)*nElx*nEly-1

            connArr[l,0]=l+(j-1)+(i-1)*(nElx+nEly+1)
            connArr[l,1]=l+(j)+(i-1)*(nElx+nEly+1)
            connArr[l,2]=l+(j)+nElx+1+(i-1)*(nElx+nEly+1)
            connArr[l,3]=l+(j)+nElx+(i-1)*(nElx+nEly+1)

            connArr[l,4]=l+(j-1)+(i)*(nElx+nEly+1)+nElx*nEly
            connArr[l,5]=l+(j)+(i)*(nElx+nEly+1)+nElx*nEly
            connArr[l,6]=l+(j)+nElx+1+(i)*(nElx+nEly+1)+nElx*nEly
            connArr[l,7]=l+(j)+nElx+(i)*(nElx+nEly+1)+nElx*nEly
            
connArr = np.int_(connArr)
NodeMatrix = np.array( NodesXYZ)
NumberNodes = nNodes
ElementMatrix = connArr
NumberElement  = nEl

#print(NodeMatrix)

#%% initialize Displacement mattrix
disp = np.zeros((NumberNodes*Dim,1))

#%% Dirichlet Boundary 
dirichletX0 = [0]
dirichletX1 = [10]

dirichletY0 = []
dirichletY1 = []
dirichletZ0 =[]
dirichletZ1=[]

dim = Dim

#%% Surface Boundary conditions Update on Nodes

SurfNode_x0 = []
SurfNodeIndex_x0 = []
for index1 in range(1,NumberNodes): 
    if (abs( NodeMatrix[index1-1, 0] - boundaryX0 ) < 1.0e-3):

        #%% For problem 1 Linear
        SurfNode_x0      = [ [SurfNode_x0], [dim * ( index1 - 2 ) + 1]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
        #%Apply Dirichlet boundary condition to respective nodes.
        if (abs( NodeMatrix[index1-1, 1 ] - boundaryY0 ) < 1.0e-3 and abs( NodeMatrix[ index1-1 , 2 ] - boundaryZ0) < 1.0e-3): 
            SurfNode_x0  = [ [SurfNode_x0], [dim * ( index1 - 2 ) + 2],[ dim * ( index1 - 2 ) + 3] ] #%u2=u3=0 at x=0 
        
        if ( abs( NodeMatrix[index1-1, 1 ] - boundaryY0 ) < 1.0e-3 and abs( NodeMatrix[ index1-1, 2 ] - boundaryZ1) < 1.0e-3):
            SurfNode_x0  = [ [SurfNode_x0] , [dim * ( index1 - 2 ) + 2]] #% u2 = 0 at x = 3 e3
                           
                       
        #%% For problem 2 Quadratic
        #SurfNode_x0      = [[SurfNode_x0], [dim * ( index1 - 2 ) + 1], [dim * ( index1 - 2 ) + 2],[dim * ( index1 - 2 ) + 3]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       
        #SurfNodeIndex_x0     = [[SurfNodeIndex_x0],[dim * ( index1 - 2 ) + 1]]
        


SurfNode_x1                  = []
SurfNodeIndex_x1             = []
for index1 in range (NumberNodes):
    if ( abs( NodeMatrix[index1 , 0 ] - boundaryX1 ) < 1.0e-3 ): 
        SurfNode_x1          = [ [SurfNode_x1],[dim * ( index1 - 1 ) + 1 ]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
        SurfNodeIndex_x1     = [ [SurfNodeIndex_x1],[ dim * ( index1 - 1 ) + 1] ] #;%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                  


SurfNode_y0                  = []
SurfNodeIndex_y0             = [] 
for index1 in range (NumberNodes):
    if ( abs( NodeMatrix[index1 , 1 ] - boundaryY0 ) < 1.0e-3): 
        SurfNode_y0          = [[SurfNode_y0],[dim * ( index1 - 1 ) + 2 ]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
        SurfNodeIndex_y0     = [[ SurfNodeIndex_y0],[dim * ( index1 - 1 ) + 2 ]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];


SurfNode_y1                  = []
SurfNodeIndex_y1             = []
for index1 in range(NumberNodes):
    if ( abs( NodeMatrix[index1 , 1 ] - boundaryY1 ) < 1.0e-3 ): 
        SurfNode_y1          = [[ SurfNode_y1],[ dim * ( index1 - 1 ) + 2 ]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
        SurfNodeIndex_y1     = [[SurfNodeIndex_y1],[dim * ( index1 - 1 ) + 2 ]] #%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                   
           
SurfNode_z0 = [] 
SurfNodeIndex_z0 = []
for index1 in range(NumberNodes):
    if ( abs( NodeMatrix[index1 , 2 ] - boundaryZ0 ) < 1.0e-3):
                       SurfNode_z0          = [ [SurfNode_z0], [dim * ( index1 - 1 ) + 3 ]] #;%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                       SurfNodeIndex_z0     = [ [SurfNodeIndex_z0], [dim * ( index1 - 1 ) + 3 ]] #;%DIMS*(index1-1)+2;DIMS*(index1-1)+3];
            
SurfNode_z1                  = []
SurfNodeIndex_z1             = [] 
for index1 in range(NumberNodes):
    if ( abs(NodeMatrix[ index1 , 2 ] - boundaryZ1 ) < 1.0e-3 ): 
        SurfNode_z1          = [[ SurfNode_z1 ],[ dim * ( index1 - 1 ) + 3] ] #;%DIMS*(index1-1)+2;DIMS*(index1-1)+3]; 
        SurfNodeIndex_z1     = [[SurfNodeIndex_z1],[dim * ( index1 - 1 ) + 3 ]] # %DIMS*(index1-1)+2;DIMS*(index1-1)+3];
                                  

#%% Exporting initial data to TXT using PlotVTK function to be exported into MATLAB
currentIncrement = 0

MAXX,MAXY,MAXZ = HW2.plotVTK(NumberNodes,dim,disp,currentIncrement,NodeMatrix,ElementMatrix)

#Disp1 = []
#disp1 = []

#for index1 in range(NumberNodes):
#    disp1 = []
#    for index2 in range(dim):
#        disp1.append(disp[dim*(index1)+index2,0])
#    Disp1.append(disp1)


#OutputVTKFileName = "HW2_vtk_"+str(currentIncrement)

#NodalCoords=NodeMatrix[:,0:2]
#h3=[1,2,3,4,5,6,7,8]
#h3=np.array([i-1 for i in h3])
#ElemConnect=ElementMatrix[:,h3]
#DispArr = np.array(Disp1)
#varargin=DispArr[:,0:dim-1]

## DISCUSS SEARCH LABIB ########################################################################
#####fem_to_vtk_Vector (OutputVTKFileName, NodalCoords, ElemConnect, varargin)

#MAXX=max(DispArr[:,0]); #%To check max
#MAXY=max(DispArr[:,1]); #%To check max
#MAXZ=max(DispArr[:,2]); #%To check max

#%% setting up FE system

## K-global F-global


for currentTime in range(1,numLoadStep):
    currentIncrement += 1
    currentIteration = 0
    res=1
    tol=10^(-13) 
    abs_tol=10^(-15) 
    initialNorm=0.0 
    currentNorm=0.0
    while 1:
        if currentIteration >= 10:
            X=print('No Convergence', currentIncrement, currentIteration)
            break
        if currentNorm > 1/(tol*tol):
            print('No Convergence', currentNorm, currentIncrement, currentIteration )
            break 
        
        #energyTotal                 = 0.0
        FGlobal                  = csr_matrix( ((NumberNodes*dim) , 1),dtype = np.int8).toarray()
        KGlobal                  = csr_matrix( ((NumberNodes*dim) , (NumberNodes*dim)),dtype = np.int8).toarray() #%Zero tHE Global Stiffness Matrix
        
        #%% Setup Elemental Matrices
        for ElementIndex in range (1 , NumberElement):                                                     #%For each element build tHE Klocal matrix,Force/Moment matrix and tHEn assemble tHEm in tHE respective Global matrix.
        #   energy                  = 0.0
            Node                    = ElementMatrix[ElementIndex , :]    #  %Element nodes
            CoordVec                = np.transpose(NodeMatrix[Node.astype(int), : ])                                                #%Element co-ordinates
            nodes                   = len( Node )
            Dispvec1                = np.reshape( disp , (dim,-1), order='F' )
            Dispvec                 = Dispvec1[ : , Node.astype(int) ]
    
            Klocal                  = np.zeros( nodes * dim )                                          # %Zero the local stiffness matrix
            Flocal                  = np.zeros( (nodes * dim , 1) )
           #%% VolumetricIntegral%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         
            if ( QuadRule == 2 ):
                val1 = np.sqrt(1/3) 
                quadpts = np.array([[val1,val1,val1],[-val1,val1,val1], [-val1,-val1,val1], [-val1,val1,-val1], [-val1,-val1,-val1],[val1,-val1,val1],[val1,-val1,-val1],[val1,val1,-val1]])
                quadwts = np.array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125])*8.0
            if (QuadRule ==3): 
                val= np.sqrt(3/5) 
                quadpts =np.array( [ [-val, val, val],[ 0.0, val, val],[ val, val, val],[-val ,val, 0], [0.0, val, 0],[val,val,0],[-val,val,-val],[0.0,val,-val],[val,val,-val],
                          [-val, 0, val],[0.0,0.0,val],[val, 0.0, val],[-val, 0.0, 0],[0.0, 0.0, 0],[val,0.0,0.0],[-val,0.0,-val],[0.0, 0.0, -val],[val,0.0,-val],
                          [-val,-val,val],[0.0,-val,val],[val,-val,val],[-val,-val,0],[0.0,-val,0],[val,-val,0],[-val,-val,-val],[0.0,-val,-val],[val,-val,-val]] )
            
                quadwts =np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])*(8.0/27.0)
            #%% Calculate Jacobian    
            for WeightIndex in range (len( quadwts )):
                N , Nfirst , Jacobian = HW2.GetShapeFunction( quadpts , WeightIndex , CoordVec )
                
                GradU = Dispvec @ np.transpose(Nfirst) 
                Klocal , Flocal = HW2.VolumePhysics(dim,mu,Lambda , nodes, Nfirst, Jacobian, quadwts, WeightIndex, GradU, Klocal, Flocal) # %Physics of Finite Strain implemented

            #%% Assemble Global K and F matrices from Local matrices    
            KGlobal , FGlobal  = HW2.AssembleMatrix( nodes , Klocal , Flocal , Node , dim , KGlobal , FGlobal)
        #%% Modifies the K and F global matrices from the boundary conditions surface and Dirichlet   
            ### CHECK here for errors
        KGlobal, FGlobal = HW2.ApplyBoundaryCondition(KGlobal,FGlobal,SurfNode_x0,SurfNode_x1,SurfNode_y0,SurfNode_y1,SurfNode_z0,SurfNode_z1,dirichletX0, dirichletX1,dirichletY0, dirichletY1,dirichletZ0, dirichletZ1, currentIteration)
        
        #%% Monitoring Resiudual Matrix convergence
        currentNorm  = np.linalg.norm(FGlobal)
        initialNorm  = max(initialNorm,currentNorm)
        res = currentNorm/initialNorm
        print("\nIncrement:",currentIncrement, "\nTime:",currentTime,"\nIteration:" ,currentIteration)
        
        if (currentIteration >1 and ( res<tol or currentNorm<abs_tol)):
            print("Residual converged.")
        
        #print(np.shape(KGlobal))
        #print(np.shape(FGlobal))
        #print(np.shape(FGlobal))
        
        
        #%% Solving for Displacement
        dU =np.linalg.solve(KGlobal,FGlobal)
        disp = disp + dU
        currentIteration += 1

    #%% Exporting Displacement data to TXT using PlotVTK function to be exported into MATLAB        
    MAXX,MAXY,MAXZ = HW2.plotVTK(NumberNodes,dim,disp,currentIncrement,NodeMatrix,ElementMatrix)
          
        