%% ME964: HW 1
% Mahmudul Hassan
% Spring 2023
%% Problem 1: Part 1, 2

% clear everything
clear all;
clc;
close all;

%% Variables

W=0.1; %m  % width
H=0.1; %m  % height
L=1; %m   % length

E=1000; % Pa  % youngs modulus
nu=0.3;   % poisson's ratio

rho=1;   % kg/m^-3 % density
delta=eye(3);



%% Create Mesh

% number of divisoin in X Y Z
nElx=4; 
nEly=4; 
nElz=40;


nNoEl=8;     %% 8 nodes in linear Hexagoan

nDof=3;        %% Degrees or Freedom 


nNodes = (nElx+1)*(nEly+1)*(nElz+1);
    
nEl = nElx*nEly*nElz;

% Node coordinates
NodesXYZ=zeros(nNodes,3);

elemL=L/nElz;
elemH=H/nEly;
elemW=W/nElx;

    for i=1:nElz+1
        for j=1:nEly+1
            for k=1:nElx+1
                NodesXYZ(k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),1)=(k-1)*elemW;  %x-coordinate
                NodesXYZ(k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),2)=(j-1)*elemH;  %y-coordinate
                NodesXYZ(k+(j-1)*(nElx+1)+(i-1)*(nElx+1)*(nEly+1),3)=(i-1)*elemL;  %z-coordinate
            end
        end
    end

%  x_co = linspace(0,0.1,nElx);
%  y_co = linspace(0,0.1,nEly);
%  z_co = linspace(0,1,nElz);
% 
%  [X_co,Y_co,Z_co] = meshgrid(x_co,y_co,z_co);
% 
%  NodesXYZ = [X_co',Y_co',Z_co'];

    
% Connectivity arrays
connArr=zeros(nEl,8);
    
    for i=1:nElz
        for j=1:nEly
            for k=1:nElx

                l=k+(j-1)*nElx+(i-1)*nElx*nEly;
                
                connArr(l,1)=l+(j-1)+(i-1)*(nElx+nEly+1);
                connArr(l,2)=l+(j)+(i-1)*(nElx+nEly+1);
                connArr(l,3)=l+(j)+nElx+1+(i-1)*(nElx+nEly+1);
                connArr(l,4)=l+(j)+nElx+(i-1)*(nElx+nEly+1);
                
                connArr(l,5)=l+(j-1)+(i)*(nElx+nEly+1)+nElx*nEly;
                connArr(l,6)=l+(j)+(i)*(nElx+nEly+1)+nElx*nEly;
                connArr(l,7)=l+(j)+nElx+1+(i)*(nElx+nEly+1)+nElx*nEly;
                connArr(l,8)=l+(j)+nElx+(i)*(nElx+nEly+1)+nElx*nEly;
                
            end
        end
    end

%% Gauss Quadrature 
% % for n=2
% 
% gmat = 2;
% gaussMatrix(1,1:2)=[-1/sqrt(3) 1];
% gaussMatrix(2,1:2)=[1/sqrt(3) 1];

% for n=3

gmatn = 3;
GaussQuadMat(1,1:2)=[-sqrt(3/5) 5/9];
GaussQuadMat(2,1:2)=[0 8/9];
GaussQuadMat(3,1:2)=[sqrt(3/5) 5/9];

%% C_ijkl 

%%% From the given parameters
lambda=E/2/(1+nu);
mu=E*nu/(1+nu)/(1-2*nu);


C=zeros(3,3,3,3);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l)=lambda*delta(i,j)*delta(k,l)  + 2*mu*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
            end
        end
    end
end


 
%%  Local and Global matrix formation

kGlobal=zeros(nDof*nNodes,nDof*nNodes);

mGlobal=zeros(nDof*nNodes,nDof*nNodes);

fGlobal=zeros(nDof*nNodes,1);

%%% Local Matrices

for e=1:nEl
    kLocal=zeros(nDof*nNoEl,nDof*nNoEl);
    mLocal=zeros(nDof*nNoEl,nDof*nNoEl);
    fLocal=zeros(nDof*nNoEl,1);
    
    for wqX=1:gmatn
        for wqY=1:gmatn
            for wqZ=1:gmatn
                
                %[N,dNdz,dNde,dNdg] = shapeFunc(gaussMatrix(wqX,1),gaussMatrix(wqY,1),gaussMatrix(wqZ,1));
                
                % shape functions
                N(1)=(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                N(2)=(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                N(3)=(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                N(4)=(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                N(5)=(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                N(6)=(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                N(7)=(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                N(8)=(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                
                % shape function derivatives
                dNdz(1)=-(1-GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNdz(2)=(1-GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNdz(3)=(1+GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNdz(4)=-(1+GaussQuadMat(wqY,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNdz(5)=-(1-GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNdz(6)=(1-GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNdz(7)=(1+GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNdz(8)=-(1+GaussQuadMat(wqY,1))*(1+GaussQuadMat(wqZ,1))/8;


                dNde(1)=-(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNde(2)=-(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNde(3)=(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNde(4)=(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqZ,1))/8;
                dNde(5)=-(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNde(6)=-(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNde(7)=(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqZ,1))/8;
                dNde(8)=(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqZ,1))/8;

                dNdg(1)=-(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))/8;
                dNdg(2)=-(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))/8;
                dNdg(3)=-(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))/8;
                dNdg(4)=-(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))/8;
                dNdg(5)=(1-GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))/8;
                dNdg(6)=(1+GaussQuadMat(wqX,1))*(1-GaussQuadMat(wqY,1))/8;
                dNdg(7)=(1+GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))/8;
                dNdg(8)=(1-GaussQuadMat(wqX,1))*(1+GaussQuadMat(wqY,1))/8;



                Jac=[dNdz*NodesXYZ(connArr(e,:),1)   dNde*NodesXYZ(connArr(e,:),1)  dNdg*NodesXYZ(connArr(e,:),1)  ;  
                    dNdz*NodesXYZ(connArr(e,:),2)   dNde*NodesXYZ(connArr(e,:),2)  dNdg*NodesXYZ(connArr(e,:),2)  ;
                    dNdz*NodesXYZ(connArr(e,:),3)   dNde*NodesXYZ(connArr(e,:),3)  dNdg*NodesXYZ(connArr(e,:),3) ];
                
                dNdX=[dNdz' dNde' dNdg']/Jac;
                
                
                %%% Stiffness matrix
                for A=1:nNoEl
                    for B=1:nNoEl
                        
                        %K=nodalStiffness(dNdX,C,A,B);

                        K=zeros(3,3);
                        dNdX_t=dNdX';

                        for i=1:3
                            for j=1:3
                                for k=1:3
                                    for l=1:3
                                        K(i,k)= K(i,k)+ dNdX(A,j)*C(i,j,k,l)*dNdX_t(l,B);
                                    end
                                end
                            end
                        end



                        %M=nodalMass(N,rho,delta,A,B);
                        kL(3*(A-1)+1:3*A,3*(B-1)+1:3*B)=K;
                        %mL(3*(A-1)+1:3*A,3*(B-1)+1:3*B)=M;
                    end
                end

                %%% Mass Matrix                

                for A=1:nNoEl
                    for B=1:nNoEl

                        %K=nodalStiffness(dNdX,C,A,B);
                        % M=nodalMass(N,rho,delta,A,B);

                        M=zeros(3,3);
                        for i=1:3
                            for k=1:3
                                M(i,k)=N(A)*rho*delta(i,k)*N(B);
                            end
                        end

                       % kL(3*(A-1)+1:3*A,3*(B-1)+1:3*B)=K;
                        mL(3*(A-1)+1:3*A,3*(B-1)+1:3*B)=M;
                    end
                end
                
                kLocal=kLocal+kL*det(Jac)*GaussQuadMat(wqX,2)*GaussQuadMat(wqY,2)*GaussQuadMat(wqZ,2);
                mLocal=mLocal + mL*det(Jac)*GaussQuadMat(wqX,2)*GaussQuadMat(wqY,2)*GaussQuadMat(wqZ,2);
                fLocal=fLocal + 0; 
            end
        end
    end

% Global Matrices

    for i=1:nNoEl
        for j=1:nNoEl
            I=3*(connArr(e,i)-1)+1:3*connArr(e,i);
            J=3*(connArr(e,j)-1)+1:3*connArr(e,j);
            kGlobal(I,J)=kGlobal(I,J) + kLocal(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
            mGlobal(I,J)=mGlobal(I,J) + mLocal(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
        end
    end
    %     fGlobal(connArr(e,:),1)=fGlobal(connArr(e,:),1)+fLocal

end


%% Boundary Conditions

%%% Dirichlet displacement at z=L

Dirich_Disp=0.05; %m


j=1;
for i=1:nNodes
    % no disp at z=0
    if(NodesXYZ(i,3)==0)     
        dirich(j,1:2)=[3*(i-1)+1,0];
        dirich(j+1,1:2)=[3*(i-1)+2,0];
        dirich(j+2,1:2)=[3*(i-1)+3,0];
        j=j+3;
    end

    % 0.05m disp at z=0
    if(NodesXYZ(i,3)==L) 
        dirich(j,1:2)=[3*(i-1)+2,Dirich_Disp];
        j=j+1;
    end
end

nDirich=length(dirich);

%%% Initial conditions disp and velocity

u0=zeros(nDof*nNodes,1);  
v0=zeros(nDof*nNodes,1); 


%% Non Dirichlet Matrices formation

kNonDirichlet=zeros(nDof*nNodes-nDirich,nDof*nNodes-nDirich);
mNonDirichlet=zeros(nDof*nNodes-nDirich,nDof*nNodes-nDirich);
fNonDirichlet=zeros(nDof*nNodes-nDirich,1);


nonDirichnodes=zeros(nDof*nNodes-nDirich,1);
j=1;

for i=1:nDof*nNodes
    check=0;
    for k=1:nDirich
        
        if(i==dirich(k,1))
            check=1;
        end
        
    end
    
    if(check==0)
        nonDirichnodes(j)=i;
        j=j+1;
    end
end

%%% Non Dirichlet Matrices

for i=1:length(nonDirichnodes)
    for j=1:length(nonDirichnodes)
        kNonDirichlet(i,j)=kGlobal(nonDirichnodes(i),nonDirichnodes(j));
        mNonDirichlet(i,j)=mGlobal(nonDirichnodes(i),nonDirichnodes(j));
    end
end

f_dash=zeros(length(nonDirichnodes),1);
for i=1:length(nonDirichnodes)
    for j=1:nDirich
        f_dash(i)=kGlobal(nonDirichnodes(i),dirich(j,1))*dirich(j,2) + f_dash(i);
    end
end

fNonDirichlet=fGlobal(nonDirichnodes)-f_dash;


%% PART 1: Elasto-static Problem

% Nodal Displacement Solution

%kNonDirichlet=sparse(kNonDirichlet);
d=kNonDirichlet\fNonDirichlet;

dGlobal=zeros(nDof*nNodes,1);

dGlobal(dirich(:,1))=dirich(:,2);

dGlobal(nonDirichnodes)=d;

    for i=1:nNodes
        x(i)=dGlobal((i-1)*nDof+1);
        y(i)=dGlobal((i-1)*nDof+2);
        z(i)=dGlobal((i-1)*nDof+3);
    end
    
    u=[x',  y', z'];
    
    fem_to_vtk_Vector ('HW1_Part1_MH', NodesXYZ, connArr, u);
    

%% PART 2: Elasto-dynamic Problem

% for Transient Conditions

% damping Parameters

a = 1; 
b = 0.001;
cDirichlet=(a*mNonDirichlet + b*kNonDirichlet);    

beta=0.255;
gamma=2*beta;

% Time step parameters
tDelta=1000;   %% deltaT
tSteps=10;    %% number of time steps

%for files = 1:tSteps

    % Initializing disp, velovity
    dGlobal(:,1)=u0;
    vGlobal(:,1)=v0;

    % Using Newmark's Scheme

    d(:,1)=dGlobal(nonDirichnodes);
    v(:,1)=vGlobal(nonDirichnodes);
    acc(:,1)=mNonDirichlet\(fNonDirichlet - kNonDirichlet*d(:,1) - cDirichlet*v(:,1) );

    for i=1:tSteps
        d_tilde= d(:,i) + tDelta*v(:,i) + (1-2*beta)/2*tDelta^2*acc(:,i);
        v_tilde= v(:,i) + (1-gamma)*tDelta*acc(:,i);

        acc(:,i+1)= (mNonDirichlet+gamma*tDelta*cDirichlet+beta*tDelta^2*kNonDirichlet)\(fNonDirichlet - kNonDirichlet*d_tilde - cDirichlet*v_tilde);

        d(:,i+1)= d_tilde + beta*tDelta^2*acc(:,i+1);
        v(:,i+1)= v_tilde + gamma*tDelta*acc(:,i+1);
    end


    %%% plotting steps
    for files = 0:tSteps
        %files = i;
        tStepPlot= files+1;
        dGlobal(dirich(:,1))=dirich(:,2);
        dGlobal(nonDirichnodes)=d(:,tStepPlot);

        for ii=1:nNodes
            x(ii)=dGlobal((ii-1)*nDof+1);
            y(ii)=dGlobal((ii-1)*nDof+2);
            z(ii)=dGlobal((ii-1)*nDof+3);
        end

        u=[x',  y', z'];

        filename = strcat('HW1_Part2_MH_frame_',num2str(files));

        fem_to_vtk_Vector (filename, NodesXYZ, connArr, u);
        %i = files;
    end
%% Functions

%    % Shape Function
%  function [N, dNdz,dNde,dNdg]=shapeFunc(zhi,eta,gamma)
% 
%     N(1)=(1-zhi)*(1-eta)*(1-gamma)/8;
%     N(2)=(1+zhi)*(1-eta)*(1-gamma)/8;
%     N(3)=(1+zhi)*(1+eta)*(1-gamma)/8;
%     N(4)=(1-zhi)*(1+eta)*(1-gamma)/8;
%     N(5)=(1-zhi)*(1-eta)*(1+gamma)/8;
%     N(6)=(1+zhi)*(1-eta)*(1+gamma)/8;
%     N(7)=(1+zhi)*(1+eta)*(1+gamma)/8;
%     N(8)=(1-zhi)*(1+eta)*(1+gamma)/8;
%     
%     dNdz(1)=-(1-eta)*(1-gamma)/8;
%     dNdz(2)=(1-eta)*(1-gamma)/8;
%     dNdz(3)=(1+eta)*(1-gamma)/8;
%     dNdz(4)=-(1+eta)*(1-gamma)/8;
%     dNdz(5)=-(1-eta)*(1+gamma)/8;
%     dNdz(6)=(1-eta)*(1+gamma)/8;
%     dNdz(7)=(1+eta)*(1+gamma)/8;
%     dNdz(8)=-(1+eta)*(1+gamma)/8;
%     
%     
%     dNde(1)=-(1-zhi)*(1-gamma)/8;
%     dNde(2)=-(1+zhi)*(1-gamma)/8;
%     dNde(3)=(1+zhi)*(1-gamma)/8;
%     dNde(4)=(1-zhi)*(1-gamma)/8;
%     dNde(5)=-(1-zhi)*(1+gamma)/8;
%     dNde(6)=-(1+zhi)*(1+gamma)/8;
%     dNde(7)=(1+zhi)*(1+gamma)/8;
%     dNde(8)=(1-zhi)*(1+gamma)/8;
%     
%     dNdg(1)=-(1-zhi)*(1-eta)/8;
%     dNdg(2)=-(1+zhi)*(1-eta)/8;
%     dNdg(3)=-(1+zhi)*(1+eta)/8;
%     dNdg(4)=-(1-zhi)*(1+eta)/8;
%     dNdg(5)=(1-zhi)*(1-eta)/8;
%     dNdg(6)=(1+zhi)*(1-eta)/8;
%     dNdg(7)=(1+zhi)*(1+eta)/8;
%     dNdg(8)=(1-zhi)*(1+eta)/8;
%  end

 %%% Nodal Stiffness

% function[K]=nodalStiffness(dNdX,C,A,B)
% 
%  K=zeros(3,3);
%  dNdX_t=dNdX';
% 
%  for i=1:3
%      for j=1:3
%          for k=1:3
%              for l=1:3
%                  K(i,k)= K(i,k)+ dNdX(A,j)*C(i,j,k,l)*dNdX_t(l,B);
%              end
%          end
%      end
%  end
% end
% 
% %%% Nodal Mass
% function [M]=nodalMass(N,rho,delta,A,B)
%  M=zeros(3,3);
%  for i=1:3
%      for k=1:3
%          M(i,k)=N(A)*rho*delta(i,k)*N(B);
%      end
%  end
% 
% end

