%% ME964: HW 3
% Mahmudul Hassan
% Spring 2023

%% clear everything
clear all;
clc;
close all;

%% %% Variables

W=0.01; %m  % width
H=0.01; %m  % height

nu=0.499;   % poisson's ratio

mu = 1.83e10; % Pa
lambda=9.15e12; % Pa

delta=eye(2);
%% Create Mesh

% number of division in X Y
nElx=4;
nEly=4;

nNoEl=4;     %% 4 nodes in a Linear Quadrilateral elements with 1x1 point quadrature.

nDof=2;        %% Degrees or Freedom in 2D

nNodes = (nElx+1)*(nEly+1);

nEl = nElx*nEly; 

% Node coordinates initiation
NodesXYZ=zeros(nNodes,2);

elemH=H/nEly;
elemW=W/nElx;


for j=1:nEly+1
    for k=1:nElx+1
        NodesXYZ(k+(j-1)*(nElx+1),1)=(k-1)*elemW;  %x-coordinate
        NodesXYZ(k+(j-1)*(nElx+1),2)=(j-1)*elemH;  %y-coordinate
    end
end

%% Part 1(b)
% 
% nEl = 16;
% nElx=4;
% nEly=4;

% Connectivity arrays
connArr =zeros(nEl,4);

k=1;
step = 0;
for r=1:nEly
     
    for c=1:nElx
        nd1 = r + c-1 +step;
        nd2 = nd1+1;
        nd4 = nd2 +nElx;
        nd3 = nd4+1;

    connArr(k,:) = [nd1 nd2 nd3 nd4];
    k = k+1;
    end

    step = step + nElx;

end

%% Gauss Quadrature

gmatn = 2;

zz = 1/sqrt(3);
ee = 1/sqrt(3);
GaussQuadMat=[-zz -ee;
               zz -ee;
              -zz ee;
               zz ee];

gsswt = 1;

[gaussLen,non] = size(GaussQuadMat);
%% C_ijkl


C=zeros(2,2,2,2);

for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                C(i,j,k,l)=lambda*delta(i,j)*delta(k,l)  + 2*mu*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
            end
        end
    end
end

%%  Local and Global matrix formation

kGlobal=zeros(nDof*nNodes,nDof*nNodes);

fGlobal=zeros(nDof*nNodes,1);

%%% Local Matrices

for e=1:nEl
    kLocal=zeros(nDof*nNoEl,nDof*nNoEl);
    fLocal=zeros(nDof*nNoEl,1);

    for gss = 1:gaussLen
    % shape funcs
    N(1)=(1-GaussQuadMat(gss,1))*(1-GaussQuadMat(gss,2))/4;
    N(2)=(1+GaussQuadMat(gss,1))*(1-GaussQuadMat(gss,2))/4;
    N(3)=(1+GaussQuadMat(gss,1))*(1+GaussQuadMat(gss,2))/4;
    N(4)=(1-GaussQuadMat(gss,1))*(1+GaussQuadMat(gss,2))/4;

    % shape derivatives
    dNdz(1)=-(1-GaussQuadMat(gss,2))/4;
    dNdz(2)= (1-GaussQuadMat(gss,2))/4;
    dNdz(3)= (1+GaussQuadMat(gss,2))/4;
    dNdz(4)=-(1+GaussQuadMat(gss,2))/4;


    dNde(1)=-(1-GaussQuadMat(gss,1))/4;
    dNde(2)=-(1+GaussQuadMat(gss,1))/4;
    dNde(3)= (1+GaussQuadMat(gss,1))/4;
    dNde(4)= (1-GaussQuadMat(gss,1))/4;

    %jacobian

    

    Jac=[dNdz*NodesXYZ(connArr(e,:),1)  dNde*NodesXYZ(connArr(e,:),1) ;
        dNdz*NodesXYZ(connArr(e,:),2)   dNde*NodesXYZ(connArr(e,:),2) ];
    
    dNdX=[dNdz' dNde']/Jac;

    % stiffness
    for A=1:nNoEl
        for B=1:nNoEl


            K=zeros(2,2);
            dNdX_t=dNdX';

            for i=1:2
                for j=1:2
                    for k=1:2
                        for l=1:2
                            K(i,k)= K(i,k)+ dNdX(A,j)*C(i,j,k,l)*dNdX_t(l,B);
                        end
                    end
                end
            end


            kL(2*(A-1)+1:2*A,2*(B-1)+1:2*B)=K;

        end
    end

    kLocal=kLocal+kL*det(Jac)*gsswt;
    fLocal=fLocal + 0;
    end
%     size(kLocal)
    for i=1:nNoEl
        for j=1:nNoEl
            I=2*(connArr(e,i)-1)+1:2*connArr(e,i);
            J=2*(connArr(e,j)-1)+1:2*connArr(e,j);
            kGlobal(I,J)=kGlobal(I,J) + kLocal(2*(i-1)+1:2*i,2*(j-1)+1:2*j);
        end
    end
end
% size(kGlobal)
%% Boundary Conditions

%%% Dirichlet displacement at x,y = 1,1

Dirich_Disp=0.001; %m
angl = 45;
slpe = cosd(angl);  
Dirich_Disp_X = slpe*Dirich_Disp;
Dirich_Disp_Y = slpe*Dirich_Disp;

j=1;
for i=1:nNodes
    % no disp at y=0
    if(NodesXYZ(i,2)==0 )
        %dirich(j,1:2)=[j,0];
        %dirich(j,1:2)=[2*(i-1)+1,0];
        %dirich(j+1,1:2)=[2*(i-1)+2,0];

        dirich(j,1:2)=[2*i-1,0];
        dirich(j+1,1:2)=[2*i,0];
        
        j=j+2;
        
     end

    % no disp at x=0
    if(NodesXYZ(i,1)==0)
        dirich(j,1:2)=[2*i-1,0];
        dirich(j+1,1:2)=[2*i,0];

        j=j+2;

    end

        
    if(NodesXYZ(i,:)==[W,H])
        dirich(j,1:2)=[2*i-1,Dirich_Disp_X];
        dirich(j+1,1:2)=[2*i,Dirich_Disp_Y];
        j=j+2;
    end
end

nDirich=length(dirich);

%%% Initial conditions disp and velocity

u0=zeros(nDof*nNodes,1);

%% Non Dirichlet Matrices formation

kNonDirichlet=zeros(nDof*nNodes-nDirich,nDof*nNodes-nDirich);
fNonDirichlet=zeros(nDof*nNodes-nDirich,1);

nonDirichnodes=zeros(nDof*nNodes-nDirich,1);
j=1;

% Seperate Non-Dirichlet Node Numbers
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
    end
end

f_dash=zeros(length(nonDirichnodes),1);
for i=1:length(nonDirichnodes)
    for j=1:nDirich
        f_dash(i)=kGlobal(nonDirichnodes(i),dirich(j,1))*dirich(j,2) + f_dash(i);
    end
end

fNonDirichlet=fGlobal(nonDirichnodes)-f_dash;

%% Nodal Displacement Solutions

d=kNonDirichlet\fNonDirichlet;

dGlobal=zeros(nDof*nNodes,1);

dGlobal(dirich(:,1))=dirich(:,2);

dGlobal(nonDirichnodes)=d;

for i=1:nNodes
    x(i)=dGlobal((i-1)*nDof+1);
    y(i)=dGlobal((i-1)*nDof+2);

end

u=[x',  y'];

fem_to_vtk_Vector ('HW3_MH_Part1-b_fig-1c', NodesXYZ, connArr, u);




