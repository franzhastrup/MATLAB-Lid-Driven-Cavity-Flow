%% << Function for computation of Laplacian matrix A >>
function [A] = NS2dLaplaceMatrix(n)
% *****************************************************************************
% NS2dLaplaceMatrix:
%   Assembles 2D Laplacian operator matrix with pure Neumann boundary conditions
%   for a structured cartesian, cell centered (finite volume) mesh with n square
%   cells in x,y of size dx. 2nd order CDS applied for the Laplace operator,
%   which is scaled with a factor $dx^2$, which thus should be multiplied on the
%   rhs-vector too.
% 
% Syntax:
%   [A] = NS2dLaplaceMatrix(n)
%
% Input  :
%   n    :  Number of cells along x,y-axis
%
% Output :
%   A    :  Laplacian operator matrix, A(1:n^2,1:n^2)
%
%   Author      : Franz Hastrup-Nielsen
%*******************************************************************************

%% Assemble Laplacian operator matrix
% << implement the Laplacian matrix with homogeneous Neumann boundary >>
% << conditions on all walls >>
    dx = 1/n;                       % cell size in x,y-direction
    D = 1; %Page 121, table 5.3
    aW = zeros(n)+D;
    aS = zeros(n)+D;
    aE = zeros(n)+D;
    aN = zeros(n)+D;
    aP = -(aW+aE+aN+aS); %Eq. 5.59
    
i=1; aP(:,i) = aP(:,i) + aW(:,i); aW(:,i)=0; %West
i=1; aP(i,:) = aP(i,:) + aS(i,:); aS(i,:)=0; %South
i=n; aP(:,i) = aP(:,i) + aE(:,i); aE(:,i)=0; %East
i=n; aP(i,:) = aP(i,:) + aN(i,:); aN(i,:)=0; %North

%% Assemble system matrix

aW = aW(:);
aN = aN(:);
aS = aS(:);
aE = aE(:);
aP = aP(:);

DiaPos = [-n -1 0 1 n];
Adiags = [[aW(n+1:end);zeros(n,1)],[aS(2:end);0],[aP],[0;aN(1:end-1)],[zeros(n,1);aE(1:end-n)]];     
A      = spdiags(Adiags,DiaPos,n^2,n^2);