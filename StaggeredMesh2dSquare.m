function [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi] = StaggeredMesh2dSquare(n)
%% *****************************************************************************
% File: StaggeredMesh2dSquare.m
%   Function for assembling of coordinate arrays for staggered mesh on 2D square
%   
% Syntax:
%   [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi] = StaggeredMesh2dSquare(n)
%
% Input        :
%   n          :  Number of cells along x,y-axis
%
% Output       :
%   dx         :  Cell size in x,y
%   Xu,Yu      :  Coord. arrays for u-grid,  Xu(1:n+2,1:n+1),  Yu(1:n+2,1:n+1)
%   Xv,Yv      :  Coord. arrays for v-grid,  Xv(1:n+1,1:n+2),  Yv(1:n+1,1:n+2)
%   Xp,Yp      :  Coord. arrays for p-grid,  Xp(1:n,1:n),      Yp(1:n,1:n)
%   Xi,Yi,     :  Coord. arrays for fd-grid, Xi(1:n+1,1:n+1),  Yi(1:n+1,1:n+1)
%   Plot of mesh arrangment for n < 12
%*******************************************************************************

%% Assemble mesh coordinate arrays
dx      = 1/n;                  % cell size in x,y
xf      = 0:dx:1;               % cell face coordinate vector, 1D
xc      = dx/2:dx:1-dx/2;       % cell center coordinate vector, 1D
xb      = [0 xc 1];             % cell center coordinate vector incl. boundaries
[Xu,Yu] = meshgrid(xf,xb);      % u-grid coordinate arrays
[Xv,Yv] = meshgrid(xb,xf);      % v-grid coordinate arrays
[Xp,Yp] = meshgrid(xc,xc);      % p-grid coordinate arrays
[Xi,Yi] = meshgrid(xf,xf);      % fd-grid coordinate arrays