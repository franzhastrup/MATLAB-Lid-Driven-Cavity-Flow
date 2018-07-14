%% << Function for computation of H1,H2 terms >>

 function [H1,H2] = NS2dHfunctions(n,dx,Re,u,v)
%% *****************************************************************************
% NS2dHfunctions:
%   Computes H1- and H2-functions for Navier-Stokes equations in 2D for 
%   simulation of translating lid flow in square cavity described by the
%   incompressible unsteady 2D Navier-Stokes equations in nondimensional form.
%   Finite Volume Method on staggered cartesian mesh with uniform square cells. 
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator). 
%   H1- and H2-expressions:
%        H1 = 1/dx*(1/Re*dudx-u*u)|_w^e + 1/dy*(1/Re*dudy-u*v)|_s^n
%        H2 = 1/dx*(1/Re*dvdx-u*v)|_w^e + 1/dy*(1/Re*dvdy-v*v)|_s^n
%
% Syntax:
%   [H1,H2] = NS2dHfunctions(n,dx,Re,u,v)
%
% Input  :
%   n    :  Number of cells along x,y-axis
%   dx   :  Cell size in x,y
%   Re   :  Global Reynolds number number
%   u    :  Horizontal velocity component array, u(1:n+2,1:n+1)
%   v    :  Vertical velocity component array, v(1:n+1,1:n+2)
%
% Output :
%   H1   :  H1-function for horizontal momentum equation, H1(1:n,1:n+1) 
%   H2   :  H2-function for vertical momentum equation,  H2(1:n+1,1:n)
%
%   Author      : Franz Hastrup-Nielsen
%   Date        : 30-10-2015
%*******************************************************************************

%% Compute H1,H2 functions
% << compute H1,H2 terms from notes p. 135 eq. (6.20)+(6.21) >>

%% H1
uP = (u(2:end-1,2:end)+u(2:end-1,1:end-1))/2; %Getting u on P grid
uw = uP(:,1:end-1);   %West
ue = uP(:,2:end);   %East
uFD = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; %Getting u on FD grid
us = uFD(1:end-1,:);
un = uFD(2:end,:);

dudx = (u(2:end-1,2:end)-u(2:end-1,1:end-1))/dx; %on P-grid
dudxe = dudx(:,2:end);
dudxw = dudx(:,1:end-1);

dudy =(u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx; %on FD-grid
dudy(1,:) = (u(2,2:end-1)-u(1,2:end-1))*2/dx;
dudy(end,:) =(u(end,2:end-1)-u(end-1,2:end-1))*2/dx;

dudys = dudy(1:end-1,:);
dudyn = dudy(2:end,:);

v_us = (v(1:end-1,3:end-1)+v(1:end-1,2:end-2))/2;
v_un = (v(2:end,3:end-1)+v(2:end,2:end-2))/2;

H1 = zeros(size(u)-[2 0]);
H1(:,2:end-1) = 1/dx*(1/Re*dudxe-ue.*ue)-1/dx*(1/Re*dudxw-uw.*uw)+(1/dx*(1/Re*dudyn-un.*v_un)-1/dx*(1/Re*dudys-us.*v_us));

%% H2 
vFD = (v(2:end-1,2:end)+v(2:end-1,1:end-1))/2;
vFD(:,1) = v(2:end-1,1); 
vFD(:,end) = v(2:end-1,end); 
vw = vFD(:,1:end-1);
ve = vFD(:,2:end);

vP = (v(1:end-1,2:end-1)+v(2:end,2:end-1))/2;
vn = vP(2:end,:);
vs = vP(1:end-1,:);

dvdx = (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dx;
dvdx(:,1) = (v(2:end-1,2)-v(2:end-1,1))*2/dx;
dvdx(:,end) = (v(2:end-1,end)-v(2:end-1,end-1))*2/dx;
dvdxe = dvdx(:,2:end);
dvdxw = dvdx(:,1:end-1);

dvdy = (v(2:end,2:end-1)-v(1:end-1,2:end-1))/dx;
dvdyn = dvdy(2:end,:);
dvdys = dvdy(1:end-1,:);

u_ve = (u(2:end-2,2:end)+u(3:end-1,2:end))/2;
u_vw = (u(2:end-2,1:end-1)+u(3:end-1,1:end-1))/2;
H2 = zeros(size(v) - [0 2]);
H2(2:end-1,:) = 1/dx*(1/Re*dvdxe-u_ve.*ve)-1/dx*(1/Re*dvdxw-u_vw.*vw)+1/dx*(1/Re*dvdyn-vn.*vn)-1/dx*(1/Re*dvdys-vs.*vs);
