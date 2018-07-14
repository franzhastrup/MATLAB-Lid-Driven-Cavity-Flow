%% *****************************************************************************
% File: NavierStokes2dMovingLidSquareCavityFlow.m
%   Script for simulation of translating lid flow in square cavity described by 
%   incompressible unsteady 2D Navier-Stokes equations in nondimensional form.
%   Finite Volume Method on staggered cartesian mesh with uniform square cells. 
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator). 
%   NS equations integrated in time using the explicit Euler method.
%   SIMPLE scheme used to derive Poisson equation for pressure. Solved directly 
%   by LU-factorizing Laplacian matrix before time-stepping, then forward- and 
%   backward-substituting using Matlab's backslash operator in each time step.  
%
% Input        :
%   n          :  Number of cells along x,y-axis
%   Re         :  Global Reynolds number 
%   Ulid       :  Lid velocity (+-1, set =-1 to agree Bruneau & Saad 2006)
%   maxstep    :  Maximum number of time steps for explicit Euler method.
%   dt         :  Time step for explicit Euler method.
%   steadytol  :  Tolerance on maximum rate of change for termination at steady 
%                 solution, thus dUmax/dt < steadytol, dUmax = max(dumax,dvmax)   
%   statstride :  stride for writting status to screen
%
% Output       :
%   u          :  Horizontal velocity component array, u(1:n+2,1:n+1)
%   v          :  Vertical velocity component array, v(1:n+1,1:n+2)
%   p          :  Pressure field array, p(1:n,1:n)
%   cmc        :  Control volume mass conservation array, cmc(1:n,1:n)
%   gmchist    :  Global mass conservation history, gmchist(tstep)=sum(cmc(:))                   
%   cmchist    :  Max. cell mass conser. history, cmchist(step)=norm(cmc(:),inf)
%   A          :  Laplacian system matrix, A(1:n^2,1:n^2)
%   s          :  Source array for Poisson equation for pressure, s(1:n,1:n)
%   Xu,Yu      :  Coord. arrays for u-grid, Xu(1:n+2,1:n+1), Yu(1:n+2,1:n+1)
%   Xv,Yv      :  Coord. arrays for v-grid, Xv(1:n+1,1:n+2), Yv(1:n+1,1:n+2)
%   Xp,Yp      :  Coord. arrays for p-grid, Xp(1:n,1:n), Yp(1:n,1:n)
%   Xi,Yi,     :  Coord. arrays for fd-grid, Xi(1:n+1,1:n+1), Yi(1:n+1,1:n+1)
%   Plots of   :  i) velocity field streamlines, ii) pressure field contours. 
%                 iii) gmchist and steadyhist
%   Results written to file 'NS2d_Re<>_n<>.mat'
%
%   Author      : Franz Hastrup-Nielsen
%   Version     : 1.1
%*******************************************************************************

%% Clear
clear all
close all
clc

%% Input
n          = 161;           % number of cells along x,y-axis
Re         = 1000;          % global Reynolds number
Ulid       = -1;            % lid velocity (+-1)
maxstep    = 10000000;      % maximum number of explicit Euler time steps

%dt         = 0.03;         % time step -- base it on stability analysis
steadytol  = 1e-5;          % tolerance on dUmax/dt to reach steady solution
statstride = 10000;         % stride for writting status to screen
Aquick = 'lu-fac';          %if 'direct' A is used as direct matrix, 
                            %if 'lu-fac' lu-factorization is used
Pois_corr = 'correction';   %if 'no-correction' no pressure correction in 
                   %Poisson, if 'correction' pressure correction in Poisson

%% Grid arrays
% Assemble staggered grid coordinate arrays
profile on;
[dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi] = StaggeredMesh2dSquare(n); 
U = 1;
dtstab     = min([dx^2*Re/4 1/(U^2*Re)]);          %Time step based on stab. analysis, Stability PDF, pg. 25
dt         = dtstab;
% Set initial flow field such that du/dx+dv/dy = 0
u = zeros(n+2,n+1);             % initial horizontal velocity component
v = zeros(n+1,n+2);             % initial vertical velocity component
u(n+2,:) = Ulid;                % set lid velocity in u-array

% Various other arrays
P          = zeros(n);          % allocate pressure array s.t. P(:) = A\s(:);
gmchist    = zeros(1,maxstep);  % allocate global mass conservation hist. vector
cmchist    = zeros(1,maxstep);  % allocate max. cell mass conserva. hist. vector
steadyhist = ones(1,maxstep);   % allocate global steady solution hist. vector

%% Assemble and factorize Laplacian operator matrix
% << Use your own function NS2dLaplaceMatrix to assemble Laplacian matrix >>
% << Once model is working, do re-ordering + lu-factorization (use colamd+lu) >>

[A] = NS2dLaplaceMatrix(n);
if  strcmp(Aquick,'lu-fac')
perm = colamd(A);
[L,U] = lu(A(perm,perm));
end
%% Time integrate NS-equations to reach steady solution
% << Make a for/while loop over the number of time steps. >>
% << In each time step, do the following sequence of computations: >>
% <<  - compute H1,H2 terms with your own function NS2dHfunctions >>
% <<  - compute source array for Poisson equation incl. hom Neumann BC >>
% <<  - solve Poisson equation for pressure - fixate it to avoid fluctuations.>>
% <<  - compute velocity time derivatives, notes eq. (6.20)+(6.21) >>
% <<  - step velocity field forward using explicit Euler method >>
% <<  - compute mass conservation on cell level, and from this on global level>>
% <<    store measure of global and cell mass conserv. in gmchist and cmchist >>
% <<  - compute the maximum rate of change (deviation from steady solution) >>
% <<    from the velocity time derivatives, and store it in steadyhist >>
% <<  - write some status message to screen if mod(step,statstride)==0 >>
% <<  - interrupt time stepping if steadyhist(step) < steadytol (use break) >>
warning('off', 'MATLAB:nearlySingularMatrix')

step = 0;
residual = 1;
tflag = tic;               % start timing using Matlab's tic-toc function
while step < maxstep
    step = step +1;
    [H1,H2] = NS2dHfunctions(n,dx,Re,u,v);
    if strcmp(Pois_corr,'no-correction' )  
    s = ((H1(:,2:end)-H1(:,1:end-1))*dx+(H2(2:end,:)-H2(1:end-1,:))*dx);
    elseif strcmp(Pois_corr,'correction')  
        du = (u(2:end-1,2:end)-u(2:end-1,1:end-1)); %on P-grid
        dv = (v(2:end,2:end-1)-v(1:end-1,2:end-1));
    s = ((H1(:,2:end)-H1(:,1:end-1))*dx+(H2(2:end,:)-H2(1:end-1,:))*dx) + dx/dt*(du+dv);  
    end
    if strcmp( Aquick,'direct') 
        P(:) = A\s(:);
    elseif strcmp( Aquick,'lu-fac') 
        s = s(:);
        P=P(:);
        P(perm,1) = U\(L\s(perm,1));
    end
    
    P = reshape(P,[n n]);
    P=(P-P(ceil(n/2),ceil(n/2)));
    
    dpdx = (P(:,2:end)-P(:,1:end-1))/dx; %eq. 6.20
    dpdy = (P(2:end,:)-P(1:end-1,:))/dx; %eq. 6.21
    udiv = H1(:,2:end-1)-dpdx;
    vdiv = H2(2:end-1,:)-dpdy;
    steadyhist(step) = (max([max(max(abs(udiv))) max(max(abs(vdiv)))]));
    
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1)+dt*udiv;
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1)+dt*vdiv;
    
    cmchist(step) = max(max(abs(u(2:end-1,2:end)-u(2:end-1,1:end-1)+v(2:end,2:end-1)-v(1:end-1,2:end-1))));
    gmchist(step) = sum(sum(u(2:end-1,2:end)-u(2:end-1,1:end-1)+v(2:end,2:end-1)-v(1:end-1,2:end-1)));
    tvec(step) = dt*step;
    if mod(step,statstride)==0
      msgbox(sprintf('residual = %-3.2e, time = %.2f [s]',steadyhist(step),tvec(step)))
    end 
    if steadyhist(step) < steadytol
     msgbox('Defined steady tolerance reached')     
     break
     tvec = tvec(1:step);
    end
end
time = toc(tflag)   % stop timing using Matlab's tic-toc function
% time_save = tvec(1:step);

%% Plots
% << produce various plots of the results, e.g streamlines, pressure field, >>
% << and gmchist, cmchist and steadyhist to show change of sol. with time >>

uP = (u(2:end-1,2:end)+u(2:end-1,1:end-1))/2; %Getting u on P grid
    vP = (v(1:end-1,2:end-1)+v(2:end,2:end-1))/2;
%save('161','Xp','Yp','P','Xu','Yu','uP','u','Xv','Yv','vP','v')
figure
streamslice(Xp,Yp,uP,vP)
axis tight
title(sprintf('Streamlines, Re = %.f, n = %i, dt = %-4.1e, time = %.2f s, Deviation from SST = %-4.1e',Re,n,dt,tvec(step),steadyhist(step)),'Interpreter','Latex','Fontsize',14)
axis([0 1 0 1])
ylabel('y','Interpreter','Latex','Fontsize',14)
xlabel('x','Interpreter','Latex','Fontsize',14)
daspect([1 1 1])

figure 
surf(Xp,Yp,P,'EdgeColor','none'),view(0,90), colorbar
% caxis([0 1])
% xlim([0 1])
% ylim([0 1])
%axis equal
title(sprintf('Pressure field, Re = %.f, n = %i, dt = %-4.1e, time = %.2f s, Deviation from SST = %-4.1e',Re,n,dt,tvec(step),steadyhist(step)),'Interpreter','Latex','Fontsize',14)
ylabel('y','Interpreter','Latex','Fontsize',14)
xlabel('x','Interpreter','Latex','Fontsize',14)

figure
semilogy(tvec(1:step),cmchist(1:step))
%title('mass error, cell')
xlim([0 tvec(step)])
title('Local continuity(mass conservation) error','Interpreter','Latex','Fontsize',14)
ylabel('Mass conservation error','Interpreter','Latex','Fontsize',14)
xlabel('time [s]','Interpreter','Latex','Fontsize',14)

figure
plot(tvec(1:step),gmchist(1:step))
%title('mass error, global')
xlim([0 tvec(step)])
title('Global continuity(mass conservation) error','Interpreter','Latex','Fontsize',14)
ylabel('Mass conservation error','Interpreter','Latex','Fontsize',14)
xlabel('time [s]','Interpreter','Latex','Fontsize',14)

figure
semilogy(tvec(1:step),steadyhist(1:step))
xlim([0 tvec(step)])
title('Deviation from steady state','Interpreter','Latex','Fontsize',14)
ylabel('Maximum rate of change','Interpreter','Latex','Fontsize',14)
xlabel('time [s]','Interpreter','Latex','Fontsize',14)

print('-f1','-depsc2',sprintf('Streamslice_n%i_Re%i.eps',n,Re))
%print('-f2','-depsc2',sprintf('Pressure_n%i_Re%i.eps',n,Re))
% print('-f3','-depsc2',sprintf('LocalMassError_n%i_Re%i.eps',n,Re))
%print('-f4','-depsc2',sprintf('GlobalMassError_n%i_Re%i.eps',n,Re))
%print('-f5','-depsc2',sprintf('Residuals_n%i_Re%i.eps',n,Re))