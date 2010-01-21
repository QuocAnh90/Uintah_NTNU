%ICE ICE Algorithm for the Shock Tube Problem in 1D (main driver) - 03/2005
%
%   This script tests full ICE time-stepping algorithm for Sod's Shocktube
%   problem in one space dimension.
%   We solve the compressible Euler equations for density (rho), velocity
%   (u), internal energy (e) and pressure (p). 
%
%   Reference:  "ICE, explicit pressure, single material, reaction model"
%               by Todd Harman, 09/24/04.
%               ICE input file: rieman_sm.ups
%
%   See also SETBOUNDARYCONDITIONS, ADVECTRHO, ADVECTQ.


function [tfinal, x_CC, delX, rho_CC, xvel_CC, press_CC, temp_CC]=iceTotalEnergy(nCells)

close all;
globalParams;                                       % Load global parameters
setenv('LD_LIBRARY_PATH', ['/usr/lib']);
[LF] = LethuyFunctions;

%================ ICE Interal Parameters, Debugging Flags ================
% Debug flags
P.compareUintah     = 0;         % Compares vs. Uintah ICE and plots results
P.debugSteps        = 0;         % Debug printout of steps (tasks) within timestep
P.debugAdvectRho    = 0;         % Debug printouts in advectRho()
P.debugAdvectQ      = 0;         % Debug printouts in advectQ()
P.printRange        = [27:32];   % Range of indices to be printed out (around the shock front at the first timestep, for testing)
P.plotInitialData   = 0;         % Plots initial data
P.plotInterval      = 100;
P.plotResults       = 0;         % plot intermediate data

%______________________________________________________________________
%     Problem Setup
%================ Set Parameters into struct P ================
% Grid
P.nCells            = nCells;    % Number of cells in each direction
P.extraCells        = 1;         % Number of ghost cells in each direction
% Time-stepping
P.initTime          = 0.0;       % Initial simulation time [sec]
P.writeData         = 1;         % output the final timestep to a .dat file
P.delt_init         = 1e-20;     % First timestep [sec]
P.maxTimeSteps      = 2000       % Maximum number of timesteps [dimensionless]
P.CFL               = 0.25;      % Courant number (~velocity*delT/delX) [dimensionless]
P.advectionOrder    = 2;         % 1=1st-order advection operator; 2=possibly-limited-2nd-order
P.gamma             = 1.4;       % gamma coefficient in the Equation of State (EOS)

Problem             = 'Test1';
switch (Problem)
  case {'Test1','Test2','Test3','Test4','Lax'}
    P.boxLower     = 0;             
    P.boxUpper     = 1;
  case {'ShuOsher'}
    P.boxLower     = -5.0;             
    P.boxUpper     = 5.0;
end

%================ Grid Struct (G) ======= ================
G.nCells      = P.nCells;                           % # interior cells
G.delX        = (P.boxUpper-P.boxLower)./G.nCells;  % Cell length
G.ghost_Left  = 1;                                  % Index of left ghost cell
G.ghost_Right = G.nCells+2*P.extraCells;            % Index of right ghost cell
G.first_CC    = 2;                                  % Index of first interior cell
G.first_FC    = 2;                                  % Index of first interior xminus face
G.last_CC     = G.nCells+1;                         % Index of last interior cell
G.last_FC     = G.ghost_Right;                      % index of last xminus face

%______________________________________________________________________
%     Allocate arrays
%================ Cell Centered (CC) ================

totCells        = P.nCells + 2*P.extraCells;              % Array size for CC vars
x_CC            = zeros(G.ghost_Left,G.ghost_Right);      % Cell centers locations (x-component)
rho_CC          = zeros(G.ghost_Left,G.ghost_Right);      % Density rho
xvel_CC         = zeros(G.ghost_Left,G.ghost_Right);      % Velocity u (x-component)
E_CC            = zeros(G.ghost_Left,G.ghost_Right);      % TotalEnergy
temp_CC         = zeros(G.ghost_Left,G.ghost_Right);      % Temperature T
press_eq_CC     = zeros(G.ghost_Left,G.ghost_Right);      % equilibraton pressure
press_CC        = zeros(G.ghost_Left,G.ghost_Right);      % Pressure p
volfrac_CC      = zeros(G.ghost_Left,G.ghost_Right);      % Volume fraction theta, = 1 for this single-material problem
spvol_CC        = zeros(G.ghost_Left,G.ghost_Right);      % Specific volume v0 = theta/rho = 1/(specific density)
mass_CC         = zeros(G.ghost_Left,G.ghost_Right);      % Mass at time n+1 (see step 9)
delPDilatate    = zeros(G.ghost_Left,G.ghost_Right);      % Pressure increment due to dilatation

%================ Cell Centered (CC), Lagrangian values (L) ================

del_mom         = zeros(G.ghost_Left,G.ghost_Right);      % Momentum accumulated sources
del_eng         = zeros(G.ghost_Left,G.ghost_Right);      % Energy accumulated sources
mass_L          = zeros(G.ghost_Left,G.ghost_Right);      % Mass ( = rho * cell volume )
mom_L           = zeros(G.ghost_Left,G.ghost_Right);      % Momentum ( = mass * velocity )
eng_L           = zeros(G.ghost_Left,G.ghost_Right);      % Energy (= mass * internal energy = mass * cv * temperature)

%================ Face Centered (FC) ================
x_FC            = zeros(G.ghost_Left,G.last_FC);          % Face centers locations (x-component)  
xvel_FC         = zeros(G.ghost_Left,G.last_FC);          % Velocity u (x-component)              
press_FC        = zeros(G.ghost_Left,G.last_FC);          % Pressure p                            
speedSound_FC   = zeros(G.ghost_Left,G.last_FC);          % c^*                                   

%================ Node Centered (NC) ================
totNodes        = P.nCells + 1;                           % Array size for NC vars
mass_vrtx_1     = zeros(1,totNodes);                      % Mass at vertices (for advection)
mass_vrtx_2     = zeros(1,totNodes);                      % --------//-------  (for advection)

%______________________________________________________________________
%     Initialization

%================ Useful constants ================

d_SMALL_NUM = 1e-100;                           % A small number (for bullet-proofing
d_TINY_RHO  = 1.0e-12;
delT        = P.delt_init;                      % Init timestep

%================ Initialize interior cells ================
% Initial data at t=0 in the interior domain.
x_CC    = ([G.ghost_Left:G.ghost_Right]-1-0.5).*G.delX + P.boxLower; 
x_FC(1) = -G.delX + P.boxLower;

for j = G.first_FC:G.last_FC   % Loop over all xminus cell faces
  x_FC(j) = (j-G.first_FC).*G.delX + P.boxLower;  
end

for j = G.first_CC:G.last_CC
    switch (Problem)
        case 'Test1'                    
            P.maxTime = 0.2;      
            if (x_CC(j) < 0.3)
                rho_CC(j)    = (1.0+d_TINY_RHO);
                xvel_CC(j)   = 0.75;
                E_CC(j)      = 2.78125;
                press_CC(j)  = 1.0;   
            else
                rho_CC(j)    = (1.0+d_TINY_RHO)*0.125;
                xvel_CC(j)   = 0.0;
                E_CC(j)      = 2.0;
                press_CC(j)  = 0.1; 
            end
        case 'Test2'
            P.maxTime = 0.15;      
            if (x_CC(j) < 0.5)
                rho_CC(j)    = (1.0+d_TINY_RHO);
                xvel_CC(j)   = -2.0;
                E_CC(j)      = 3.0;
                press_CC(j)  = 0.4;   
            else
                rho_CC(j)    = (1.0+d_TINY_RHO);
                xvel_CC(j)   = 2.0;
                E_CC(j)      = 3.0;
                press_CC(j)  = 0.4; 
            end
        case 'Test3'
            P.maxTime = 0.011;      
            if (x_CC(j) < 0.5)
                rho_CC(j)    = (1.0+d_TINY_RHO);
                xvel_CC(j)   = 0.0;
                E_CC(j)      = 2500.0;
                press_CC(j)  = 1000.0;   
            else
                rho_CC(j)    = (1.0+d_TINY_RHO);
                xvel_CC(j)   = 0.0;
                E_CC(j)      = 0.25;
                press_CC(j)  = 0.01; 
            end
        case 'Test4'
            P.maxTime = 0.034;      
            if (x_CC(j) < 0.4)
                rho_CC(j)    = (1.0+d_TINY_RHO)*5.99924;
                xvel_CC(j)   = 19.5975;
                E_CC(j)      = 384.0944978343299;
                press_CC(j)  = 460.894;  
            else
                rho_CC(j)    = (1.0+d_TINY_RHO)*5.99242;
                xvel_CC(j)   = -6.19633;
                E_CC(j)      = 38.405929430605944;
                press_CC(j)  = 46.0950; 
            end
        case 'ShuOsher'
            P.maxTime = 1.8;      
            if (x_CC(j) < -4.0)
                rho_CC(j)    = (1.0+d_TINY_RHO)*3.85714;
                xvel_CC(j)   = 2.62936;                
                press_CC(j)  = 10.33333; 
                E_CC(j)      = press_CC(j)/(0.4*rho_CC(j)) + 0.5*xvel_CC(j)*xvel_CC(j);
            else
                rho_CC(j)    = (1.0+d_TINY_RHO)*(1.0+0.2*sin(5*x_CC(j)));
                xvel_CC(j)   = 0.0;
                press_CC(j)  = 1.0; 
                E_CC(j)      = press_CC(j)/(0.4*rho_CC(j)) + 0.5*xvel_CC(j)*xvel_CC(j);
                
            end
        case 'Lax'
            P.maxTime = 0.16;  
            if (x_CC(j) < 0.5)
                rho_CC(j)    = (1.0+d_TINY_RHO)*0.445;
                xvel_CC(j)   = 0.698;                
                press_CC(j)  = 3.528; 
                E_CC(j)      = 19.5766;
            else
                rho_CC(j)    = (1.0+d_TINY_RHO)*0.5;
                xvel_CC(j)   = 0.0;
                press_CC(j)  = 0.571; 
                E_CC(j)      = 2.855;                
            end
    end
end

%================ Initialize ghost cells ================
% Impose boundary conditions (determine ghost cell values from
% interior cell values).
rho_CC      = setBoundaryConditions(rho_CC  ,'rho_CC',   G);
xvel_CC     = setBoundaryConditions(xvel_CC ,'xvel_CC',  G);
E_CC        = setBoundaryConditions(E_CC ,   'e_CC',     G);
press_CC    = setBoundaryConditions(press_CC,'press_CC', G);

%================ Initialize graphics ================
if (P.plotInitialData)
  figure(1);
  set(gcf,'position',[100,1000,1000,400]);
  %================ Plot results ================

  subplot(2,2,1), plot(rho_CC);
  %xlim([P.boxLower(1) P.boxUpper(1)]);
  legend('\rho');
  grid on;

  subplot(2,2,2), plot(x_CC,xvel_CC);
  xlim([P.boxLower(1) P.boxUpper(1)]);
  legend('u_1');
  grid on;

  subplot(2,2,3), plot(x_CC,E_CC);
  xlim([P.boxLower(1) P.boxUpper(1)]);
  legend('E_CC');
  grid on;

  subplot(2,2,4), plot(x_CC,press_CC);
  xlim([P.boxLower(1) P.boxUpper(1)]);
  legend('p');
  grid on;

  %M(tstep) = getframe(gcf);
  print -depsc iceInitialData.eps
  pause
end

%______________________________________________________________________
%     Time integration loop
t = P.initTime + delT;
tstep = 0;
for tstep = 1:P.maxTimeSteps
  fprintf('\n_____________________________tstep=%d, t=%e, prev. delT=%e\n', tstep, t, delT);
  
  %_____________________________________________________
  % 0. Dummy setting for a single-material problem
  % Set the volume fraction and specific volume.

  if (P.debugSteps)
    fprintf('Step 0: dummy setting for single-material\n');
  end

  volfrac_CC = ones(G.ghost_Left,G.ghost_Right);             % Single material ==> covers 100% of each cell (volfrac=1)
  spvol_CC   = volfrac_CC ./ (rho_CC + d_SMALL_NUM);         % Specific volume, here also = 1/rho

  %_____________________________________________________
  % 1. Compute thremodynamic/transport properties
  % These are constants in this application and were already computed in the
  % initialization stage.

  %_____________________________________________________
  % 2. Compute the equilibration pressure
  % Using an equation of State for ideal gas.
  if (P.debugSteps)
    fprintf('Step 2: compute equilibration pressure\n');
  end

  % Update pressure from EOS
  press_eq_CC      = (P.gamma-1.0).*rho_CC.*(E_CC-0.5.*xvel_CC.*xvel_CC);      % Evaluate p from EOS

  % Compute speed of sound at cell centers
  %DpDrho          = (P.gamma-1.0).*(E_CC-0.5.*xvel_CC.*xvel_CC);         % d P / d rho
  %DpDe            = (P.gamma-1.0).*rho_CC;       % d P / d e
  %tmp             = DpDrho + ( DpDe .* (press_eq_CC./rho_CC.^2));
  
  tmp             = P.gamma*press_eq_CC./ rho_CC;
  speedSound_CC   = sqrt(tmp);                   % Speed of sound


  % Set boundary conditions on p
  press_eq_CC     =  setBoundaryConditions(press_eq_CC,'press_CC',G);

  %_____________________________________________________
  % 3. Compute sources of energy due to chemical reactions
  % Not applicable to this model.


  %_____________________________________________________
  % 4. Compute the face-centered velocities
  if (P.debugSteps)
    fprintf('Step 4: compute face-centered velocities\n');
  end

if(0)
  for j = G.first_FC:G.last_FC   % Loop over all xminus cell faces
    L = j-1;
    R = j;
    term1   = (rho_CC(L)*xvel_CC(L) + rho_CC(R)*xvel_CC(R)) / (rho_CC(L) + rho_CC(R));
    term2a  = (2.0*spvol_CC(L)*spvol_CC(R)) / (spvol_CC(L) + spvol_CC(R));
    term2b  = (press_eq_CC(R) - press_eq_CC(L))/G.delX;    
    xvel_FC(j)  = term1 - delT*(term2a*term2b);
  end
end

if(1)  
  xvel_FC = LF.computeVel_FC(rho_CC, xvel_CC, spvol_CC, E_CC, press_eq_CC, delT, P, G);
end
  % Set boundary conditions on u @ face centers
  xvel_FC = setBoundaryConditions(xvel_FC,'xvel_FC', G);

  %_____________________________________________________
  % 5. Compute delta p
  % p satisfies a differential equation that can be
  % derived from the other ones plus the EOS. It looks
  % something like: p_t = -div(theta*u). Here we advance in time
  % p from tn to t_{n+1} using the advection operator of
  % theta*u.
  if (P.debugSteps)
    fprintf('Step 5: compute delta p\n');
  end

  % Compute the advection operator of theta*u
  [ofs, rx] = outFluxVol(xvel_CC, xvel_FC, delT, G);                  % Compute the outflux volumes ofs and centroid vector rx. Needs to be done once for all calls of advectQ in this timestep.
  
  [q_advected, gradLim, grad_x, mass_slab, mass_vrtx_1, mass_vrtx_2] = ...
    advectRho(volfrac_CC, ofs, rx, xvel_FC, G);                       % Treat vol frac (theta) advection like a mass advection (van Leer limiter, etc.)

  for j = G.first_CC:G.last_CC
    term1 = ((speedSound_CC(j)^2) / spvol_CC(j) );
    delPDilatate(j) = term1 * q_advected(j);
    press_CC(j) = press_eq_CC(j) + delPDilatate(j);                   % Advance pressure in time
  end
  press_CC =  setBoundaryConditions(press_CC,'press_CC',G);           % Set boundary conditions on the new pressure


  %_____________________________________________________
  % 6. Compute the face-centered pressure
  if (P.debugSteps)
    fprintf('Step 6: compute the face-centered pressure\n');
  end

  for j = G.first_FC:G.last_FC                                   
    L  = j-1; 
    R = j;
    press_FC(j) = (press_CC(L)* rho_CC(R) + press_CC(R) * rho_CC(L)) / ...
                   (rho_CC(L) + rho_CC(R) + d_SMALL_NUM);
  end


  %_____________________________________________________
  % 7. Accumulate sources
  % We neglect tortion forces, gravity effects, and the
  % temperature gradient term.
  if (P.debugSteps)
    fprintf('Step 7: accumulate sources\n');
  end
  
  for j = G.first_CC:G.last_CC
    L = j;
    R = j+1;
    
    del_mom(j) = -delT * (press_FC(R) - press_FC(L));
    del_eng(j) = -delT * (press_FC(R).*xvel_FC(R) - press_FC(L).*xvel_FC(L));
  end


  %_____________________________________________________
  % 8. Compute Lagrangian quantities
  if (P.debugSteps)
    fprintf('Step 8: compute Lagrangian quantities\n');
  end

  mass_L    = rho_CC .* G.delX;                  % m^L = rho * cell_volume;
  mom_L     = mass_L .* xvel_CC + del_mom;       % (mu)^L = (mu) + del_(mu)
  eng_L     = mass_L .* E_CC + del_eng;          % (mE)^L = (mE) + del_(me)

  % Translate fluxed quantities to primitive variables
  rho_CC    = mass_L ./ G.delX;
  xvel_CC   = mom_L ./ (mass_L + d_SMALL_NUM);
  E_CC      = eng_L ./ (mass_L + d_SMALL_NUM);

  % Set boundary conditions for primitive variables
  rho_CC    = setBoundaryConditions(rho_CC, 'rho_CC',  G);
  xvel_CC   = setBoundaryConditions(xvel_CC,'xvel_CC', G);
  E_CC      = setBoundaryConditions(E_CC,'e_CC', G);

  % Translate primitives into fluxed quantities, now with the "B.C." set for them too
  mass_L    = rho_CC .* G.delX;
  mom_L     = mass_L .* xvel_CC;
  eng_L     = mass_L .* E_CC;

  %_____________________________________________________
  % 9. Advect and advance in time
  if (P.debugSteps)
    fprintf('Step 9: advect and advance in time\n');
  end

  %==================================
  % M A S S
  % Uses van Leer limiter
  if (P.debugSteps)
    fprintf ('Advecting density\n');
  end
  [q_advected, gradLim, grad_x, mass_slab, mass_vrtx_1, mass_vrtx_2] = ...
    advectRho(mass_L, ofs, rx, xvel_FC, G);
  
  mass_CC     = mass_L + q_advected;                        % Advection of rho*ofs (the advected volume) = advection correction to the mass
  rho_CC      = mass_CC ./ G.delX;                          % Updated density
  rho_CC      = setBoundaryConditions(rho_CC,'rho_CC', G);  % We need to set B.C. on rho,T,u

  %==================================
  % M O M E N T U M
  % Uses compatible flux limiter
  if (P.debugSteps)
    fprintf ('Advecting momentum\n');
  end
  [q_advected, gradLim, grad_x] = ...
    advectQ(mom_L, mass_L, mass_slab, mass_vrtx_1, mass_vrtx_2, ofs, rx, xvel_FC, G);
    
  xvel_CC     = (mom_L + q_advected) ./ (mass_CC);           % Updated velocity
  xvel_CC     = setBoundaryConditions(xvel_CC,'xvel_CC', G);

  %==================================
  % E N E R G Y
  % Uses compatible flux limiter
  if (P.debugSteps)
    fprintf ('Advecting energy\n');
  end
  [q_advected, gradLim, grad_x] = ...
    advectQ(eng_L, mass_L, mass_slab, mass_vrtx_1, mass_vrtx_2, ofs, rx, xvel_FC, G);
  
  E_CC     = (eng_L + q_advected)./(mass_CC);        % Updated energy
  E_CC     = setBoundaryConditions(E_CC,'e_CC', G);
  
  e_CC     = E_CC -  (0.5 .* xvel_CC .* xvel_CC);
  temp_CC  = e_CC./P.cv;

  %_____________________________________________________
  % 10. End of the timestep

  %================ Compute del_T ================
  % A G G R E S S I V E
  
  if (t >= P.maxTime)   
    fprintf('Reached maximum time\n');
    break;
  end

  delt_CFL        = 1e+30;
  for j = G.first_CC:G.last_CC
    speed_Sound = speedSound_CC(j);
    A           = P.CFL*G.delX/(speed_Sound + abs(xvel_CC(j)));
    delt_CFL    = min(A, delt_CFL);
  end
  
  if (P.debugSteps)
    fprintf('Aggressive delT Based on currant number CFL = %.3e\n',delt_CFL);
  end

  %================ Compare with Uintah ICE and Plot Results ================
  if (P.compareUintah)
    loadUintah;                             % Load Uintah results
  end
  plotResults;                              % Plot results (and also against Uintah, if compareUintah flag is on)

  %================ Various breaks ================

  delT    = delt_CFL;                                             % Compute delT - "agressively" small
 % delT    = 5e-4;
 % fprintf ('-------------WARNING:   delT has been hard coded to compare with Uintah results\n');
  
  t       = t + delT;                                             % Advance time
  if (t >= P.maxTime)
    fprintf('Reached maximum time\n');
    break;
  end

end

tfinal = t - delT;
delX   = G.delX;

if (P.writeData == 1)
  fname = sprintf('matlab_CC_%g.dat', P.nCells);
  fid = fopen(fname, 'w');
  fprintf(fid,'time %15.16E\n',t-delT)
  fprintf(fid,'X_CC \t press_eq \t delP \t press_CC \t xvel_CC \t temp_CC \t rho_CC\n');

  for c=1:length(x_CC)
    fprintf(fid,'%16.15E %16.15E %16.15E %16.15E %16.15E %16.15E %16.15E\n',x_CC(c), press_eq_CC(c), delPDilatate(c), press_CC(c), xvel_CC(c), temp_CC(c), rho_CC(c));
  end
  fclose(fid);
  
  fname = sprintf('matlab_FC_%g.dat', P.nCells);
  fid = fopen(fname, 'w');
  fprintf(fid,'time %15.16E\n',t-delT)
  fprintf(fid,'X_FC \t xvel_FC \t press_FC\n');

  for c=1:length(x_FC)
    fprintf(fid,'%16.15E %16.15E %16.15E \n',x_FC(c), xvel_FC(c), press_FC(c) );
  end
  fclose(fid);
end

if(P.plotResults)
  plotResults;
  figure(1);
  print -depsc iceResult1.eps
  figure(2);
  print -depsc iceResult2.eps
end
% Show a movie of the results
%hFig = figure(2);
%movie(hFig,M,1,10)
