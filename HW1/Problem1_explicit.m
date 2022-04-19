%% MAE259B HW1
%% Jiahui Xi
%% UID：005730084
% 
% 
% Problem 1：Simulation of the motion of a sphere falling inside viscous fluid 
% (see Section 4.2)
% 
% 1: (i) Explicit Method

clear all; clc; 
format short;

%% Physical parameters
N = 3;                                        % Number of nodes
dt = 1e-5;                                    % Time step size second
R1 = 0.005; R2 = 0.025; R3 = 0.005;           % Radius of spheres

% Density
rho_metal = 7000;                            
rho_f = 1000;
rho = rho_metal - rho_f;

%Rod Properties
RodLength = 0.1;                             % Rod length [m]
r0 = 0.001;                                  % Rod radius
deltaL = RodLength / (N-1);                   % Discrete length
Y = 1e9;                                     % Using Y instead of E to avoid ambiguity
g = 9.8;                                     % Gravity m/s^2
visc = 1000;                                 % Viscosity Pa-s
totalTime = 10;                              % Total time seconds
ne = N - 1;                                  % Number of edges
EI = Y * pi * r0^4 / 4;                      % Stretching Stiffness
EA = Y * pi * r0^2;                          %  Bending Stiffiness

% Geometry
nodes = zeros(N, 2);                        
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
%     nodes(c,2) = 0;
end

% Mass matrix
% M(1,1) = 4/3*pi*R1^3*rho_metal;
M = zeros(2*N,2*N);                         
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

 % Viscous damping matrix
C = zeros(6,6);                           
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity Vector
W = zeros(2*N,1);                          
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1);                       % x coordinate
    q0 ( 2*c ) = nodes(c,2);                           % y coordinate
end

% New position and velocity
q = q0                                              % DOF vector
u = (q - q0) / dt                                 % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros(6,Nsteps);                     % y-position of R2
all_mid_v = zeros(6,Nsteps);                     % y-velocity of R2

all_mid_y(:,1) = q;
all_mid_v(:,1) = u;

% Tolerance
%tol = EI / RodLength^2 * 1e-3;

%Time marching scheme
for c=2:Nsteps
    
    fprintf('Time = %f\n', (c-1) * dt );
    f = zeros(6,1);
 
   
    % Compute q(tk+1) directly by explicit method
    % Elastic forces
    % Linear spring 1 between nodes 1 and 2
    xk = q(1);
    yk = q(2);
    xkp1 = q(3);
    ykp1 = q(4);
    dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    %dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
    f(1:4) = f(1:4) - dF;
    %J(1:4,1:4) = J(1:4,1:4) + dJ;
    
    % Linear spring 2 between nodes 2 and 3
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    %dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
    f(3:6) = f(3:6) - dF;
    %J(3:6,3:6) = J(3:6,3:6) + dJ;
    
    % Bending spring between nodes 1, 2, and 3
    xkm1 = q(1);
    ykm1 = q(2);
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    curvature0 = 0;
    dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1,curvature0, deltaL, EI);
    %dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1,curvature0, deltaL, EI);
    f(1:6) = f(1:6) - dF;
    %J(1:6,1:6) = J(1:6,1:6) + dJ;
    
    % Viscous force
    visc = C*u; 
    f(1:6) = f(1:6)-visc;

    % Weight
    f = f+W;

    % Intermidia step 1
    step1 = dt*inv(M) *f +u;

    % Intermidia step 2
    q = step1*dt + q;

    % update
    u = (q - q0)/dt;     % velocity
    q0 = q;             % Position

    % storing value
    all_mid_y(:,c) = q;   %Storing the value for old positon & Velocity
    all_mid_v(:,c) = u;
end

%%
figure(1);

for i = 1:10000:length(all_mid_y)
    title('Simulation of position vs. Time')
    plot(all_mid_y(1:2:end,i), all_mid_y(2:2:end,i), 'ro-');
    axis equal
    drawnow
end

figure(2);
timeArray = dt*(1:Nsteps);
plot(timeArray, all_mid_v(N+1,:), 'k-');
title('Velocity of sphere vs. Time')
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');