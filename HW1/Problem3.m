% Problem 3

%% Physical parameters
clear all; clc; clf
% Number of vertices
N = 50;

% Time step size
dt = 1e-3;                             % second
% Total time
totalTime = 1;                         % seconds

% Rod length
L = 1; % meter
d = 0.75; % load applied at

% Load
P = 2000; % Newtons (specified values of 2,000 and 20,000)

% Discrete length
deltaL = L / (N-1);
L_space = linspace(0,L,50);
L_vec = abs(L_space - d);
% loc = min(L_vec);
% ska = find(L_vec == loc);

%Density
rho = 2.7e3;

% Rod radius
R = 0.013;
r = 0.011;

% Young's modulus
Y = 70e9;                              % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8;                               % m/s^2



% Utility quantities
A = (R^2-r^2)*pi;                      % Cross-section area of beam     
I = (R^4-r^4)*pi/4;                    % Moment inertia of cross section 
%ne = N-1;                              % Number of edges
EI = Y*I;                              % Stretching Stiffness
EA = Y*A;                              %  Bending Stiffiness

% Geometry
nodes = zeros(N, 2);                   %Record the initial position
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
end

%Mass matrix
m = (pi*(R^2-r^2)*L*rho)/(N-1);
M = zeros(2*N, 2*N);
for i = 1:2:2*N
    M(i,i) = m;
    M(i+1,i+1) = m;
end
 
%Pressure location
Px = zeros(2*N,1);
Px(round(N*3/4)*2) = -P;

%Initial DOF vector
q0 = zeros(2*N,1);
for c = 1:N
    q0 (2*c - 1) = nodes(c,1);        % x coordinate
    q0 (2*c) = nodes(c,2);            % y coordinate
end

%position and velocity
q = q0;                                 % Old Position vector
u = (q - q0)/dt;                      % Velocity vector

%Number of time steps
Nsteps = round(totalTime / dt);
ymax = zeros(Nsteps, 1)

%Tolerance
tol = EI / L^2 * 1e-3;

%Time marching scheme
for i=2:Nsteps
    
    fprintf('Time = %f\n', (i-1) * dt );
    
    q = q0; % Guess
    
    %Newton Raphson
    err = 10 * tol;
    while err > 1e-5
        %Inertia
        f = M / dt * ((q-q0)/dt - u );
        J = M / dt^2;
        
        %Elastic forces
        %Linear spring 1 between nodes 1 and N-1
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
            J(2*k-1:2*k+2,2*k-1:2*k+2)=J(2*k-1:2*k+2,2*k-1:2*k+2) + dJ;
        end
        
        %Bending spring between nodes 
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
            J(2*k-3:2*k+2,2*k-3:2*k+2) = J(2*k-3:2*k+2,2*k-3:2*k+2) + dJ;
        end
        
        %Load
        f = f - Px;
        
        %Update
        q(3:2*N-1) = q(3:2*N-1) - J((3:2*N-1),(3:2*N-1)) \ f(3:2*N-1);
        
        err = sum(abs(f(3:2*N-1)));
    end

    %Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
    
    %simulation of beam
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis([0 1 -0.3 0])
    drawnow
    grid on
    
    %Store
    ymax(i) = min(q);
end
% Plotting max vertical displacement vs time
figure(2);
timeArray = linspace(0,totalTime,Nsteps);
plot(timeArray, ymax, 'b-');
xlabel('Time [s]');
ylabel('Maximum Vertical Displacement, y_max [m]');
title('Maximum Vertical Displacement vs. Time');
grid on
%% Euler beam theory
c = min(d, L-d);
y_euler = (P*c*(L^2 - c^2)^1.5) / (9*sqrt(3)*EI*L);
fprintf('y_max from Euler beam theory is -%0.4f [m] \n' ,y_euler);
% Difference between simulation and Euler beam theory
max_disp = max(abs(ymax));
diff = abs(max_disp - y_euler);