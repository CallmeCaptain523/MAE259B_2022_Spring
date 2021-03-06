function [all_mid_y, all_mid_v, all_yq,all_yv] = P2_simu(N, dt, totalTime)

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
    ne = N - 1;                                  % Number of edges
    EI = Y * pi * r0^4 / 4;                      % Stretching Stiffness
    EA = Y * pi * r0^2;                          %  Bending Stiffiness
    
    %Sphere Properties
    R_mid = 0.025;
    R_rest = deltaL/10;
%     m_rest = 4/3*pi*rho_metal*R_rest^3;
%     m_mid = 4/3*pi*rho_metal*R_mid^3;
    
    % Geometry
    nodes = zeros(N, 2);                        
    for c = 1:N
        nodes(c,1) = (c-1) * deltaL;
    %   nodes(c,2) = 0;
    end

    % Mass Matrix
    M = zeros(2*N, 2*N);
    for i = 1:2:2*N
        M(i,i) = 4/3*pi*R_rest^3*rho_metal;
        %M(i+1,i+1) = 4/3*pi()*R_rest^3*densMetal;
        M(i+1,i+1) = M(i,i);
    end
    M(N,N) = 4/3*pi()*R_mid^3*rho_metal;
%   M(N+1,N+1) = 4/3*pi()*R_mid^3*rho_metal;
    M(N+1,N+1) = M(N,N);
    
    % Viscous damping matrix
    C = zeros(2*N,2*N);
    C_rest = 6*pi*visc*R_rest;
    C_mid = 6*pi*visc*R_mid;
    
    for i = 1:2*N
        C(i,i) = C_rest;
    end
    C(N,N) = C_mid;
    C(N+1,N+1) = C_mid;
    
    % Gravity
    W = zeros(2*N,1);
    for i = 2:2:2*N
        W(i,1) = -4/3*pi*R_rest^3*rho*g;
    end
    W(N+1)= -4/3*pi*R_mid^3*rho*g;
    
    % Initial DOF vector
    q0 = zeros(2*N,1);
    for c=1:N
        q0 (2*c - 1) = nodes(c,1);        %Initial Shape
        q0 (2*c) = nodes(c,2);
    end
    
    % Position and Velocity Record
    q = q0;                               % Old Position vector
    u = (q - q0) / dt;                    % Velocity Vector
    
    % Time steps
    Nsteps = round(totalTime / dt);
    all_mid_y = zeros(Nsteps, 1);      % y - position
    all_mid_v = zeros(Nsteps, 1);      % y - velocity
    
    all_mid_y(1) = q(N+1);
    all_mid_v(1) = u(N+1);
    
    tol = EI / RodLength^2 * 1e-3;

    % Time Matching
    for c = 2:Nsteps
        fprintf('Time = %f\n', (c-1)*dt);
        q = q0;     % Guess
        
        % Newton Raphson
        err = 10 *tol;
        while err > tol
            % Inertia
            f = M/dt * ( ((q-q0)/dt) -u );
            J = M /dt^2;
            
            % Elastic Forces/ Protential
            % spring # 1 b/t 1-2
            for k = 1:(N-1)
                xk = q(2*k-1);
                yk = q(2*k);
                xkpl = q(2*k+1);
                ykpl = q(2*k+2);
                % l_k = deltaL;
                dF = gradEs (xk, yk, xkpl, ykpl, deltaL, EA);
                dJ = hessEs (xk, yk, xkpl, ykpl, deltaL, EA);
                f((2*k-1):(2*k+2)) = f(2*k-1:2*k+2)+dF;
                J(2*k-1:2*k+2, 2*k-1:2*k+2) = J(2*k-1:2*k+2, 2*k-1:2*k+2) + dJ;
            end
        
            % bending spring b/t node 1, 2, 3
            for k = 2:(N-1)
                xkml = q(2*k-3);
                ykml = q(2*k-2);
                xk = q(2*k-1);
                yk = q(2*k);
                xkpl = q(2*k+1);
                ykpl = q(2*k+2);
                curvature0 = 0;
                dF = gradEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
                dJ = hessEb (xkml, ykml, xk, yk, xkpl, ykpl, curvature0, deltaL, EI);
                f((2*k-3):(2*k+2)) = f((2*k-3):(2*k+2))+dF;
                J((2*k-3):(2*k+2), (2*k-3):(2*k+2)) = J((2*k-3):(2*k+2), (2*k-3):(2*k+2))+dJ;
            end
        
            % Calculate q using N-R
            % Viscous force
            f = f+ C* (q-q0)/dt;
            J = J +C/dt;
            % Weight
            f = f-W;
    
            % Update
            q = q-J\f;
            err = sum(abs(f));
        end
        % update u and q0
        u = (q -q0)/dt;     %velocity
        q0 = q;
    
        % storing value
        all_mid_y(c) = q(N+1);
        all_mid_v(c) = u(N+1);

        all_yq(:,c) = q;
        all_yv(:,c) = u;
    end

end