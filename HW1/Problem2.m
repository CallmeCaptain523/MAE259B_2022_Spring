% MAE259B HW1
% Jiahui Xi
%  UID：005730084
% % Problem 2 ：Simulation of the motion of 21 sphere drop in viscous fluid (Section 4.3)

clear all; clc;clf;

N = 21;
dt = 1e-2;      % second
totalTime = 10;
Nsteps = round(totalTime / dt);

% Simulation
[all_mid_y, all_mid_v, all_yq,all_yv] = P2_simu(N, dt, totalTime);    %Use Implicit method to calculate the position and velocity
figure(1)
for i = 1:1:length(all_yq)
    plot(all_yq(1:2:end,i), all_yq(2:2:end,i), 'ro-');               
    axis equal
    grid on;
    drawnow    
end

%% Q1: Vertical position and velocity of middle node
figure (2)
timeArray =  dt*(1:Nsteps);
plot(timeArray, all_mid_y, 'b-');
xlabel('Time [s]');
ylabel('Vertical Position of Mid-node [m]');
title('Position of Mid-node vs. Time');
grid on

figure (3)
plot(timeArray, all_mid_v, 'b-');
xlabel('Time [s]');
ylabel('Velocity of Mid-Node vs. [m/s]');
title('Velocity of Mid-Node vs. Time');
grid on

%% Q2: final deformed shape of the beam
figure (4)
plot(all_yq(1:2:end,length(all_yq)), all_yq(2:2:end,length(all_yq)), 'ro-');
xlabel('X [m]');
ylabel('Y [m]');
%axis equal
title('Final Deformed Shape');
% hold on 
grid on

%% Q3: terminal velocity vs. the number of nodes
num_Node = 25 ;
% dt = 1e-2;
all_yq_c_n = cell(num_Node,1);
all_yv_c_n = cell(num_Node,1);

for N = 3:2:num_Node*2+1
    [~,~,all_yq_c_n{N-2}, all_yv_c_n{N-2}] = P2_simu(N, dt, totalTime);
end

Nodes_cell = 3:2:num_Node*2+1;
terminal_velo = zeros(num_Node,1);

for i = 1:2:num_Node*2
    [wid,len] = size(all_yq_c_n{i});
    terminal_velo(i) = all_yv_c_n{i}(wid/2+1,len);
end
terminal_velo = nonzeros(terminal_velo);

figure (5)
axis equal;
plot(Nodes_cell,terminal_velo, 'b*-');
xlabel('Number of Nodes');
ylabel('Terminal Velocity [m/s]');
title('Terminal Velocity vs. Number of Nodes')
grid on

%% Q3: terminal velocity vs. Step_size
num_step = 6;
N = 21;
all_yq_c_s = cell(num_step,1);
all_yv_c_s = cell(num_step,1);

for i = 1:num_step
    dt = 1 * 10^(-i);
    [~,~,all_yq_c_s{i}, all_yv_c_s{i}] = P2_simu(N, dt, totalTime);
end

Nodes_cell = 1:num_step;
terminal_velo = zeros(num_step,1);

for i = 1:num_step
    [wid,len] = size(all_yq_c_s{i});
    terminal_velo(i) = all_yv_c_s{i}(wid/2+1,len);
end

figure (6)
plot(Nodes_cell,terminal_velo, 'b*-');
xlabel('Step size (10^(^-^N^))');
ylabel('Terminal Velocity [m/s]');
title('Terminal Velocity vs. StepSize')
grid on