%% Clear workspace and Command Window Screen 
clear all 
close all 
clc
%% PSO Parameters
swarm_size = 8;                                    % Swarm Size
no_design_variable = 4;                             % # of Design Variables

x_max = 3;                                          % Maximum Position

iter_max = 30;                                     % Maximum no. of Iterations
iter = 0;                                           % Current Iterataion
particle_best = cell(1,swarm_size);                 % Particle Best Position
particle_best_objective = ones(1,swarm_size)*1E50;  % Particle Best Position Objective Value
global_best_objective = 1E50;                       % Global Best Position Objective Value
global_best = zeros(1,no_design_variable);          % Global Best Position
particle_position = cell(1,swarm_size);             % Particle Position
obj_fun_val_particle = zeros(1,swarm_size);         % Objective Function Value of the Particle
particle_velocity = cell(1,swarm_size);             % Particle_Velocity
dummy = zeros(1,length(no_design_variable));        % Dummy List
dt = 1;                                             % Time Step
%% Pressures Vessel Design Problem Definition  

PCONST = 1000000;                                                               % Penalty Function Constant


objective_function = @(x) 0.6224.*x(1).*x(3).*x(4) + 1.7781.*x(2).*(x(3)^2) + 3.1661 .*(x(1)^2).*(x(4)) + 19.84 .*(x(1)^2).*(x(3));      % Objective Function f(x)


g1 = @(x) -x(1) + 0.0193*x(3);                                             
g2 = @(x) -x(2) + 0.00954*x(3) ;                                                 
g3 = @(x) -pi*x(3)^2*x(4)^2 - (4/3)*pi*x(3)^3 + 1296000;
g4 = @(x) x(4) - 240;

penalty_function = @(x) objective_function(x) + PCONST*(max(0,g1(x))^2 + max(0,g1(x))^2 + max(0,g2(x))^2 + max(0,g3(x))^2 + max(0,g4(x))^2);                                             % Penalty Function
%% PSO Initialize Particles' Positions & Velocities
for i = 1:swarm_size
    for j = 1:no_design_variable
        dummy(j) = rand()*x_max;
    end
    particle_position{i} = dummy;
    particle_velocity{i} = particle_position{i}/dt;
end
%% PSO Main Loop
while iter < iter_max
    
    C1 = 1.8;   % personal learning factor
    C2 = 1.8;   % social learning factor
    W = 0.8;    % inertia factor
    
    iter = iter + 1;
    for i = 1:swarm_size
        obj_fun_val_particle(i) = penalty_function(particle_position{i});
        if obj_fun_val_particle(i) < particle_best_objective(i) && obj_fun_val_particle(i) >= 0  % Best Local
            particle_best{i} = particle_position{i};
            particle_best_objective(i) = obj_fun_val_particle(i);
        end
    end
    if min(obj_fun_val_particle) < global_best_objective && min(obj_fun_val_particle) >= 0       % Best Global
        global_best = particle_position{obj_fun_val_particle == min(obj_fun_val_particle)};
        global_best_objective = min(obj_fun_val_particle);
    end
    for i = 1:swarm_size    % Update Particle Positon
        R1 = rand();        % Random Number
        R2 = rand();        % Random Number
        particle_velocity{i} = W*particle_velocity{i}+C1*R1*(particle_best{i}-particle_position{i})/dt+C2*R2*(global_best-particle_position{i})/dt; % Update Particle Velocity
        particle_position{i} = particle_position{i}+particle_velocity{i}*dt; % Update Particle Position
    end
    disp(['BEST PARTICLE VALUE >> ' num2str(global_best_objective)]);
end
disp(['BEST PARTICLE POSITION >> ' num2str(global_best)]);
