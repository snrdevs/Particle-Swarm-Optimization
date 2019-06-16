%% Clear workspace and Command Window Screen 
clear all 
close all 
clc
%% PSO Parameters
swarm_size = 8;                                    % Swarm Size
no_design_variable = 7;                             % # of Design Variables

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
%% Speed Reducer Design Problem Definition  

PCONST = 1000000;                                                               % Penalty Function Constant

objective_function = @(x) 0.7854*x(1)*x(2)^2 * (3.3333*x(3)^2  + 14.9335*x(3) -  43.0934)...
		-1.508*x(1)*(x(6)^2 + x(7)^2) + 7.4777*(x(6)^3 + x(7)^3)...
		+0.7854*(x(4)*x(6)^2 + x(5)*x(7)^2);                                    % Objective Function f(x)

g1  = @(x) 	27 / (x(1)*x(2)^2*x(3)); 
g2  = @(x) 	(397.5 / (x(1) * x(2)^2 * x(3)^2) )- 1; 
g3 	= @(x) 	((1.93 * x(4)^4) / (x(2) * x(3) * x(6)^4) ) -1 ; 
g4 	= @(x) 	((1.93 * x(5)^3) /(x(2) * x(3) * x(7)^4)) - 1; 

g5 	= @(x)	((1/(110*x(6)^3)) *(sqrt((745*x(4)  / (x(2)*x(3))  )^2 + 16.9e6)) ) - 1; 
g6 	= @(x)	((1/(85*x(7)^3)) *(sqrt((745*x(5)  / (x(2)*x(3))  )^2 + 157.5e6)) ) - 1; 

g7 	= @(x) 	(x(2)*x(3) / 40 )- 1;  
g8 	= @(x)	(5*x(2) / x(1)) -1; 
g9 	= @(x)	(x(1) / (12 * x(2))) - 1;
g10 = @(x)	((1.5*x(6) + 1.9) / (x(4))) -1 ; 
g11 = @(x) 	((1.1*x(7) + 1.9) / (x(5))) -1; 



penalty_function = @(x) objective_function(x) + PCONST*...
   (max(0,g1(x))^2+max(0,g2(x))^2+...
    max(0,g3(x))^2+max(0,g4(x))^2+max(0,g5(x))^2+...
    max(0,g6(x))^2+max(0,g7(x))^2+max(0,g8(x))^2+...
    max(0,g9(x))^2+max(0,g10(x))^2+max(0,g11(x))^2); % Penalty Function
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
