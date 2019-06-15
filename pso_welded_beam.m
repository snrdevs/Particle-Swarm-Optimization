%% Clear workspace and Command Window Screen 
clear all 
close all 
clc
%% PSO Parameters
swarm_size = 60;                                    % Swarm Size
no_design_variable = 4;                             % # of Design Variables
c1_final = 0.5;                                     % Final Self Confidence Parameter
c1_initial = 2.5;                                   % Initial Self Confidence Parameter
c2_final = 2.5;                                     % Final Swarm Confidence Parameter 
c2_initial = 0.5;                                   % Initial Swarm Confidence Parameter
weight_final = 0.4;                                 % Final Inertia Weight
weight_initial = 0.9;                               % Initial Inertia Weight
x_max = 3;                                          % Maximum Position
iter_max = 300;                                     % Maximum no. of Iterations
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
%% Welded Beam Design Problem Definition  

P = 6000;                                                                       % Applied Tip Load
E = 30e6;                                                                       % Young Modulus
G = 12e6;                                                                       % Shear Modulus
L = 14;                                                                         % Length of the Beam

PCONST = 1000000;                                                               % Penalty Function Constant
TAUMAX = 13600;                                                                 % Max Allowed Shear Stress
SIGMAX = 30000;                                                                 % Max Allowed Bending Stress
DELTMAX = 0.25;                                                                 % Max Allowed Tip Perfection

M = @(x) P*(L+x(2)/2);                                                          % Bending Moment at Welding Point
R = @(x) sqrt((x(2)^2)/4+((x(1)+x(3))/2)^2);                                    % Constant
J = @(x) 2*(sqrt(2)*x(1)*x(2)*((x(2)^2)/12+((x(1)+x(3))/2)^2));                 % Polar Moment of Inertia

objective_function = @(x) 1.10471*x(1)^2*x(2)+0.04811*x(3)*x(4)*(14+x(2));      % Objective Function f(x)
sigma = @(x) (6*P*L)/(x(4)*x(3)^2);                                             % Bending Stress
delta = @(x) (4*P*L^3)/(E*x(4)*x(3)^3);                                         % Tip Deflection

Pc = @(x) 4.013*E*sqrt((x(3)^2*x(4)^6)/36)*(1-x(3)*sqrt(E/(4*G))/(2*L))/(L^2);  % Buckling Load
tau_p = @(x) P/(sqrt(2)*x(1)*x(2));                                             % Tau_prime
tau_pp = @(x) (M(x)*R(x))/J(x);                                                 % Tau_double_prime
tau = @(x) sqrt(tau_p(x)^2+2*tau_p(x)*tau_pp(x)*x(2)/(2*R(x))+tau_pp(x)^2);     % Tau (Shear Stress) 

g1 = @(x) tau(x)-TAUMAX;                                                    
g2 = @(x) sigma(x)-SIGMAX;                                                  
g3 = @(x) x(1)-x(4); 
g4 = @(x) 0.10471*x(1)^2+0.04811*x(3)*x(4)*(14+x(2))-5; 
g5 = @(x) 0.125-x(1); 
g6 = @(x) delta(x)-DELTMAX; 
g7 = @(x) P-Pc(x); 

penalty_function = @(x) objective_function(x) + PCONST*(max(0,g1(x))^2+max(0,g1(x))^2+max(0,g2(x))^2+...
    max(0,g3(x))^2+max(0,g4(x))^2+max(0,g5(x))^2+...
    max(0,g6(x))^2+max(0,g7(x))^2);                                             % Penalty Function
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
    
    C1 = (c1_final-c1_initial)*(iter/iter_max)+c1_initial;                      % Change self confidence
    C2 = (c2_final-c2_initial)*(iter/iter_max)+c2_initial;                      % Change swarm confidence 
    W = (weight_final-weight_initial)*(iter/iter_max)+weight_initial;           % Change inertia weight
    iter = iter + 1;
    for i = 1:swarm_size
        obj_fun_val_particle(i) = penalty_function(particle_position{i});
        if obj_fun_val_particle(i) < particle_best_objective(i) & obj_fun_val_particle(i) >= 0  % Best Local
            particle_best{i} = particle_position{i};
            particle_best_objective(i) = obj_fun_val_particle(i);
        end
    end
    if min(obj_fun_val_particle) < global_best_objective & min(obj_fun_val_particle) >= 0       % Best Global
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
