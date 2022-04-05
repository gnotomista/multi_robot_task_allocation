% Basic code associated with the paper "A Resilient and Energy-Aware Task
% Allocation Framework for Heterogeneous Multi-Robot Systems"
% 
% https://arxiv.org/abs/2105.05586
% 
% Written by Gennaro Notomista, Siddharth Mayya, and Yousef Emam, 2020

clc
clear
close all

mqtt_interface = MqttInterface('central_unit', 'localhost', 1883, 10);
mqtt_interface.subscribe('robots_state');

% numbers and dimensions
n_r = 5;
n_t = 2; % 1: transport and 2: perimeter defense
n_c = 2; % 1: locomotion and 2: monitoring
n_f = 3; % 1: wheels, 2: propellers, 3: camera
n_x = 3;
n_u = 3;

% robots' features and capabilities
A = zeros(n_f, n_r);
A(1,1) = 1; A(2,1) = 0; A(3,1) = 1; % Robot 1's features
A(1,2) = 1; A(2,2) = 0; A(3,2) = 1; % Robot 2's features
A(1,3) = 1; A(2,3) = 0; A(3,3) = 1; % Robot 3's features
A(1,4) = 1; A(2,4) = 0; A(3,4) = 1; % Robot 4's features
A(1,5) = 0; A(2,5) = 1; A(3,5) = 1; % Robot 5's features
T = zeros(n_t, n_c);
T(1,1) = 1; T(1,2) = 0; % Task 1 (tranport)'s capabilities
T(2,1) = 3; T(2,2) = 3; % Task 2 (perimeter defense)'s capabilities
Hs = cell(1, n_c); % make sure to normalize each row of each cell
Hs{1} = zeros(2, n_f); Hs{1}(1,1) = 1; Hs{1}(2,2) = 1;  % Capability 1 (locomotion)'s features
Hs{2} = zeros(1, n_f); Hs{2}(1,3) = 1;                  % Capability 2 (monitoring)'s features
scenario_params.A = A;
scenario_params.T = T;
scenario_params.Hs = Hs;
scenario_params.ws = [];
updated_scenario_params.A = A;

% robot model
f = @(x) 0*x;
g = @(x) eye(n_x);
sys_dyn = @(x,u) f(x)+g(x)*u;
scenario_params.robot_dyn.f = f;
scenario_params.robot_dyn.g = g;
scenario_params.robot_dyn.n_x = n_x;
scenario_params.robot_dyn.n_u = n_u;

% tasks
scenario_params.tasks = cell(1,n_t);
global p_start p_goal t_start delta_t p_transport_t poi s G
% transport
p_start = [1;-0.6];
p_goal = [-1;0.6];
t_start = 2;
delta_t = 60;
scenario_params.tasks{1}.function = @transport_function;
scenario_params.tasks{1}.gradient = @transport_gradient;
scenario_params.tasks{1}.time_derivative = @transport_time_derivative;
% perimeter defense
poi = [0;1];
robots = cell(1,n_r);
for i = 1 : n_r
    robots{i} = SingleIntegrator();
end
environment = [1.8 1.2; -1.8 1.2; -1.8 -1.2; 1.8 -1.2]';
x_perimeter_ring = @(t) p_transport(t)'*[1;0] + [0.3*cos(linspace(0,2*pi,36)) 0.5*cos(linspace(2*pi,0,36))];
y_perimeter_ring = @(t) p_transport(t)'*[0;1] + [0.3*sin(linspace(0,2*pi,36)) 0.5*sin(linspace(2*pi,0,36))];
phi = @phi_perimeter;
s = Swarm('robots', robots, ...
    'environment', environment, ...
    'densityFunction', phi);
scenario_params.tasks{2}.function = @coverage_control_task_function;
scenario_params.tasks{2}.gradient = @coverage_control_task_gradient;
scenario_params.tasks{2}.time_derivative = @coverage_control_task_time_derivative;

% optimization parameters
opt_params.l = 1e-6; % relative weight delta/u in the cost
opt_params.kappa = 1e6; % scale between tasks with different priorities
opt_params.delta_max = 1e3; % delta_max
opt_params.n_r_bounds = [1 1; 3 3]; % row i is min and max number of robots for task i (n_r_bounds = n_t x 2)
opt_params.gamma = @(x) 5*x; % class K function for task execution

rta = RTA(scenario_params, opt_params);

iter = 0;
avg_time_to_solve_MIQP = 0;

while true
    iter = iter + 1;
    msg_recv = mqtt_interface.receive_json('robots_state');
    if isempty(msg_recv)
        continue
    end
    if isfield(msg_recv,'scenario_params')
        rta.set_scenario_params(msg_recv.scenario_params);
    end
    if isfield(msg_recv,'specialization')
        rta.set_specializations(msg_recv.specialization);
    end
    if isfield(msg_recv,'x') && isfield(msg_recv,'t')
        p_transport_t = p_transport(msg_recv.t);
        s.setPoses([msg_recv.x(1:2,:); zeros(1,n_r)]);
        [G,~,VC] = s.coverageControl();
        tic
        [alpha, ~, ~, ~, ~] = rta.solve_miqp(msg_recv.x, msg_recv.t);
        t_miqp = toc;
        avg_time_to_solve_MIQP = (avg_time_to_solve_MIQP*(iter-1)+t_miqp)/iter;
        disp(['time to solve MIQP: ', num2str(t_miqp)])
        disp(['average time to solve MIQP: ', num2str(avg_time_to_solve_MIQP)])
        msg_send = [];
        msg_send.alpha = alpha;
        mqtt_interface.send_json('allocation', msg_send);
    end
end

% task 1 (transport)
function x = clamp(x,m,M)
x = min(max(x,m),M);
end

function p = p_transport(t)
global p_start p_goal t_start delta_t
if t<t_start
    p = p_start;
elseif t>t_start+delta_t
    p = p_goal;
else
    p = [clamp(1-(t-t_start)/delta_t,0,1)*p_start(1) + clamp((t-t_start)/delta_t,0,1)*p_goal(1);
        clamp(1-(t-t_start)^8/delta_t^8,0,1)*p_start(2) + clamp((t-t_start)^8/delta_t^8,0,1)*p_goal(2)];
end
end

function dp_dt = p_transport_time_derivative(t)
global p_start p_goal t_start delta_t
if t<t_start || t>t_start+delta_t
    dp_dt = [0;0];
else
    dp_dt = [(p_goal(1)-p_start(1))/delta_t;
        8*(t-t_start)^7/delta_t^8 * (p_goal(2)-p_start(2))];
end
end

function hi = transport_function(x_i, t, i)
hi = -norm(x_i(1:2)-p_transport(t))^2;
end

function dhi_dxi = transport_gradient(x_i, t, i)
dhi_dxi = [-2*(x_i(1:2)-p_transport(t))', 0];
end

function dhi_dt = transport_time_derivative(x_i, t, i)
dhi_dt = -2*(x_i(1:2)-p_transport(t))'*p_transport_time_derivative(t);
end

% task 2 (perimeter defense)
function phi = phi_perimeter(x, y)
global p_transport_t
phi = exp(-100*((x-p_transport_t'*[1;0])^2+(y-p_transport_t'*[0;1])^2-0.4^2)^2);
end

function ci = coverage_control_task_function(xi, t, i)
global poi G
ci = -norm(xi(1:2) - G(:,i))^2 - (xi(3)-atan2(poi(2)-xi(2),poi(1)-xi(1)))^2;
end

function dci_dxi = coverage_control_task_gradient(xi, t, i)
global poi G
dci_dxi = [-2*(xi(1:2) - G(:,i))', -2*(xi(3)-atan2(poi(2)-xi(2),poi(1)-xi(1)))];
end

function dci_dt = coverage_control_task_time_derivative(xi, t, i)
dci_dt = 0;
end


