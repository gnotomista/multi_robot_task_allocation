% Basic code associated with the paper "A Resilient and Energy-Aware Task
% Allocation Framework for Heterogeneous Multi-Robot Systems"
% 
% https://arxiv.org/abs/2105.05586
% 
% Written by Gennaro Notomista, Siddharth Mayya, and Yousef Emam, 2020

clc
clear
close all

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

% disturbances parameters
global robot_exo_dist task_exo_dist x_mud y_mud
load('data/mud')
t_endogenous = 15;
robot_exo_dist = 4;
task_exo_dist = 2;

% initialize simulation
DT = 0.1;
x = [1.65*ones(1,n_r);
    linspace(-1,1,n_r);
    zeros(1,n_r)];
rta = RTA(scenario_params, opt_params);

% initialize figure
figure('units','pixels','Position',[0 0 1920 1080]), hold on
axis equal, axis([-1.8 1.8 -1.2 1.2]), set(gca, 'Visible', 'off')
rectangle('Position', [-1.8 -1.2 3.6 2.4], 'EdgeColor', [.5 .5 .5], 'LineWidth', 5)
scatter(poi(1), poi(2), 750, [.5 0 0], 'p', 'LineWidth', 3)
patch(x_mud, y_mud, [0.5 0.3 0.05], 'EdgeColor', 'none', 'FaceAlpha', 0.5) % low-friction zone
h_perimeter = patch(x_perimeter_ring(0), y_perimeter_ring(0), [0 .75 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
% h_G = scatter(zeros(1,n_r), zeros(1,n_r), 700, [0 .5 0], 'o', 'LineWidth', 3);
scatter(p_goal(1), p_goal(2), 750, [.5 0 0], 'x', 'LineWidth', 3)
h_P = scatter(nan(1,n_r), nan(1,n_r), 750, [.5 0 0], 'o', 'LineWidth', 3);
h_P_traj = line(nan, nan, 'LineStyle', '--', 'Color', [.5 0 0], 'LineWidth', 3);
h_r = scatter(x(1,:), x(2,:), 5000, [0.2 0.4 0.6], '.'); % robots
drawnow

% loop and plot variables
max_iter = 800;
x_traj = zeros(n_x, n_r, max_iter+1); x_traj(:,:,1) = x;
P_traj = nan(2,1);
% u_traj = zeros(n_u*n_r, max_iter);
% V_traj = zeros(1, max_iter);
h_quad = [];
h_fov = [];

for iter = 1 : max_iter
    t = iter * DT;
    
    p_transport_t = p_transport(t);
    P_traj = [P_traj, p_transport_t];
    s.setPoses([x(1:2,:); zeros(1,n_r)]);
    [G,~,VC] = s.coverageControl();
    
    [alpha, u, delta, time_to_synthesize_controller, opt_sol_info] = rta.solve_miqp(x, t);
    
    alpha = reshape(alpha,n_t,n_r);
    u = reshape(u,n_u,n_r);
    delta = reshape(delta,n_t,n_r);
    task_assignment = get_task_assignment(alpha);
    
    for i = 1 : n_r
        % endogenous disturbance
        if t >= t_endogenous && updated_scenario_params.A(3,1) ~= 0
            updated_scenario_params.A(3,1) = 0;
            rta.set_scenario_params(updated_scenario_params);
        end
        
        % exogenous disturbance
        x_sim_i = x(:,i) + (sys_dyn(x(:,i),u(:,i))) * DT;
        if exogenous_disturbance(x,alpha,i)
            x(1:2,i) = x(1:2,i); % integration step for the position
            x(3,i) = x(3,i) + (sys_dyn(x(:,i),u(:,i)))'*[0;0;1] * DT; % integration step for the camera
        else
            x(:,i) = x(:,i) + (sys_dyn(x(:,i),u(:,i))) * DT; % integration step
        end
        S = rta.get_specializations();
        for j = 1 : n_t
            Dh_ij = scenario_params.tasks{j}.function(x(:,i),t,i)-scenario_params.tasks{j}.function(x_sim_i,t,i);
            S(j,i) = max(0, S(j,i) + 10*alpha(j,i)*Dh_ij);
        end
        rta.set_specializations(S);
    end
    
    x_traj(:,:,iter+1) = x;
    % u_traj(:,iter) = reshape(u,n_u*n_r,1);
    h = zeros(2,1);
    for i = 1 : n_r
        if task_assignment(i) ~= 0
            h(task_assignment(i)) = h(task_assignment(i)) + scenario_params.tasks{task_assignment(i)}.function(x(:,i),t,i);
        end
    end
    % V_traj(:,iter) = norm(opt_params.gamma(h))^2;
    
    h_r.XData = x(1,:);
    h_r.YData = x(2,:);
    h_P.XData = p_transport_t(1);
    h_P.YData = p_transport_t(2);
    h_P_traj.XData = P_traj(1,:);
    h_P_traj.YData = P_traj(2,:);
    h_perimeter.XData = x_perimeter_ring(t);
    h_perimeter.YData = y_perimeter_ring(t);
    % h_G.XData = G(1,task_assignment==2);
    % h_G.YData = G(2,task_assignment==2);
    % s.plotVoronoiCells(VC,'Color',[0.25 0.25 0.25],'LineWidth',2)
    h_quad = plot_quad(h_quad, x_traj(:,:,1:iter+1));
    h_fov = plot_fov(h_fov, x, task_assignment);
    drawnow
    
    print_info(t, time_to_synthesize_controller, task_assignment, u, delta, S, t_endogenous)
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

% exogenous disturbance condition
function cond = exogenous_disturbance(x,alpha,i)
global robot_exo_dist task_exo_dist x_mud y_mud
cond = (i == robot_exo_dist) && ... % robot robot_exo_dist and ...
    alpha(task_exo_dist,i) > 0 && ... % ... assigned to task task_exo_dist and ...
    inpolygon(x(1,i),x(2,i),x_mud,y_mud); % ... in the mud
end

% utils
function task_assignment = get_task_assignment(alpha)
[~,task_assignment] = max(alpha,[],1);
for i = 1 : size(alpha,2)
    if norm(alpha(:,i)) < 1e-3
        task_assignment(i) = 0;
    end
end
end

% prints
function print_info(t, time_to_synthesize_controller, task_assignment, u, delta, S, t_endogenous)
clc
disp('time')
disp(t)
disp('time to synthesize controller')
disp(time_to_synthesize_controller)
disp('task assignment')
disp(task_assignment)
disp('input')
disp(u)
disp('delta')
disp(delta)
disp('specialization')
disp(S)
if t >= t_endogenous
    disp('Robot 1 lost feature (camera broke)')
end
end

% plots
function h_quad = plot_quad(h_quad, x_traj)
x = x_traj(1,5,end);
y = x_traj(2,5,end);
dir = x_traj(:,5,end)-x_traj(:,5,end-1);
th = atan2(dir(2),dir(1))+pi/4;
l = 0.18;
r = 0.04;
q1 = [x;y] + l/2*[cos(th); sin(th)];
q2 = [x;y] + l/2*[cos(th+pi/2); sin(th+pi/2)];
q3 = [x;y] + l/2*[cos(th+pi); sin(th+pi)];
q4 = [x;y] + l/2*[cos(th+3/2*pi); sin(th+3/2*pi)];
prop = r*[cos(linspace(0,2*pi,36)); sin(linspace(0,2*pi,36))];
if isempty(h_quad)
    h_quad = cell(1,6);
    h_quad{1} = line([q1(1) q3(1)], [q1(2) q3(2)], 'Color', 'k', 'LineWidth', 2);
    h_quad{2} = line([q2(1) q4(1)], [q2(2) q4(2)], 'Color', 'k', 'LineWidth', 2);
    h_quad{3} = patch(q1(1)+prop(1,:), q1(2)+prop(2,:), [.5 .5 .5], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.5);
    h_quad{4} = patch(q2(1)+prop(1,:), q2(2)+prop(2,:), [.5 .5 .5], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.5);
    h_quad{5} = patch(q3(1)+prop(1,:), q3(2)+prop(2,:), [.5 .5 .5], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.5);
    h_quad{6} = patch(q4(1)+prop(1,:), q4(2)+prop(2,:), [.5 .5 .5], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.5);
else
    h_quad{1}.XData = [q1(1) q3(1)];
    h_quad{1}.YData = [q1(2) q3(2)];
    h_quad{2}.XData = [q2(1) q4(1)];
    h_quad{2}.YData = [q2(2) q4(2)];
    h_quad{3}.XData = q1(1)+prop(1,:);
    h_quad{3}.YData = q1(2)+prop(2,:);
    h_quad{4}.XData = q2(1)+prop(1,:);
    h_quad{4}.YData = q2(2)+prop(2,:);
    h_quad{5}.XData = q3(1)+prop(1,:);
    h_quad{5}.YData = q3(2)+prop(2,:);
    h_quad{6}.XData = q4(1)+prop(1,:);
    h_quad{6}.YData = q4(2)+prop(2,:);
end
end


function h_fov = plot_fov(h_fov, x, task_assignment)
idx_task_2 = find(task_assignment==2);
xy_robot = x(1:2,idx_task_2);
th_cam = x(3,idx_task_2);
L = 100;
TH = pi/72;
if isempty(h_fov)
    h_fov = cell(3);
    for i = 1 : numel(idx_task_2)
        h_fov{i} = patch(xy_robot(1,i)+[0,L*cos(th_cam(i)-TH),L*cos(th_cam(i)+TH)], ...
            xy_robot(2,i)+[0,L*sin(th_cam(i)-TH),L*sin(th_cam(i)+TH)], ...
            [1 1 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
else
    for i = 1 : numel(idx_task_2)
        h_fov{i}.XData = xy_robot(1,i)+[0,L*cos(th_cam(i)-TH),L*cos(th_cam(i)+TH)];
        h_fov{i}.YData = xy_robot(2,i)+[0,L*sin(th_cam(i)-TH),L*sin(th_cam(i)+TH)];
    end
end
end
