% Basic code associated with the paper "An Optimal Task Allocation Strategy
% for Heterogeneous Multi-Robot Systems"
% 
% https://arxiv.org/abs/1903.08641
% 
% Written by Gennaro Notomista and Siddharth Mayya, 2019

clc
clear
close all

% constants
DT = 0.01;
T = 60;
N = 3;
M = 3;
% PI_star = [1, 0, 0]';
% PI_star = [0, 1, 0]';
% PI_star = [0, 0, 1]';
% PI_star = [1/2, 1/2, 0]';
% PI_star = [1/2, 0, 1/2]';
% PI_star = [0, 1/2, 1/2]';
PI_star = [1/3, 1/3, 1/3]';

ltaSingle = LTASingle('N', N, ...
    'M', M, ...
    'PI_star', PI_star, ...
    'C', 1e3, ...
    'l', 1e-3, ...
    'K', 1e6, ...
    'delta_max', 10, ...
    's', [1 0 0; 0 1 0; 0 0 1]);

% multi-robot system and environment
w = 1;
h = 1;
environment = [0 0;w 0;w h;0 h]';
robots = cell(N,1);
for i = 1 : N
    robots{i} = SingleIntegrator('initialState', [w;h].*rand(2,1), ...
        'width', .04, ...
        'simulationTimeStep', DT);
end
s = Swarm('robots',robots,...
    'environment',environment);

% loop variables initialition
counter = 1;
speedUpFactor = 2;
printSpeedUpFactor = 20;

% other variables initialition
p1 = [0.35; 0.5];
p2 = [0.75; 0.75];
p3 = [0.85; 0.35];
optimization_data.p1 = p1;
optimization_data.p2 = p2;
optimization_data.p3 = p3;

% plot
s.plotFigure()
plot([p1(1) p2(1) p3(1)], [p1(2) p2(2) p3(2)], '.', 'MarkerSize', 50)
text(p1(1)-.2,p1(2),'$p_1$','interpreter','latex','FontSize',40)
text(p2(1)-.2,p2(2),'$p_2$','interpreter','latex','FontSize',40)
text(p3(1)-.2,p3(2),'$p_3$','interpreter','latex','FontSize',40)

% main loop
for t = 0 : DT : T
    
    counter = counter + 1;
    
    q = s.getPoses();
    
    optimization_data.p = q(1:2,:);
    
    % [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaSingle.solve_miqp(optimization_data);
    [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaSingle.solve_qp_relax(optimization_data);
    % [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaSingle.solve_qp_sdp_relax(optimization_data);
    u = reshape(u,2,N);
    
    s.moveSingleIntegrators(u);
    
    s.plotRobots([0.933,0.698,0.067],'EdgeColor','none')
    drawnow limitrate
    
    if mod(counter, round(printSpeedUpFactor)) == 0
        clc
        disp('================================================')
        disp('QP:')
        disp(['    (*) time to solve: ', num2str(time_to_solve_qp)])
        switch exit_flag
            case 1
                disp('    (*) converged')
            case 0
                disp('    (*) max iter reached')
            case -2
                disp('    (*) INFEASIBLE')
        end
        disp(['    (*) lambda: ',num2str((lambda>1e-6)')])
        disp('')
        disp('alpha:')
        disp(alpha')
        disp('u:')
        disp(u)
        disp('delta:')
        disp(delta')
        disp('================================================')
    end
    
end
