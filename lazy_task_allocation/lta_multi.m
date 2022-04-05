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
N = 6;
M = 2;
% PI_star = [1, 0]';
% PI_star = [0, 1]';
PI_star = [1/2, 1/2]';

ltaMulti = LTAMulti('N', N, ...
    'M', M, ...
    'PI_star', PI_star, ...
    'C', 1e3, ...
    'l', 1e-3, ...
    'K', 1e3, ...
    'delta_max', 20, ...
    's', [1 1; 1 1; 1 1; 1 1; 1 1; 1 1]);

% multi-robot system and environment
w = 3;
h = 2;
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
optimization_data.y = [linspace(0.2,w-0.2,N); 0.2*ones(1,N)];

% plot
s.plotFigure()
plot(optimization_data.y(1,:), optimization_data.y(2,:), '.', ...
    'MarkerSize', 50, 'MarkerEdgeColor', [.5 0 0], 'MarkerFaceColor', [.5 0 0])

% main loop
for t = 0 : DT : T
    
    counter = counter + 1;
    
    q = s.getPoses();
    [G, A, VC] = s.coverageControl();
    
    optimization_data.p = q(1:2,:);
    optimization_data.A = A;
    optimization_data.G = G;
    for i = 1 : N
        optimization_data.coverage_cost(i) = s.evaluateCoverageCost(q,VC,i);
    end
    
    % [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaMulti.solve_miqp(optimization_data);
    [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaMulti.solve_qp_relax(optimization_data);
    % [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = ltaMulti.solve_qp_sdp_relax(optimization_data);
    u = reshape(u,2,N);
    
    s.moveSingleIntegrators(u);
    
    if mod(counter, round(speedUpFactor)) == 0
        s.plotRobots([0.933,0.698,0.067],'EdgeColor','none')
        s.plotVoronoiCells(VC,'Color',[0.25 0.25 0.25],'LineWidth',2)
        s.plotCentroids(G,'.','Color',[0.5 0.5 0.5],'MarkerSize',20)
        drawnow limitrate
    end
    
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