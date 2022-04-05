% Basic code associated with the paper "Constraint-Driven Coordinated
% Control of Multi-Robot Systems"
% 
% https://arxiv.org/abs/1811.02465
% 
% Written by Gennaro Notomista, 2019

clc
clear
close all

% flags
PLOT = true;

% constants
N = 6;
DT = 0.01;
T = 10;

% formation control parameters
L = [3 -1 -1 -1 0 -1 ; ...
    -1 3 -1 0 -1 0 ; ...
    -1 -1 3 -1 0 -1 ; ...
    -1 0 -1 3 -1 0 ; ...
    0 -1 0 -1 3 -1 ; ...
    -1 0 -1 0 -1 3];
l = 0.4;
weights = [0 l sqrt(3)*l 2*l 0 l; ...
    l 0 l 0 2*l 0; ...
    sqrt(3)*l l 0 l 0 2*l; ...
    2*l 0 l 0 l 0; ...
    0 2*l 0 l 0 l; ...
    l 0 2*l 0 l 0];

% multi-robot system
robots = cell(1,N);
for i = 1 : N
    robots{i} = SingleIntegrator('width',0.05,...
        'initialState',-1+2*rand(2,1),...
        'simulationTimeStep',DT);
end
s = Swarm('robots',robots,'L',L);

% optimization parameters
opts = optimoptions(@quadprog,'Display','off');
H = eye(2*N);
f = zeros(1,2*N);
Kformation = 1;
gamma = 1; % increase to get faster task convergence

% init plots
s.plotFigure()

% main loop
for t = 0 : DT : T
    tic
    
    % get robot poses
    q = s.getPoses();
    q = eye(2,3) * q; % extract positions only
    
    % evaluate controller
    A = zeros(1,2*N);
    b = 0;
    for i = 1 : N
        for j = s.getNeighbors(i)
            A(1,2*(i-1)+1:2*(i-1)+2) = A(1,2*(i-1)+1:2*(i-1)+2) + Kformation*2*(norm(q(:,i)-q(:,j))^2-weights(i,j)^2)*(q(:,i)-q(:,j))';
            b = b - gamma*Kformation*(norm(q(:,i)-q(:,j))^2-weights(i,j)^2)^2;
        end
    end
    u = quadprog(H, f, A, b, [],[], [], [], [], opts);
    u = reshape(u, 2, N);
    
    % integration step
    s.moveSingleIntegrators(u)
    
    % update plots
    if PLOT
        s.plotGraph('Color',[0,0.145,0.298],'LineWidth',5)
    end
    s.plotRobots([0.933,0.698,0.067],'EdgeColor','none')
    
    drawnow limitrate
    pause(DT-toc)
end