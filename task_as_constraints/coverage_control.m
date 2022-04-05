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
PLOT_VORONOI = true;
PLOT_CENTROID = true;
PLOT_DENSITY = true;
DENSITY = ''; % '' or 'uniform'

% constants
N = 5;
DT = 0.01;
T = 10;

% multi-robot system and environment with density field
robots = cell(1,N);
for i = 1 : N
    robots{i} = Unicycle('width',0.1,...
        'length',0.1,...
        'initialState',[-0.5;-0.5;0]+1*rand(3,1),...
        'simulationTimeStep',DT);
end
environment = 2*[cos(linspace(0,2*pi,6)); sin(linspace(0,2*pi,6))];
if strcmp(DENSITY, 'uniform')
    phi = 'uniform';
else
    phi = @(x,y) exp(-((x-1).^2+(y-0.6).^2)/0.3) + 0.5*exp(-((x-0.2).^2+(y+0.2).^2)/0.1);
end
s = Swarm('robots',robots,...
    'environment',environment,...
    'densityFunction',phi);
coverageCost = @(p,g) sum(diag((p-g)'*(p-g)));

% optimization parameters
opts = optimoptions(@quadprog,'Display','off');
kappa = 1e2; % weight of the slack variable
H = blkdiag(eye(2*N), kappa);
f = zeros(1,2*N+1);
gamma = 1e2; % increase to get faster task convergence

% init plots
s.plotFigure()
if PLOT_DENSITY
    s.plotDensity(linspace(0,1.5,24), 'LineWidth', 2)
end
s.plotEnvironment('LineWidth', 5, 'Color', [0 0 0])

% main loop
for t = 0 : DT : T
    tic
    
    % get robot poses
    q = s.getPoses();
    [G, area, VC] = s.coverageControl();
    
    % evaluate controller
    A = zeros(1,2*N+1);
    A(end) = -1;
    b = -gamma * coverageCost(q(1:2,:), G);
    for i = 1 : N
        si = s;
        A(1,2*(i-1)+1:2*(i-1)+2) = 2*(q(1:2,i)-G(:,i))'*(eye(2) - dGidpi(si, q(1:2,:), G, i));
    end
    u = quadprog(H, f, A, b, [],[], [], [], [], opts);
    u = reshape(u(1:2*N), 2, N);
    for i = 1 : N
        uNorm = norm(u(:,i));
        if uNorm >= 0.5
            u(:,i) = u(:,i) / uNorm * 0.5;
        end
    end
    
    % integration step
    s.moveSingleIntegrators(u)
    
    % update plots
    s.plotRobots([0.933,0.698,0.067],'EdgeColor','none')
    if PLOT_VORONOI
        s.plotVoronoiCells(VC,'Color',[0.25 0.25 0.25],'LineWidth',2)
    end
    if PLOT_CENTROID
        s.plotCentroids(G,'.','Color',[0.5 0.5 0.5],'MarkerSize',20)
    end
    
    drawnow limitrate
    pause(DT-toc)
end

function dGidpiVal = dGidpi(s, p, G, i)
epsilon = 1e-4;
dGidpiVal = zeros(2);
px = p; px(:,i) = px(:,i) + epsilon * [1;0];
py = p; py(:,i) = py(:,i) + epsilon * [0;1];
[Gx, ~, ~] = s.coverageControl(px);
[Gy, ~, ~] = s.coverageControl(py);
dGidpiVal(:,1) = (Gx(:,i)-G(:,i)) / epsilon;
dGidpiVal(:,2) = (Gy(:,i)-G(:,i)) / epsilon;
dGidpiVal = dGidpiVal';
end