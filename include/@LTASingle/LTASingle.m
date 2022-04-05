classdef LTASingle < handle
    properties (Access = public)
        N
        M
        d
        C
        PI_star
        l
        K
        delta_max
        JM
        SM
        optimization_data
    end
    methods (Access = public)
        function obj = LTASingle(varargin)
            ip = inputParser;
            ip.CaseSensitive = true;
            addParameter(ip, 'N', 3)
            addParameter(ip, 'M', 2)
            addParameter(ip, 'd', 2)
            addParameter(ip, 'C', 1)
            addParameter(ip, 'PI_star', [0.1, 0.2]')
            addParameter(ip, 'l', 1)
            addParameter(ip, 'K', 1)
            addParameter(ip, 'delta_max', 1)
            addParameter(ip, 's', [])
            parse(ip,varargin{:})
            
            obj.N = ip.Results.N;                  % # robots
            obj.M = ip.Results.M;                  % # tasks
            obj.d = ip.Results.d;                  % # dims task space (2: planar robots, 3: aerial robots)
            obj.C = ip.Results.C;                  % weight of alpha w.r.t. u
            obj.PI_star = ip.Results.PI_star;      % desired average tasks violation
            obj.l = ip.Results.l;                  % weight of delta w.r.t. u
            obj.K = ip.Results.K;                  % ratio between tasks at different priorities
            obj.delta_max = ip.Results.delta_max;  % upper bound for delta
            s = ip.Results.s;                      % robots' preferences (NxM matrix, s_{ij} = weight that robot i puts on task j)
            
            obj.JM = repmat(eye(obj.M),1,obj.N);
            obj.SM = [];
            for i = 1 : obj.N
                if isempty(s)
                    obj.SM = blkdiag(obj.SM, eye(obj.M));
                else
                    obj.SM = blkdiag(obj.SM, diag(s(i,:)));
                    Si = diag(s(i,:));
                    % Si(:,~any(Si,1)) = [];
                    % obj.JM(:,(i-1)*obj.M+1:i*obj.M) = Si/(Si'*Si)*Si';
                    obj.JM(:,(i-1)*obj.M+1:i*obj.M) = Si*pinv(Si);
                end
            end
        end
        function [alpha, u, delta, exit_flag, lambda, time_to_solve_qp, cost] = solve_miqp(obj, optimization_data)
            obj.optimization_data = optimization_data;
            [A,b,Aeq,beq,lb,ub] = obj.get_constraints('centralized');
            tic
            cvx_begin quiet
            binary variable alpha_var(obj.N*obj.M)
            variable u_var(obj.N*obj.d)
            variable delta_var(obj.N*obj.M)
            minimize(obj.C*(alpha_var'*(obj.JM'*obj.JM)*alpha_var/obj.N^2 - 2*obj.PI_star'*obj.JM/obj.N*alpha_var) +...
                norm(u_var) + obj.l*delta_var'*obj.SM'*obj.SM*delta_var)
            A*[alpha_var; u_var; delta_var] <= b
            Aeq*[alpha_var; u_var; delta_var] == beq
            alpha_var >= lb(1:obj.N*obj.M)
            alpha_var <= ub(1:obj.N*obj.M)
            delta_var >= lb(obj.N*obj.M+obj.N*obj.d+1:end) % redundant: deltas are always positive because of the cost
            delta_var <= ub(obj.N*obj.M+obj.N*obj.d+1:end)
            cvx_end
            time_to_solve_qp = toc;
            alpha = alpha_var;
            u = u_var;
            delta = delta_var;
            exit_flag = -99;
            switch cvx_status
                case 'Solved'
                    exit_flag = 1;
                case 'Infeasible'
                    exit_flag = -2;
            end
            lambda = nan;
            cost = obj.C*(alpha'*(obj.JM'*obj.JM)*alpha/obj.N^2 - 2*obj.PI_star'*obj.JM/obj.N*alpha) +...
                norm(u) + obj.l*delta'*obj.SM'*obj.SM*delta;
        end
        function [alpha, u, delta, exit_flag, lambda, time_to_solve_qp] = solve_qp_relax(obj, optimization_data)
            obj.optimization_data = optimization_data;
            H = blkdiag(2*obj.C*(obj.JM'*obj.JM)/obj.N^2, 2*eye(obj.N*obj.d), 2*obj.l*obj.SM);
            f = [-2*obj.C*obj.PI_star'*obj.JM/obj.N, zeros(1,obj.N*obj.d), zeros(1,obj.N*obj.M)]';
            [A,b,Aeq,beq,lb,ub] = obj.get_constraints('centralized');
            tic
            [alpha_u_delta, ~, exit_flag, ~, y] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], optimoptions(@quadprog,'Display','off'));
            time_to_solve_qp = toc;
            alpha = alpha_u_delta(1:obj.N*obj.M);
            u = alpha_u_delta(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.d);
            delta = alpha_u_delta(obj.N*obj.M+obj.N*obj.d+1:end);
            lambda = y.ineqlin;
        end
        function [u, delta, exit_flag, lambda, time_to_solve_qp] = solve_reduced_qp(obj, optimization_data, alpha)
            obj.optimization_data = optimization_data;
            H = blkdiag(2*eye(obj.N*obj.d), 2*obj.l*obj.SM);
            f = [zeros(1,obj.N*obj.d), zeros(1,obj.N*obj.M)]';
            [A,b,Aeq,beq,lb,ub] = obj.get_reduced_constraints(alpha);
            tic
            [u_delta, ~, exit_flag, ~, y] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], optimoptions(@quadprog,'Display','off'));
            time_to_solve_qp = toc;
            u = u_delta(1:obj.N*obj.d);
            delta = u_delta(obj.N*obj.d+1:end);
            lambda = y.ineqlin;
        end
        function [alpha, u, delta, exit_flag, lambda, time_to_solve_qp, ALPHA_VAR] = solve_qp_sdp_relax(obj, optimization_data)
            obj.optimization_data = optimization_data;
            [A,b,Aeq,beq,lb,ub] = obj.get_constraints('centralized');
            tic
            cvx_begin quiet
            variable alpha_var(obj.N*obj.M)
            variable ALPHA_VAR(obj.N*obj.M,obj.N*obj.M)
            variable u_var(obj.N*obj.d)
            variable delta_var(obj.N*obj.M)
            dual variable lambda
            minimize(trace(obj.C*(obj.JM'*obj.JM)/obj.N^2*ALPHA_VAR)-2*obj.C*obj.PI_star'*obj.JM/obj.N*alpha_var+norm(u_var)+obj.l*norm(delta_var))
            lambda : A*[alpha_var; u_var; delta_var] <= b
            Aeq*[alpha_var; u_var; delta_var] == beq
            alpha_var >= lb(1:obj.N*obj.M)
            alpha_var <= ub(1:obj.N*obj.M)
            delta_var >= lb(obj.N*obj.M+obj.N*obj.d+1:end) % not needed: because of the cost, deltas are always positive
            delta_var <= ub(obj.N*obj.M+obj.N*obj.d+1:end)
            diag(ALPHA_VAR) == alpha_var
            [1, alpha_var'; alpha_var, ALPHA_VAR] == semidefinite(obj.N*obj.M+1)
            % strengthening SDP
            % for i = 1 : size(ALPHA_VAR,1)
            %     for j = 1 : size(ALPHA_VAR,2)
            %         ALPHA_VAR(i,j) >= 0
            %         ALPHA_VAR(i,j) >= alpha_var(i)+alpha_var(j)-1
            %         ALPHA_VAR(i,j) <= alpha_var(i)
            %         ALPHA_VAR(i,j) <= alpha_var(j)
            %     end
            % end
            cvx_end
            time_to_solve_qp = toc;
            alpha = alpha_var;
            u = u_var;
            delta = delta_var;
            exit_flag = -99;
            switch cvx_status
                case 'Solved'
                    exit_flag = 1;
                case 'Infeasible'
                    exit_flag = -2;
            end
        end
        function h = get_h(obj)
            h = zeros(obj.N*obj.M,1);
            for i = 1 : obj.N
                for j = 1 : obj.M
                    h((i-1)*obj.M+j) = obj.gamma_of_h_of_x(i,j);
                end
            end
        end
    end
    
    methods (Access = protected)
        function h = gamma_of_h_of_x(obj, i, j)
            switch j
                case 1 % go to point 1
                    J = norm(obj.optimization_data.p(:,i)-obj.optimization_data.p1)^2;
                case 2 % go to point 2
                    J = norm(obj.optimization_data.p(:,i)-obj.optimization_data.p2)^2;
                case 3 % go to point 3
                    J = norm(obj.optimization_data.p(:,i)-obj.optimization_data.p3)^2;
            end
            h = -J;
        end
        function dh = dh_dx(obj, i, j)
            switch j
                case 1 % go to point 1
                    dJ = 2*(obj.optimization_data.p(:,i)-obj.optimization_data.p1);
                case 2 % go to point 2
                    dJ = 2*(obj.optimization_data.p(:,i)-obj.optimization_data.p2);
                case 3 % go to point 3
                    dJ = 2*(obj.optimization_data.p(:,i)-obj.optimization_data.p3);
            end
            dh = -dJ;
        end
        function H = get_cost_quadratic_term_matrix(obj)
            H = blkdiag(2*obj.C*(obj.JM'*obj.JM)/obj.N^2, 2*eye(obj.N*obj.d), 2*obj.l*obj.SM);
        end
        function [A,b,Aeq,beq,lb,ub] = get_constraints(obj, mode, varargin)
            switch mode
                case 'centralized'
                    A = zeros(obj.N*obj.M+obj.N*obj.M^2, 2*obj.N*obj.M+obj.N*obj.d);
                    b = zeros(obj.N*obj.M+obj.N*obj.M^2, 1);
                    Aeq = zeros(obj.N, 2*obj.N*obj.M+obj.N*obj.d);
                    % beq = zeros(obj.N, 1);
                    lb = -inf(2*obj.N*obj.M+obj.N*obj.d, 1);
                    ub = inf(2*obj.N*obj.M+obj.N*obj.d, 1);
                    for i = 1 : obj.N
                        for j = 1 : obj.M
                            % CBFs for tasks
                            A((i-1)*obj.M+j, obj.N*obj.M+(i-1)*obj.d+1:obj.N*obj.M+i*obj.d) = -obj.dh_dx(i,j);
                            b((i-1)*obj.M+j) = 10 * nthroot(obj.gamma_of_h_of_x(i,j),3);
                            % deltas-alphas
                            A(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+1:obj.N*obj.M+(i-1)*obj.M^2+j*obj.M,(i-1)*obj.M+1:i*obj.M) = obj.delta_max*obj.onec(obj.M,j);
                            A(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+1:obj.N*obj.M+(i-1)*obj.M^2+j*obj.M,obj.N*obj.M+obj.N*obj.d+(i-1)*obj.M+1:obj.N*obj.M+obj.N*obj.d+i*obj.M) = -1/obj.K*eye(obj.M) + obj.onec(obj.M,j);
                        end
                        % alpha_i sum up to 1
                        Aeq(i,(i-1)*obj.M+1:i*obj.M) = ones(1,obj.M);
                    end
                    % CBFs for tasks
                    A(1:obj.N*obj.M, obj.N*obj.M+obj.N*obj.d+1:end) = -eye(obj.N*obj.M);
                    % deltas-alphas
                    b(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.M^2) = obj.delta_max*ones(obj.N*obj.M^2,1);
                    % alpha_i sum up to 1
                    beq = ones(obj.N,1);
                    % bounds for alpha
                    lb(1:obj.N*obj.M) = zeros(obj.N*obj.M,1);
                    ub(1:obj.N*obj.M) = ones(obj.N*obj.M,1);
                    % bounds for u
                    lb(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.d) = -0.1*ones(obj.N*obj.d,1);
                    ub(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.d) = 0.1*ones(obj.N*obj.d,1);
                    % bounds for delta: zero is not feasible, but with -delta_max the inequality constraints should be switched around
                    lb(obj.N*obj.M+obj.N*obj.d+1:end) = zeros(obj.N*obj.M,1); % redundant: deltas are always positive because of the cost
                    ub(obj.N*obj.M+obj.N*obj.d+1:end) = obj.delta_max*ones(obj.N*obj.M,1);
                    % remove constraints between a task and itself
                    counter = 0;
                    for i = 1 : obj.N
                        for j = 1 : obj.M
                            for k = 1 : obj.M
                                if j == k
                                    A(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+k-counter,:) = [];
                                    b(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+k-counter) = [];
                                    counter = counter + 1;
                                end
                            end
                        end
                    end
                case 'decentralized'
                    A = nan(obj.M+obj.M^2, 2*obj.M+obj.d);
                    b = nan(obj.M+obj.M^2, 1);
                    Aeq = nan(1, 2*obj.M+obj.d);
                    beq = nan(1, 1);
                    lb = -inf(2*obj.M+obj.d, 1);
                    ub = inf(2*obj.M+obj.d, 1);
                    i = varargin{1};
                otherwise
                    error('Mode not recognized (only ''centralized'' or ''decentralized'' allowed).')
            end
        end
        function [A,b,Aeq,beq,lb,ub] = get_reduced_constraints(obj, alpha)
            A = zeros(obj.N*obj.M+obj.N*obj.M^2, obj.N*obj.M+obj.N*obj.d);
            b = zeros(obj.N*obj.M+obj.N*obj.M^2, 1);
            Aeq = [];
            beq = [];
            lb = -inf(obj.N*obj.M+obj.N*obj.d, 1);
            ub = inf(obj.N*obj.M+obj.N*obj.d, 1);
            for i = 1 : obj.N
                for j = 1 : obj.M
                    % CBFs for tasks
                    A((i-1)*obj.M+j, (i-1)*obj.d+1:i*obj.d) = -obj.dh_dx(i,j);
                    b((i-1)*obj.M+j) = 100*nthroot(obj.gamma_of_h_of_x(i,j),3);
                    % deltas-alphas
                    A(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+1:obj.N*obj.M+(i-1)*obj.M^2+j*obj.M, obj.N*obj.d+(i-1)*obj.M+1:obj.N*obj.d+i*obj.M) = -1/obj.K*eye(obj.M) + obj.onec(obj.M,j);
                    b(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+1:obj.N*obj.M+(i-1)*obj.M^2+j*obj.M) =  -obj.delta_max*obj.onec(obj.M,j)*alpha((i-1)*obj.M+1:i*obj.M);
                end
            end
            % CBFs for tasks
            A(1:obj.N*obj.M, obj.N*obj.d+1:end) = -eye(obj.N*obj.M);
            % deltas-alphas
            b(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.M^2) = b(obj.N*obj.M+1:obj.N*obj.M+obj.N*obj.M^2) + obj.delta_max*ones(obj.N*obj.M^2,1);
            % bounds for u
            lb(1:obj.N*obj.d) = -0.1*ones(obj.N*obj.d,1);
            ub(1:obj.N*obj.d) = 0.1*ones(obj.N*obj.d,1);
            % bounds for delta: zero is not feasible, but with -delta_max the inequality constraints should be switched around
            lb(obj.N*obj.d+1:end) = zeros(obj.N*obj.M,1); % redundant: deltas are always positive because of the cost
            ub(obj.N*obj.d+1:end) = obj.delta_max*ones(obj.N*obj.M,1);
            % remove constraints between a task and itself
            counter = 0;
            for i = 1 : obj.N
                for j = 1 : obj.M
                    for k = 1 : obj.M
                        if j == k
                            A(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+k-counter,:) = [];
                            b(obj.N*obj.M+(i-1)*obj.M^2+(j-1)*obj.M+k-counter) = [];
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    end
    
    methods (Static)
        function m = onec(dim,col_idx)
            m = zeros(dim);
            m(:,col_idx) = ones(dim,1);
        end
    end
end
