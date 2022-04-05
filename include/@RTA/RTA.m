classdef RTA < handle
    
    properties %(Access=private)
        scenario_params_
        opt_params_
        dim_
        P_
        constraints_
    end
    
    methods
        % constructor
        function this = RTA(scenario_params, opt_params)
            % scenario_params.A         : Matrix containing feature to robot mapping (size(A) = n_f * n_r)
            % scenario_params.Hs        : Each cell in Hs corresponds to a capability.
            %                             Inside each cell is a matrix of size [c_nhe x n_f] is stored:
            %                                 c_nhe: the number of hyper-edges for capability k
            %                                 n_f  : the number of features.
            % scenario_params.T         : Matrix containg task to capability mapping. (size(T) = n_t * n_c)
            % scenario_params.ws        : Each cell in ws corresponds to a capability.
            %                             Inside each cell is a vector containing the weights for the hyper-edges.
            % scenario_params.robot_dyn : Cell array of function handles for f and g in $\dot x = f(x) + g(x) u$,
            %                             and $n_x$ and $n_u$, where $x\in\mathbb R^{n_x}$ and $u\in\mathbb R^{n_u}$
            % scenario_params.tasks     : Cell array of tasks.(size(tasks) = n_t * 2)
            %                             On row i, first and the second cells contain the function handles
            %                             for h and its gradient, respectively, corresponding to task i.
            %                             If gradient is [], it will be calculated numerically.
            % opt_params.l              : Scaling constant for the cost
            % opt_params.kappa          : Scaling constant for alpha's
            % opt_params.gamma          : Class K function handle
            % opt_params.n_r_bounds     : Matrix of min and max number of robots for each task (size(n_r_bounds) = n_t * 2)
            % opt_params.delta_max      : Max oo-norm for each delta_i
            assert(all(isfield(scenario_params, {'A','Hs','T','ws','robot_dyn','tasks'})), 'Missing scenario parameters.')
            assert(all(isfield(opt_params, {'l','kappa','gamma','n_r_bounds','delta_max'})), 'Missing optimization parameters.')
            this.scenario_params_ = scenario_params;
            this.opt_params_ = opt_params;
            this.dim_.n_r = size(scenario_params.A, 2); % Number of Robots
            this.dim_.n_t = size(scenario_params.T, 1); % Number of Tasks
            this.dim_.n_c = size(scenario_params.T, 2); % Number of Capabilities
            this.dim_.n_f = size(scenario_params.A, 1); % Number of Features
            this.dim_.n_x = this.scenario_params_.robot_dyn.n_x; % State dimension (copied to dim_ property for convenience)
            this.dim_.n_u = this.scenario_params_.robot_dyn.n_u; % Input dimension (copied to dim_ property for convenience)
            this.evaluate_mappings_and_specializations()
            this.check_tasks()
        end
        
        % optimization problem
        function [alpha, u, delta, time_to_solve_miqp, opt_sol_info] = solve_miqp(this, x, t)
            this.build_constraints(x, t)
            alpha_dim = this.dim_.n_r*this.dim_.n_t;
            u_dim = this.dim_.n_r*this.dim_.n_u;
            delta_dim = this.dim_.n_r*this.dim_.n_t;
            tic
            cvx_begin quiet
            binary variable alpha_var(alpha_dim)
            variable u_var(u_dim)
            variable delta_var(delta_dim)
            minimize(1e6*max(1,this.opt_params_.l)*alpha_var'*this.P_'*this.P_*alpha_var + u_var'*u_var + this.opt_params_.l * delta_var'*diag(reshape(this.scenario_params_.S, 1, numel(this.scenario_params_.S)))*delta_var)
            this.constraints_.A_ineq * [alpha_var; u_var; delta_var] <= this.constraints_.b_ineq
            this.constraints_.A_eq * [alpha_var; u_var; delta_var] <= this.constraints_.b_eq % 1^T alpha == 1 is now 1^T alpha <= 1 (from "exactly 1" to "at most 1" task assigned to each robot)
            alpha_var >= this.constraints_.lb(1:this.dim_.n_r*this.dim_.n_t)
            alpha_var <= this.constraints_.ub(1:this.dim_.n_r*this.dim_.n_t)
            delta_var >= this.constraints_.lb(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end)
            delta_var <= this.constraints_.ub(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end)
            cvx_end
            time_to_solve_miqp = toc;
            
            alpha = alpha_var;
            u = u_var;
            delta = delta_var;
            opt_sol_info = cvx_status;
        end
        
        function [u, delta, time_to_solve_qp, opt_sol_info] = solve_reduced_qp(this, x, alpha, t)
            this.build_reduced_constraints(x, t)
            alpha_dim = this.dim_.n_r*this.dim_.n_t;
            u_dim = this.dim_.n_r*this.dim_.n_u;
            delta_dim = this.dim_.n_r*this.dim_.n_t;
            H = blkdiag(2*eye(u_dim), 2*this.opt_params_.l*diag(reshape(this.scenario_params_.S, 1, numel(this.scenario_params_.S))));
            f = zeros(u_dim+delta_dim,1);
            A = this.constraints_.A_ineq(:,alpha_dim+1:end);
            b = this.constraints_.b_ineq - this.constraints_.A_ineq(:,1:alpha_dim)*alpha;
            tic
            [u_delta, ~, opt_sol_info, ~] = quadprog(H,f,A,b,[],[],this.constraints_.lb(alpha_dim+1:end),this.constraints_.ub(alpha_dim+1:end),[],optimoptions(@quadprog,'Display','off'));
            time_to_solve_qp = toc;
            u = u_delta(1:u_dim);
            delta = u_delta(u_dim+1:end);
        end
        
        % setters
        function set_scenario_params(this, scenario_params)
            fn = fieldnames(this.scenario_params_);
            for i = 1 : numel(fn)
                if isfield(scenario_params, fn{i})
                    assert(isfield(this.scenario_params_, fn{i}), [fn{i}, ' is not a field in scenario_params.'])
                    assert(~(strcmp(fn{i},'F')||strcmp(fn{i},'S')), 'Matrices F and S cannot be set (they are automatically evaluated). To update S, use set_specializations')
                    this.scenario_params_.(fn{i}) = scenario_params.(fn{i});
                end
            end
            this.evaluate_mappings_and_specializations()
        end
        
        function set_opt_params(this, opt_params)
            fn = fieldnames(this.opt_params_);
            for i = 1 : numel(fn)
                if isfield(opt_params, fn{i})
                    assert(isfield(this.opt_params_, fn{i}), [fn{i}, ' is not a field in opt_params.'])
                    this.opt_params_.(fn{i}) = opt_params.(fn{i});
                end
            end
        end
        
        function set_specializations(this, S)
            this.scenario_params_.S = S;
            this.build_projector()
        end
        
        % getters
        function S = get_specializations(this)
            S = this.scenario_params_.S;
        end
        
    end
    
    methods (Access=private)
        function build_constraints(this, x, t)
            this.constraints_.A_ineq = zeros(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+2*this.dim_.n_t, 2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u);
            this.constraints_.b_ineq = zeros(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+2*this.dim_.n_t, 1);
            this.constraints_.A_eq = zeros(this.dim_.n_r, 2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u);
            this.constraints_.b_eq = zeros(this.dim_.n_r, 1);
            this.constraints_.lb = -inf(2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u, 1);
            this.constraints_.ub = inf(2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u, 1);
            for i = 1 : this.dim_.n_r
                for j = 1 : this.dim_.n_t
                    % CBFs for tasks
                    this.constraints_.A_ineq((i-1)*this.dim_.n_t+j, this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_u+1:this.dim_.n_r*this.dim_.n_t+i*this.dim_.n_u) = -this.scenario_params_.tasks{j}.gradient(x(:,i),t,i)*this.scenario_params_.robot_dyn.g(x(:,i));
                    this.constraints_.b_ineq((i-1)*this.dim_.n_t+j) = this.scenario_params_.tasks{j}.gradient(x(:,i),t,i)*this.scenario_params_.robot_dyn.f(x(:,i)) + this.scenario_params_.tasks{j}.time_derivative(x(:,i),t,i) + this.opt_params_.gamma(this.scenario_params_.tasks{j}.function(x(:,i),t,i));
                    % deltas-alphas
                    this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+j*this.dim_.n_t,(i-1)*this.dim_.n_t+1:i*this.dim_.n_t) = this.opt_params_.delta_max*this.onec(this.dim_.n_t, j);
                    this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+j*this.dim_.n_t,this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+(i-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+i*this.dim_.n_t) = -1/this.opt_params_.kappa*eye(this.dim_.n_t) + this.onec(this.dim_.n_t, j);
                end
                % alpha_i sum up to 1
                this.constraints_.A_eq(i,(i-1)*this.dim_.n_t+1:i*this.dim_.n_t) = ones(1, this.dim_.n_t);
            end
            % CBFs for tasks
            this.constraints_.A_ineq(1:this.dim_.n_r*this.dim_.n_t, this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = -eye(this.dim_.n_r*this.dim_.n_t);
            % deltas-alphas
            this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2) = this.opt_params_.delta_max*ones(this.dim_.n_r*this.dim_.n_t^2,1);
            % alpha_i sum up to 1
            this.constraints_.b_eq = ones(this.dim_.n_r, 1);
            % bounds for alpha
            this.constraints_.lb(1:this.dim_.n_r*this.dim_.n_t) = zeros(this.dim_.n_r*this.dim_.n_t, 1);
            this.constraints_.ub(1:this.dim_.n_r*this.dim_.n_t) = ones(this.dim_.n_r*this.dim_.n_t, 1);
            % bounds for u
            % this.constraints_.lb(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u) = -0.1*ones(this.dim_.n_r*this.dim_.n_u,1);
            % this.constraints_.ub(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u) = 0.1*ones(this.dim_.n_r*this.dim_.n_u,1);
            % bounds for delta: zero is not feasible, but with -delta_max the inequality constraints should be switched around
            this.constraints_.lb(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = zeros(this.dim_.n_r*this.dim_.n_t, 1);
            this.constraints_.ub(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = this.opt_params_.delta_max*ones(this.dim_.n_r*this.dim_.n_t, 1);
            for j = 1 : this.dim_.n_t
                % F alpha >= T (minimum amount of capabilities for each task)
                this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+(j-1)*this.dim_.n_c+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+j*this.dim_.n_c,j:this.dim_.n_t:this.dim_.n_r*this.dim_.n_t) = -this.scenario_params_.F;
                this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+(j-1)*this.dim_.n_c+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+j*this.dim_.n_c) = -this.scenario_params_.T(j,:)';
                % maximum and minimnum numbers of robots for each task
                this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+j,j:this.dim_.n_t:this.dim_.n_r*this.dim_.n_t) = ones(1,this.dim_.n_r);
                this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+j) = this.opt_params_.n_r_bounds(j,2);
                this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+this.dim_.n_t+j,j:this.dim_.n_t:this.dim_.n_r*this.dim_.n_t) = -ones(1,this.dim_.n_r);
                this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2+this.dim_.n_t*this.dim_.n_c+this.dim_.n_t+j) = -this.opt_params_.n_r_bounds(j,1);
            end
            % remove constraints between a task and itself
            counter = 0;
            for i = 1 : this.dim_.n_r
                for j = 1 : this.dim_.n_t
                    for k = 1 : this.dim_.n_t
                        if j == k
                            this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+k-counter,:) = [];
                            this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+k-counter) = [];
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
        
        function build_reduced_constraints(this, x, t)
            this.constraints_.A_ineq = zeros(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2, 2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u);
            this.constraints_.b_ineq = zeros(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2, 1);
            this.constraints_.A_eq = zeros(this.dim_.n_r, 2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u);
            this.constraints_.b_eq = zeros(this.dim_.n_r, 1);
            this.constraints_.lb = -inf(2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u, 1);
            this.constraints_.ub = inf(2*this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u, 1);
            for i = 1 : this.dim_.n_r
                for j = 1 : this.dim_.n_t
                    % CBFs for tasks
                    this.constraints_.A_ineq((i-1)*this.dim_.n_t+j, this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_u+1:this.dim_.n_r*this.dim_.n_t+i*this.dim_.n_u) = -this.scenario_params_.tasks{j}.gradient(x(:,i),t,i)*this.scenario_params_.robot_dyn.g(x(:,i));
                    this.constraints_.b_ineq((i-1)*this.dim_.n_t+j) = this.scenario_params_.tasks{j}.gradient(x(:,i),t,i)*this.scenario_params_.robot_dyn.f(x(:,i)) + this.scenario_params_.tasks{j}.time_derivative(x(:,i),t,i) + this.opt_params_.gamma(this.scenario_params_.tasks{j}.function(x(:,i),t,i));
                    % deltas-alphas
                    this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+j*this.dim_.n_t,(i-1)*this.dim_.n_t+1:i*this.dim_.n_t) = this.opt_params_.delta_max*this.onec(this.dim_.n_t, j);
                    this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+j*this.dim_.n_t,this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+(i-1)*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+i*this.dim_.n_t) = -1/this.opt_params_.kappa*eye(this.dim_.n_t) + this.onec(this.dim_.n_t, j);
                end
            end
            % CBFs for tasks
            this.constraints_.A_ineq(1:this.dim_.n_r*this.dim_.n_t, this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = -eye(this.dim_.n_r*this.dim_.n_t);
            % deltas-alphas
            this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_t^2) = this.opt_params_.delta_max*ones(this.dim_.n_r*this.dim_.n_t^2,1);
            % bounds for u
            % this.constraints_.lb(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u) = -0.1*ones(this.dim_.n_r*this.dim_.n_u,1);
            % this.constraints_.ub(this.dim_.n_r*this.dim_.n_t+1:this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u) = 0.1*ones(this.dim_.n_r*this.dim_.n_u,1);
            % bounds for delta: zero is not feasible, but with -delta_max the inequality constraints should be switched around
            this.constraints_.lb(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = zeros(this.dim_.n_r*this.dim_.n_t, 1);
            this.constraints_.ub(this.dim_.n_r*this.dim_.n_t+this.dim_.n_r*this.dim_.n_u+1:end) = this.opt_params_.delta_max*ones(this.dim_.n_r*this.dim_.n_t, 1);
            % remove constraints between a task and itself
            counter = 0;
            for i = 1 : this.dim_.n_r
                for j = 1 : this.dim_.n_t
                    for k = 1 : this.dim_.n_t
                        if j == k
                            this.constraints_.A_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+k-counter,:) = [];
                            this.constraints_.b_ineq(this.dim_.n_r*this.dim_.n_t+(i-1)*this.dim_.n_t^2+(j-1)*this.dim_.n_t+k-counter) = [];
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
        
        function evaluate_mappings_and_specializations(this)
            this.scenario_params_.F = zeros(this.dim_.n_c, this.dim_.n_r);
            for k = 1 : this.dim_.n_c
                if ~isempty(this.scenario_params_.ws)
                    W_k = diag(this.scenario_params_.ws{k}); % size(W_k) = n_k * n_k
                    this.scenario_params_.F(k,:) = max(W_k * ((this.scenario_params_.Hs{k} * this.scenario_params_.A) > 0.999),[],1);
                else
                    this.scenario_params_.F(k,:) = max((this.scenario_params_.Hs{k} * this.scenario_params_.A) > 0.999, [], 1); % Use > 0.999 instead of floor for numerical accuracy issues.
                end
            end
            this.scenario_params_.S = double((this.scenario_params_.T * this.scenario_params_.F) > 0.999); % size(S) = n_t * n_r
            this.build_projector()
        end
        
        function build_projector(this)
            this.P_ = repmat(eye(this.dim_.n_t), 1, this.dim_.n_r);
            for i = 1 : this.dim_.n_r
                Si = diag(this.scenario_params_.S(:,i));
                this.P_(:,(i-1)*this.dim_.n_t+1:i*this.dim_.n_t) = this.P_(:,(i-1)*this.dim_.n_t+1:i*this.dim_.n_t) - Si*pinv(Si);
            end
        end
        
        function check_tasks(this)
            for i = 1 : this.dim_.n_t
                if isempty(this.scenario_params_.tasks{i}.gradient)
                    dh_dx_handle = get_dh_dx_handle(i);
                    this.scenario_params_.tasks{i}.gradient = dh_dx_handle;
                end
                if isempty(this.scenario_params_.tasks{i}.time_derivative)
                    dh_dt_handle = get_dh_dt_handle(i);
                    this.scenario_params_.tasks{i}.time_derivative = dh_dt_handle;
                end
            end
            function dh_dx_handle = get_dh_dx_handle(task_idx)
                dh_dx_handle = @dh_dx;
                function dh_dx_value = dh_dx(x_value, t_value, i)
                    n = numel(x_value);
                    dh_dx_value = zeros(1,n);
                    for j = 1 : n
                        ej = zeros(n,1);
                        ej(j) = 1;
                        dh_dx_value(j) = this.scenario_params_.tasks{task_idx}.function(x_value+1e-3*ej, t_value, i)-this.scenario_params_.tasks{task_idx}.function(x_value-1e-3*ej, t_value, i);
                    end
                    dh_dx_value = dh_dx_value / 2e-3;
                end
            end
            function dh_dt_handle = get_dh_dt_handle(task_idx)
                dh_dt_handle = @dh_dt;
                function dh_dt_value = dh_dt(x_value, t_value, i)
                    dh_dt_value = (this.scenario_params_.tasks{task_idx}.function(x_value, t_value+1e-3, i)-this.scenario_params_.tasks{task_idx}.function(x_value, t_value-1e-3, i)) / 2e-3;
                end
            end
        end
    end
    
    methods (Access=private, Static)
        function m = onec(dim, col_idx)
            m = zeros(dim);
            m(:,col_idx) = ones(dim,1);
        end
    end
end

