classdef Grid
    properties
        base_MVA
        base_V
        A
        A_t
        C
        Csh
        Q
        Y
        n
        m
        K
        i_limit
        penalty_factor_i
        v_limit
        penalty_factor_v
        p_ref_upper_limit_base
        p_ref_lower_limit_base
        penalty_factor_S
        S_limit
        f_upper_limit
        f_lower_limit
        penalty_factor_f
        cost_vector_p_g
        cost_vector_q_g
        E_short
    end
    
    methods
        function obj = Grid()
            %definition of p.u. base units
            base_MVA = 100;
            base_MV = 0.132;
            A_onedirec=[-1 1 0 0;-1 0 1 0; 0 -1 1 0; -1 0 0 1; 0 -1 0 1];
            obj.n = size(A_onedirec,2); %number of nodes
            obj.m = size(A_onedirec,1); %number of lines
            A = [A_onedirec;-A_onedirec];
            A_t = (A-A.^2)/2;
            R_half=[0.0192 0.0452 0.057 0.0132 0.0472]; %line resistance in p.u.
            X_half=[0.0575 0.1652 0.1737 0.0379 0.1983]; %line reactance in p.u.
            C_half=diag(1./(R_half + 1i*X_half));
            C = [C_half,zeros(length(C_half));zeros(length(C_half)),C_half];
            obj.Y=A'*C*A;
            obj.C = C;
            obj.A = A;
            obj.A_t = A_t;
            
            obj.K=[-1.0 -1.0 0 0.01]'; %stiffness of the primary frequency controller at each generator or load in p.u./Hz
            
            Gsh = [0 0 0 0 0]; %shunt conductance in MW at V = 1 p.u.
            Bsh = [0 0 0 19 0]; % shunt susceptance in MVar at V = 1 p.u.
            Csh_half= (1/1^2)*(1/base_MVA)*diag(Gsh + 1i*Bsh); % admittance matrix in p.u.
            obj.Csh = [Csh_half,zeros(length(Csh_half));zeros(length(Csh_half)),Csh_half];
            
            obj.Q = obj.C*obj.A-obj.Csh*obj.A_t;
            
            %limits
            
            obj.i_limit = 1.06*ones(2*obj.m,1); % line current limit in p.u.
            obj.penalty_factor_i = 10;
            
            obj.v_limit = 1.06*ones(obj.n,1); %voltage limits at each node in p.u.
            obj.penalty_factor_v = 5;
            
            obj.p_ref_upper_limit_base = [ 1.0;  0.6; -0.2; -0.36];
            obj.p_ref_lower_limit_base = [   0;    0; -0.3; -0.36];
            obj.penalty_factor_S =  2;
            obj.S_limit =           [ 1.0;  1.0;  Inf;   Inf];
            
            obj.f_upper_limit = 0.0; % in Hz
            obj.f_lower_limit = -0.0; % in Hz
            obj.penalty_factor_f = 20.5;
            
            obj.cost_vector_p_g = [  2 1.1  0   0]';
            obj.cost_vector_q_g = [0.3 0.3  0   0]';
            
            %this is a selector matrix selecting the parts of the state
            %which are fixed
            E = diag([zeros(1,obj.n),...   %v
                [1,zeros(1,obj.n-1)],...     %theta
                zeros(1,obj.n),...           %p_g
                ones(1,obj.n),...            %q_g
                ones(1,obj.n),...            %p_ref
                zeros(1,2*obj.m),...         %i
                0]);                            %f
            
            obj.E_short = E(logical(sum(E)'),:);
            assert(min(sum(obj.E_short'))==1);            
        end
        
        %marginal cost of generation
        function c = cost_of_generation(obj, mystate)
            
            c = obj.cost_vector_p_g'*mystate.p_g + obj.cost_vector_q_g'*mystate.q_g.^2;
        end
        
        function cD = cost_of_generationD(obj, mystate)
            
            cD = [zeros(1, 2*obj.n), obj.cost_vector_p_g', 2*(obj.cost_vector_q_g.*mystate.q_g)', zeros(1,obj.n), zeros(1,2*obj.m), 0];
        end
        
        function upper = p_ref_upper_limit(obj,time)
            upper = min(obj.p_ref_upper_limit_base,0) + max(obj.p_ref_upper_limit_base,0);
        end
        function lower = p_ref_lower_limit(obj,time)
            lower = min(obj.p_ref_lower_limit_base,0) + max(obj.p_ref_lower_limit_base,0);
        end
        
        function [lb, ub] = bounds(obj, state, time)
            
            v_max_rel = Inf*ones(obj.n,1);
            v_min_rel = - state.v;
            %    v_min_rel = -Inf*ones(mygrid.n,1);
            
            theta_max_rel = [0; Inf*ones(obj.n-1,1)];
            theta_min_rel = [0; -Inf*ones(obj.n-1,1)];
            
            p_g_max_rel = Inf*ones(obj.n,1);
            p_g_min_rel = -Inf*ones(obj.n,1);
            
            q_g_max_rel = Inf*ones(obj.n,1);
            q_g_min_rel = -Inf*ones(obj.n,1);
            
            p_ref_max_rel = obj.p_ref_upper_limit(time) - state.p_ref;
            p_ref_min_rel = obj.p_ref_lower_limit(time) - state.p_ref;
            
            %     p_ref_max_rel = Inf*abs((obj.p_ref_upper_limit - state.p_ref));
            %     p_ref_min_rel = -Inf*abs((obj.p_ref_lower_limit - state.p_ref));
            
            i_max_rel = Inf*ones(2*obj.m,1);
            i_min_rel = - state.i;
            %    i_min_rel = - Inf*ones(2*obj.m,1);
            
            f_max_rel = Inf;
            f_min_rel = -Inf;
            
            ub = [v_max_rel; theta_max_rel; p_g_max_rel; q_g_max_rel; p_ref_max_rel; i_max_rel; f_max_rel];
            lb = [v_min_rel; theta_min_rel; p_g_min_rel; q_g_min_rel; p_ref_min_rel; i_min_rel; f_min_rel];
            
        end

    end
end