classdef grid
    properties
        base_MVA
        base_V
        A
        A_t
        C
        Csh
        Y
        n
        m
        K
        i_limit
        v_limit
        f_upper_limit
        f_lower_limit
        cost_vector_p_g
        cost_vector_q_g
    end
    
    methods
        function obj = grid()
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
            
            obj.K=[0.1 -0.1 0 0.1]'; %stiffness of the primary frequency controller at each generator or load
            
            Gsh = [0 0 0 0 0]; %shunt conductance in MW at V = 1 p.u.
            Bsh = [0 0 0 19 0]; % shunt susceptance in MVar at V = 1 p.u.
            Csh_half= (1/1^2)*(1/base_MVA)*diag(Gsh + 1i*Bsh); % admittance matrix in p.u.
            obj.Csh = [Csh_half,zeros(length(Csh_half));zeros(length(Csh_half)),Csh_half];
            
            %limits
            
            obj.i_limit = 1.06*ones(1,2*obj.m); % line current limit in p.u.
            
            obj.v_limit = 1.06*ones(1,obj.n); %voltage limits at each node in p.u.
            
            obj.f_upper_limit = 0.5; % in Hz
            obj.f_lower_limit = -0.5; % in Hz
            
            obj.cost_vector_p_g = [1.2 1.1 1 1]';
            obj.cost_vector_q_g = [0.1 0.1 0.1 0.2]';
        end
        
        function c = cost_of_generation(obj, mystate)
            
            c = obj.cost_vector_p_g'*mystate.p_g + obj.cost_vector_q_g'*mystate.q_g;
        end
        
        function cD = cost_of_generationD(obj, state)
            
            cD = [zeros(1, 2*obj.n), obj.cost_vector_p_g', obj.cost_vector_q_g', zeros(1,obj.n), 0];
        end
    end
end