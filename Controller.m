classdef Controller < handle
    properties(Constant)
        penalty_factor_i = 10
        penalty_factor_v = 5
        penalty_factor_S = 2
        penalty_factor_f = 3
        step_size = 0.15
    end
    
    methods(Static)
        function v_limit = v_limit_reached(x, mygrid)
            v_limit = (0 ~= Controller.getPenaltyVD(x, mygrid)');
        end
        
        function i_limit = i_limit_reached(x,mygrid)
            i_limit = (0 ~= Controller.getPenaltyID(x, mygrid)');
        end
        
        function f_limit = f_limit_reached(x,mygrid)
            %f_limit = (0 ~= Controller.getPenaltyFD(x, mygrid));
            f_limit = (abs(x.f) > 10e-3);
        end
        
        function S_limit = S_limit_reached(x,mygrid)
            S_limit = (0 ~= Controller.getPenaltySD(x, mygrid)');
        end
        
        function d = getStep(x, mygrid, k)
            H=eye(5*mygrid.n + 2*mygrid.m + 1);
            ff = Controller.step_size*Controller.n_Jt(x,mygrid);
            Aeq = Physics.n_h(x,mygrid);
            beq = zeros(size(Aeq,1),1);
            [lb, ub] = mygrid.bounds(x,k);
            d = quadprog(H,ff,[],[],Aeq,beq,lb,ub,[],optimoptions('quadprog','Display','off'));  
            %     assert(max(Aeq*d) < 1e-10, 'next step is not in the tangent plane');
        end
        
        %%cost function f
        function val = Jt(state, mygrid)
            val = mygrid.cost_of_generation(state)...
                + Controller.getPenaltyI(state, mygrid) ...
                + Controller.getPenaltyV(state, mygrid) ...
                + Controller.getPenaltyF(state, mygrid) ...
                + Controller.getPenaltyS(state, mygrid);
        end
        
        function val = n_Jt(state, mygrid)
            val = mygrid.cost_of_generationD(state) ...
                + Controller.getPenaltyID(state, mygrid) ...
                + Controller.getPenaltyVD(state, mygrid) ...
                + Controller.getPenaltyFD(state, mygrid) ...
                + Controller.getPenaltySD(state, mygrid);
            
        end
        
        function val = getPenaltyS(state, mygrid)
            val = Controller.penalty_factor_S*sum(max(state.p_g.^2 + state.q_g.^2 - mygrid.S_limit.^2,0));
        end
        
        function val = getPenaltySD(state, mygrid)
            within_limits = (state.p_g.^2 + state.q_g.^2 <= mygrid.S_limit);
            if min(within_limits) == 0
                disp( 'S limit reached for at least one generator');
            end
            pq_D = [2*(1-within_limits').*state.p_g', 2*(1-within_limits').*state.q_g'];
            val = Controller.penalty_factor_S*[zeros(1,2*mygrid.n), pq_D, zeros(1,mygrid.n + 2*mygrid.m + 1)];
        end
        
        function val = getPenaltyI(mystate, mygrid)
            
            val = Controller.penalty_factor_i*sum(max(mystate.i - mygrid.i_limit, 0).^2);
        end
        
        function val = getPenaltyID(mystate, mygrid)
            val = Controller.penalty_factor_i*[zeros(1,5*mygrid.n), 2*max(mystate.i - mygrid.i_limit, 0)', 0];
        end
        
        function val = getPenaltyV(mystate,mygrid)
            val = Controller.penalty_factor_v*sum(max(mystate.v - mygrid.v_limit,0).^2);
        end
        
        function val = getPenaltyVD(mystate,mygrid)
            val = Controller.penalty_factor_v*[2*max(mystate.v - mygrid.v_limit,0)',zeros(1,4*mygrid.n+2*mygrid.m+1)];
        end
        
        function val = getPenaltyF(mystate,mygrid)
            val = Controller.penalty_factor_f*(max(mystate.f - mygrid.f_upper_limit,0)^2 + max(-mystate.f + mygrid.f_lower_limit,0)^2);
        end
        
        function val = getPenaltyFD(mystate,mygrid)
            val = Controller.penalty_factor_f*[zeros(1,5*mygrid.n + 2*mygrid.m), 2*(max(mystate.f - mygrid.f_upper_limit,0) + (-1)*max(-mystate.f + mygrid.f_lower_limit,0))];
        end
        
    end
end