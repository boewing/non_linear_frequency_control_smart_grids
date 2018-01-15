classdef Controller < handle
    properties(Constant)
        penalty_factor_i = 0.5
        penalty_factor_v = 5
        penalty_factor_S = 2
        penalty_factor_f = 100
        step_size = 0.25
    end
    
    methods(Static)
        function v_limit = v_limit_reached(x, mygrid)
            v_limit = (0 ~= Controller.getPenaltyVD(x, mygrid)');
        end
        
        function i_limit = i_limit_reached(x,mygrid)
            i_re_im = (0 ~= Controller.getPenaltyID(x, mygrid)');
            i_limit = i_re_im(5*mygrid.n+1:5*mygrid.n+2*mygrid.m) | i_re_im(1+5*mygrid.n+2*mygrid.m:5*mygrid.n+4*mygrid.m);
        end
        
        function f_limit = f_limit_reached(x,mygrid)
            %f_limit = (0 ~= Controller.getPenaltyFD(x, mygrid));
            f_limit = (abs(x.f) > 10e-3);
        end
        
        function S_limit = S_limit_reached(mystate,mygrid)
            S_limit = (mystate.p_g.^2 + mystate.q_g.^2 - mygrid.S_limit.^2 >= -10e-4);
        end
        
        % matrix to extend Aeq in quadprog (the inner optimization)
        % matrix has dimensions n x (5*n + 4*m + 1)
        function [A, b] = S_hard_constraint(x, mygrid)            
%             k = 2;
%             x_circle = mygrid.S_limit(k)*cos(linspace(0,2*pi,1000));
%             y_circle = mygrid.S_limit(k)*sin(linspace(0,2*pi,1000));
%             gamma = mygrid.S_limit(k)/sqrt(x.p_g(k)^2 + x.q_g(k)^2);
%             x_line = gamma*x.p_g(k) + linspace(-3,+3,10).*x.q_g(k);
%             y_line = gamma*x.q_g(k) - linspace(-3,+3,10).*x.p_g(k);
%             figure(3);
%             plot(x_circle, y_circle, x_line, y_line, x.p_g(k)*[1, gamma], x.q_g(k)*[1 gamma]);
%             axis equal;
%             drawnow;
%                 
            b = mygrid.S_limit.*sqrt(x.p_g.^2 + x.q_g.^2) - x.p_g.^2 - x.q_g.^2;
            A = [zeros(mygrid.n,2*mygrid.n), diag(x.p_g), diag(x.q_g), zeros(mygrid.n,mygrid.n + 4*mygrid.m + 1)];
            bounded = (b ~= Inf);
            
            A = A(bounded,:);
            b = b(bounded);
        end
        
        function d = getStep(x, mygrid, k)
            H=eye(5*mygrid.n + 4*mygrid.m + 1);
            ff = Controller.step_size*Controller.n_Jt(x,mygrid);
            [A, b] = Controller.S_hard_constraint(x, mygrid);
            Aeq = Physics.n_h(x,mygrid);
            beq = zeros(3*mygrid.n + 4*mygrid.m,1);
            [lb, ub] = mygrid.bounds(x,k);
            %d = quadprog(H,ff,[],[],Aeq,beq,lb,ub,[],optimoptions('quadprog'));              
            %assert(min(b - A*d) >= 0)                      %condition to
            %check if there is an active constraint
            d = quadprog(H,ff,A,b,Aeq,beq,lb,ub,[],optimoptions('quadprog','Display','off'));
            assert(length(d) == 5*mygrid.n + 4*mygrid.m + 1);
            %     assert(max(Aeq*d) < 1e-10, 'next step is not in the tangent plane');
        end
        
        %%cost function f
        function val = Jt(state, mygrid)
            val = mygrid.cost_of_generation(state)...
                + Controller.getPenaltyI(state, mygrid) ...
                + Controller.getPenaltyV(state, mygrid) ...
                + Controller.getPenaltyF(state, mygrid);
        end
        
        function val = n_Jt(state, mygrid)
            val = mygrid.cost_of_generationD(state) ...
                + Controller.getPenaltyID(state, mygrid) ...
                + Controller.getPenaltyVD(state, mygrid) ...
                + Controller.getPenaltyFD(state, mygrid);
            
        end
        
        function val = getPenaltyS(mystate, mygrid)
            val = Controller.penalty_factor_S*sum(max(sqrt(mystate.p_g.^2 + mystate.q_g.^2) - mygrid.S_limit, 0).^2);
        end
        
        function val = getPenaltySD(mystate, mygrid)
            within_limits = (mystate.p_g.^2 + mystate.q_g.^2 <= mygrid.S_limit);
            if min(within_limits) == 0
                disp( 'S limit reached for at least one generator');
            end
            
            generic_block = (Controller.penalty_factor_S*2*max(1-mygrid.S_limit./sqrt(mystate.p_g.^2 + mystate.q_g.^2),0))';
            pq_D = [generic_block.*mystate.p_g', generic_block.*mystate.q_g'];
            %pq_D = [2*(1-within_limits').*state.p_g', 2*(1-within_limits').*state.q_g'];
            val = [zeros(1,2*mygrid.n), pq_D, zeros(1,mygrid.n + 4*mygrid.m + 1)];
        end
        
        function val = getPenaltyI(mystate, mygrid)
            
            val = Controller.penalty_factor_i*sum(max(sqrt(mystate.i_re.^2 + mystate.i_im.^2) - mygrid.i_limit, 0).^2);
        end
        
        function val = getPenaltyID(mystate, mygrid)
            %off_limits = (mystate.i_re.^2 + mystate.i_im.^2 >= mygrid.i_limit.^2)';
            %i_D_old = Controller.penalty_factor_i*[2*off_limits.*((1 - mygrid.i_limit./sqrt(mystate.i_re.^2 + mystate.i_im.^2)).*mystate.i_re)',... 
            %    2*off_limits.*((1 - mygrid.i_limit./sqrt(mystate.i_re.^2 + mystate.i_im.^2)).*mystate.i_im)'];
            
            generic_block = (Controller.penalty_factor_i*2*max(1-mygrid.i_limit./sqrt(mystate.i_re.^2 + mystate.i_im.^2),0))';
            i_D = [generic_block.*mystate.i_re', generic_block.*mystate.i_im'];
            val = Controller.penalty_factor_i*[zeros(1,5*mygrid.n), i_D, 0];
        end
        
        function val = getPenaltyV(mystate,mygrid)
            val = Controller.penalty_factor_v*sum(max(mystate.v - mygrid.v_limit,0).^2);
        end
        
        function val = getPenaltyVD(mystate,mygrid)
            val = Controller.penalty_factor_v*[2*max(mystate.v - mygrid.v_limit,0)',zeros(1,4*mygrid.n+4*mygrid.m+1)];
        end
        
        function val = getPenaltyF(mystate,mygrid)
            val = Controller.penalty_factor_f*(max(mystate.f - mygrid.f_upper_limit,0)^2 + max(-mystate.f + mygrid.f_lower_limit,0)^2);
        end
        
        function val = getPenaltyFD(mystate,mygrid)
            val = Controller.penalty_factor_f*[zeros(1,5*mygrid.n + 4*mygrid.m), 2*(max(mystate.f - mygrid.f_upper_limit,0) + (-1)*max(-mystate.f + mygrid.f_lower_limit,0))];
        end
        
    end
end