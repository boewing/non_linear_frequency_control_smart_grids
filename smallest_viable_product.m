function smallest_viable_product()

close
%%building up the network from constants
mygrid = grid();

%%state initialisation
x = state(mygrid);

%calculate the projected gradient
alpha=0.2;
iterations=500;
myrec = Recorder(x, Jt(x,mygrid),norm(Physics.h(x,mygrid)),iterations);
for k=1:iterations
    
    %calculate a feasible gradient
    %d must be in the tangent plane defined by nabla_h(x), nabla_h(x)*d=0
    %d must be as close to the negative gradient of the cost function as possible:
    %min|d-nabla_Jt(x)|^2
    %analytic version without constraints: d = n_Jt(x,mygrid)' - n_h(x,mygrid)'*inv(n_h(x,mygrid)*n_h(x,mygrid)')*n_h(x,mygrid)*n_Jt(x,mygrid)';
    
    H=eye(5*mygrid.n + 2*mygrid.m + 1);
    ff = n_Jt(x,mygrid);
    Aeq = Physics.n_h(x,mygrid);
    beq = zeros(size(Aeq,1),1);
    [lb, ub] = bounds(x,mygrid);
    d = quadprog(H,ff,[],[],Aeq,beq,lb,ub,[],optimoptions('quadprog','Display','off'));
    
%     assert(max(Aeq*d) < 1e-10, 'next step is not in the tangent plane');
    
    %make a step
    x.setx(x.getx() + alpha*d);
    
    %retraction
    assert(x.theta(1)==0, 'reference bus has non-zero voltage angle');              %check that the angle reference stays 0
    x = Physics.retraction(x, mygrid);
%     assert(max(abs(h(x,mygrid))) < 1e-3, 'retraction failed');    %check that the retraction works makes h(x)=0
    
    %recording
    myrec.store(x,Jt(x,mygrid),norm(Physics.h(x,mygrid)));
    
end

myrec.plotAll();

end

    function [lb, ub] = bounds(state,mygrid)

    v_max_rel = Inf*ones(mygrid.n,1);
    v_min_rel = - state.v;
%    v_min_rel = -Inf*ones(mygrid.n,1);
    
    theta_max_rel = [0; Inf*ones(mygrid.n-1,1)];
    theta_min_rel = [0; -Inf*ones(mygrid.n-1,1)];
    
    p_g_max_rel = Inf*ones(mygrid.n,1);
    p_g_min_rel = -Inf*ones(mygrid.n,1);
    
    q_g_max_rel = Inf*ones(mygrid.n,1);
    q_g_min_rel = -Inf*ones(mygrid.n,1);
    
    p_ref_max_rel = mygrid.p_ref_upper_limit - state.p_ref;
    p_ref_min_rel = mygrid.p_ref_lower_limit - state.p_ref;
    
%     p_ref_max_rel = Inf*abs((mygrid.p_ref_upper_limit - state.p_ref));
%     p_ref_min_rel = -Inf*abs((mygrid.p_ref_lower_limit - state.p_ref));
    
    i_max_rel = Inf*ones(2*mygrid.m,1);
    i_min_rel = - state.i;
%    i_min_rel = - Inf*ones(2*mygrid.m,1);
    
    f_max_rel = Inf;
    f_min_rel = -Inf;
        
    ub = [v_max_rel; theta_max_rel; p_g_max_rel; q_g_max_rel; p_ref_max_rel; i_max_rel; f_max_rel];
    lb = [v_min_rel; theta_min_rel; p_g_min_rel; q_g_min_rel; p_ref_min_rel; i_min_rel; f_min_rel];
    
    end

%%cost function f
    function val = Jt(state, mygrid)
        val = mygrid.cost_of_generation(state)...
        + getPenaltyI(state, mygrid) ...
        + getPenaltyV(state, mygrid) ...
        + getPenaltyF(state, mygrid) ...
        + getPenaltyS(state, mygrid);
    end

    function val = n_Jt(state, mygrid)
        val = mygrid.cost_of_generationD(state) ...
            + getPenaltyID(state, mygrid) ...
            + getPenaltyVD(state, mygrid) ...
            + getPenaltyFD(state, mygrid) ...
            + getPenaltySD(state, mygrid);
            
    end

    function val = getPenaltyS(state, mygrid)
        val = mygrid.penalty_factor_S*sum(max(state.p_g.^2 + state.q_g.^2 - mygrid.S_limit.^2,0));
    end
    
    function val = getPenaltySD(state, mygrid)
        within_limits = (state.p_g.^2 + state.q_g.^2 <= mygrid.S_limit);
        if min(within_limits) == 0
            disp( 'S limit reached for at least one generator');
        end
        pq_D = [2*(1-within_limits').*state.p_g', 2*(1-within_limits').*state.q_g'];
        val = mygrid.penalty_factor_S*[zeros(1,2*mygrid.n), pq_D, zeros(1,mygrid.n + 2*mygrid.m + 1)];
    end

    function val = getPenaltyI(mystate, mygrid)

        val = mygrid.penalty_factor_i*sum(max(mystate.i - mygrid.i_limit, 0).^2);
    end

    function val = getPenaltyID(mystate, mygrid)
        val = mygrid.penalty_factor_i*[zeros(1,5*mygrid.n), 2*max(mystate.i - mygrid.i_limit, 0)', 0];
    end

    function val = getPenaltyV(mystate,mygrid)
        val = mygrid.penalty_factor_v*sum(max(mystate.v - mygrid.v_limit,0).^2);
    end
        
    function val = getPenaltyVD(mystate,mygrid)
        val = mygrid.penalty_factor_v*[2*max(mystate.v - mygrid.v_limit,0)',zeros(1,4*mygrid.n+2*mygrid.m+1)];
    end

    function val = getPenaltyF(mystate,mygrid)
        val = mygrid.penalty_factor_f*(max(mystate.f - mygrid.f_upper_limit,0)^2 + max(-mystate.f + mygrid.f_lower_limit,0)^2);
    end

    function val = getPenaltyFD(mystate,mygrid)
        val = mygrid.penalty_factor_f*[zeros(1,5*mygrid.n + 2*mygrid.m), 2*(max(mystate.f - mygrid.f_upper_limit,0) + (-1)*max(-mystate.f + mygrid.f_lower_limit,0))];
    end
