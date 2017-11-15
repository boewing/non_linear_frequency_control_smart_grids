function smallest_viable_product()

close
%%building up the network from constants
mygrid = grid();

%%state initialisation
x = state(mygrid);

%calculate the projected gradient
alpha=0.2;
iterations=500;
myrec = Recorder(x, Jt(x,mygrid),norm(h(x,mygrid)),iterations);
for k=1:iterations
    
    %calculate a feasible gradient
    %d must be in the tangent plane defined by nabla_h(x), nabla_h(x)*d=0
    %d must be as close to the negative gradient of the cost function as possible:
    %min|d-nabla_Jt(x)|^2
    %analytic version without constraints: d = n_Jt(x,mygrid)' - n_h(x,mygrid)'*inv(n_h(x,mygrid)*n_h(x,mygrid)')*n_h(x,mygrid)*n_Jt(x,mygrid)';
    
    H=eye(5*mygrid.n + 2*mygrid.m + 1);
    ff = n_Jt(x,mygrid);
    Aeq = n_h(x,mygrid);
    beq = zeros(size(Aeq,1),1);
    [lb, ub] = bounds(x,mygrid);
    d = quadprog(H,ff,[],[],Aeq,beq,lb,ub,[],optimoptions('quadprog','Display','off'));
    
%     assert(max(Aeq*d) < 1e-10, 'next step is not in the tangent plane');
    
    %make a step
    x.setx(x.getx() + alpha*d);
    
    %retraction
    assert(x.theta(1)==0, 'reference bus has non-zero voltage angle');              %check that the angle reference stays 0
    x = retraction(x, mygrid);
%     assert(max(abs(h(x,mygrid))) < 1e-3, 'retraction failed');    %check that the retraction works makes h(x)=0
    
    %angle correction
    ctrl_angle_correction(x);
    
    %recording
    myrec.store(x,Jt(x,mygrid),norm(h(x,mygrid)));
    
end

myrec.plotAll();

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
        val = mygrid.penalty_factor_f*[zeros(1,5*mygrid.n + 2*mygrid.m), 2*(max(mystate.f - mygrid.f_upper_limit,0) + max(-mystate.f + mygrid.f_lower_limit,0))];
    end

%%equality constraints function h
    function h = h(mystate,mygrid)
        u = mystate.v.*exp(1i*mystate.theta);
        temp1 = diag(u)*conj(mygrid.Y*u);
        temp3 = (1 + mygrid.K*mystate.f).*mystate.p_ref - mystate.p_g;
        temp4 = abs((mygrid.C*mygrid.A-mygrid.Csh*mygrid.A_t)*u).^2 - mystate.i.^2;
        
        h=[real(temp1)-mystate.p_g; imag(temp1)-mystate.q_g; temp3; temp4];
        %h=[real(temp1)-mystate.p_g; imag(temp1)-mystate.q_g; temp3];
    end

%this is the h function with a split in a fixed part and a part to
%change variables
%fix: p_ref ; variable: v, theta(2:n), p_g, i, and f
    function var = h_fix(mystate, mygrid, x_var)
        n=mygrid.n;
        m=mygrid.m;
        tempstate = mystate;
        tempstate.v = x_var(1:n);
        assert(tempstate.theta(1)==0);
        tempstate.theta = [0;x_var(1+n:2*n-1)];
        tempstate.p_g = x_var(1+2*n-1:3*n-1);
        tempstate.i = x_var(1+3*n-1:2*m+3*n-1);
        tempstate.f = x_var(2*m+3*n);
        var = h(tempstate,mygrid);
    end
%%nabla h
    function n_h = n_h(mystate,mygrid) % has a reduced form: contains only the equations for p_g, q_g and p_ref
        V0 = mystate.v;
        A0 = mystate.theta;
        
        J0 = mygrid.Y * (V0.*exp(1j*A0));
        UU = bracket(diag(V0.*exp(1j*A0)));
        JJ = bracket(diag(conj(J0)));
        NN = Nmatrix(2*mygrid.n);
        YY = bracket(mygrid.Y);
        PP = Rmatrix(mystate.v, mystate.theta);
        lu_block = [(JJ + UU*NN*YY)*PP, -eye(2*mygrid.n)];
        
        Q = mygrid.C*mygrid.A - mygrid.Csh*mygrid.A_t;
        u = V0.*exp(1j*A0);
        
        line1 = [lu_block, zeros(2*mygrid.n,mygrid.n),zeros(2*mygrid.n,2*mygrid.m),zeros(2*mygrid.n,1)];
        line2 = [zeros(mygrid.n,2*mygrid.n),-eye(mygrid.n),zeros(mygrid.n),(eye(mygrid.n)+mystate.f*diag(mygrid.K)), zeros(mygrid.n,2*mygrid.m),diag(mygrid.K)*mystate.p_ref];
        line3 = [2*real(diag(conj(Q*u))*Q*diag(exp(1j*A0))), 2*imag(diag(conj(Q*u))*Q*diag(u)), zeros(2*mygrid.m,3*mygrid.n), -2*diag(mystate.i), zeros(2*mygrid.m,1)];
        
        n_h=[line1; line2; line3];
        
    end

    function R = Rmatrix(V,TH)
        VV = diag(V);
        COSTH = diag(cos(TH));
        SINTH = diag(sin(TH));
        R = [COSTH, -VV*SINTH; SINTH, VV*COSTH];
    end

    function BRX = bracket(X)
        BRX = [real(X), -imag(X); imag(X), real(X)];
    end

    function N = Nmatrix(d)
        N = [eye(d/2), zeros(d/2); zeros(d/2),  -eye(d/2)];
    end

%%retraction to a physical meaningful state by adjusting v, theta(2:n), p_g,
%%q_g, i, and f such that ||h(x)|| = 0 again.
%v:       n    vars
%theta: n-1    vars
%p_g:     n    vars
%i:     2*m    vars
%f:       1    var
%sum:  3*n+2*m vars

    function newstate = retraction(mystate, mygrid)
        
        xx = mystate.getx;
        %                   v               theta(2:n)                   p_g
        initial_x_var = [xx(1:mygrid.n); xx(2+mygrid.n:2*mygrid.n); xx(1+2*mygrid.n:3*mygrid.n); xx(1+5*mygrid.n:2*mygrid.m+5*mygrid.n);xx(2*mygrid.m+5*mygrid.n + 1)];
        assert(0==norm(h(mystate,mygrid) - h_fix(mystate,mygrid, initial_x_var)));
        new_x_var = fsolve(@(y) h_fix(mystate, mygrid, y), initial_x_var,optimoptions('fsolve','Display', 'off', 'FunctionTolerance',1e-8));
        
        mystate.v = new_x_var(1:mygrid.n);
        mystate.theta = [xx(1+mygrid.n); new_x_var(mygrid.n+1:2*mygrid.n-1)];
        mystate.p_g = new_x_var(1+2*mygrid.n-1:3*mygrid.n-1);
        mystate.i = new_x_var(1+3*mygrid.n-1:2*mygrid.m+3*mygrid.n-1);
        mystate.f = new_x_var(1+2*mygrid.m+3*mygrid.n-1);
        newstate = mystate;
        
    end

    function ctrl_angle_correction(mystate)
        % correct voltage angles to lie in [-pi, pi]
        mystate.theta = mod(mystate.theta + pi/2, pi) - pi/2;
    end

end