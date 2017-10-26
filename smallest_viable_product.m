function smallest_viable_product()

%%building up the network from constants
mygrid = grid();

%%state initialisation
x = state(mygrid);

%test all functions
h(x,mygrid);
f(x,mygrid);
n_h(x,mygrid);
n_f(x,mygrid);


%calculate the projected gradient
alpha=1;
iterations=10;
for k=1:iterations
    
    %d must be in the tangent plane defined by nabla_h(x), nabla_h(x)*d=0
    %d must be as close to the negative gradient of the cost function:
    %min|d-nabla_f(x)|^2
    d = n_f(x,mygrid)' - n_h(x,mygrid)'*inv(n_h(x,mygrid)*n_h(x,mygrid)')*n_h(x,mygrid)*n_f(x,mygrid)';
    
    norm(h(x,mygrid))
    x.setx(x.getx() + alpha*d);
    
    %retraction
    x = physics(x, mygrid);

end
%%cost function f
    function f = f(state, mygrid)
        f = mygrid.cost_of_generation(state);
           % + getPenaltyI(state, mygrid) ...
           %+ norm(max(state.v - mygrid.v_limit, 0))...
           % + norm(max(state.f - mygrid.f_upper_limit,0)) + norm(max(-state.f + mygrid.f_lower_limit,0));
    end

    function n_f = n_f(state, mygrid)
        n_f = mygrid.cost_of_generationD(state);
    end

    function val = getPenaltyI(mystate, mygrid)
        val = norm(max(mystate.i - mygrid.i_limit, 0));
    end

%%equality constraints function h
    function h = h(mystate,grid)
        u = mystate.v.*exp(1i*mystate.theta);
        temp1 = diag(u)*conj(grid.Y*u);
        temp3 = (1 + grid.K*mystate.f).*mystate.p_ref - mystate.p_g;
        %temp4 = abs((grid.C*grid.A-grid.Csh*grid.A_t)*u).^2 - state.i.^2;
        
        %h=[real(temp1)-state.p_g; imag(temp1)-state.q_g; temp3; temp4];
        h=[real(temp1)-mystate.p_g; imag(temp1)-mystate.q_g; temp3];
    end
    
    %this is the h function with a split in a fixed part and a part to
    %change variables
    %fix: p_g, q_g, p_ref; variable: v, theta, f
    function var = h_fix(mystate, mygrid, x_var)
        n=mygrid.n;
        tempstate = state(mygrid);
        tempstate.v = x_var(1:n);
        tempstate.theta = x_var(n+1:2*n);
        tempstate.f = x_var(end);
        var=h(tempstate,mygrid);
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
        
        line1 = [lu_block, zeros(2*mygrid.n,mygrid.n + 1)];
        line2 = [zeros(mygrid.n,2*mygrid.n),-eye(mygrid.n),zeros(mygrid.n),(eye(mygrid.n)+mystate.f*diag(mygrid.K)), diag(mygrid.K)*mystate.p_ref];
        
        n_h=[line1; line2];
        
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

%%physics calculation from p_ref

    function newstate = physics(mystate, mygrid)
        xx = mystate.getx;
        initial_x_var = [xx(1:mygrid.n); xx(mygrid.n+1:2*mygrid.n); xx(end)];
        new_x_var = fsolve(@(y) h_fix(mystate, mygrid, y), initial_x_var);
        newstate = %... there has to be created a new state object
            
    end

end