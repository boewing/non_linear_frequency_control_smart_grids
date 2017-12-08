classdef Physics < handle
    %this class is able to calculate the system response according to the
    %controller input and the grid properties. more specific it is able to
    %bring an infeasible state back to a feasible one by retraction.
    
    methods(Static)
        
        %%equality constraints function h
        function h = h(mystate,mygrid)
            if isa(mystate,'State')
                u = mystate.v.*exp(1i*mystate.theta);
                temp1 = diag(u)*conj(mygrid.Y*u);
                temp3 = (mygrid.K*mystate.f + mystate.p_ref) - mystate.p_g;
                temp4 = abs(mygrid.Q*u).^2 - mystate.i.^2;
                
                h=[real(temp1)-mystate.p_g; imag(temp1)-mystate.q_g; temp3; temp4];
                %h=[real(temp1)-mystate.p_g; imag(temp1)-mystate.q_g; temp3];
            else
                n=mygrid.n;
                m=mygrid.m;
                u = mystate(1:n).*exp(1i*mystate(1+n:2*n));
                temp1 = u.*conj(mygrid.Y*u);
                temp3 = mygrid.K*mystate(2*m+5*n+1) + mystate(4*n+1:5*n) - mystate(2*n+1:3*n);
                temp4 = abs(mygrid.Q*u).^2 - mystate(5*n+1:5*n+2*m).^2;
                
                h=[real(temp1)-mystate(1+2*n:3*n); imag(temp1)-mystate(1+3*n:4*n); temp3; temp4];
            end
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
            var = Physics.h(tempstate,mygrid);
        end
        
        %%nabla h
        function n_h = n_h(mystate,mygrid) % has a reduced form: contains only the equations for p_g, q_g and p_ref
            V0 = mystate.v;
            A0 = mystate.theta;
            u = V0.*exp(1j*A0);
            
            J0 = mygrid.Y * (u);
            UU = bracket(diag(u));
            JJ = bracket(diag(conj(J0)));
            NN = Nmatrix(2*mygrid.n);
            YY = bracket(mygrid.Y);
            PP = Rmatrix(mystate.v, mystate.theta);
            lu_block = [(JJ + UU*NN*YY)*PP, -eye(2*mygrid.n)];
            
            Q = mygrid.Q;
            
            line1 = [lu_block, zeros(2*mygrid.n,mygrid.n),zeros(2*mygrid.n,2*mygrid.m),zeros(2*mygrid.n,1)];
            line2 = [zeros(mygrid.n,2*mygrid.n), -eye(mygrid.n), zeros(mygrid.n), eye(mygrid.n), zeros(mygrid.n,2*mygrid.m), mygrid.K];
            line3 = [2*real(diag(conj(Q*u))*Q*diag(exp(1j*A0))), 2*imag(diag(conj(Q*u))*Q*diag(u)), zeros(2*mygrid.m,3*mygrid.n), -2*diag(mystate.i), zeros(2*mygrid.m,1)];
            
            n_h=[line1; line2; line3];
            
            
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
            %fsolve(@(y) Physics.h(mystate, mygrid
            xx = mystate.getx;
            %                   v               theta(2:n)                   p_g
            %initial_x_var = [xx(1:mygrid.n); xx(2+mygrid.n:2*mygrid.n); xx(1+2*mygrid.n:3*mygrid.n); xx(1+5*mygrid.n:2*mygrid.m+5*mygrid.n);xx(2*mygrid.m+5*mygrid.n + 1)];
            %assert(0==norm(Physics.h(mystate,mygrid) - Physics.h_fix(mystate,mygrid, initial_x_var)));
            %rank_of_gradient = min(abs(eig([Physics.n_h(mystate,mygrid);mygrid.E_short])));
            %condition = norm([Physics.n_h(mystate,mygrid);mygrid.E_short])*norm(inv([Physics.n_h(mystate,mygrid);mygrid.E_short]))
            new_x_var = fsolve(@(y) Physics.h_aug(mystate, mygrid, y), xx, optimoptions('fsolve','Algorithm','levenberg-marquardt','Display', 'off', 'FunctionTolerance',1e-10));
            %new_x_var = fsolve(@(y) Physics.h_fix(mystate, mygrid, y), initial_x_var,optimoptions('fsolve','Display', 'off', 'FunctionTolerance',1e-8));
            %new_x_var = lsqnonlin(@(y) Physics.h_fix(mystate, mygrid, y), initial_x_var, [], [], optimoptions('lsqnonlin','Display', 'off', 'FunctionTolerance',1e-8));
            
            mystate.v      = new_x_var(1:mygrid.n);
            mystate.theta  = new_x_var(1+mygrid.n:2*mygrid.n);
            mystate.p_g    = new_x_var(1+2*mygrid.n:3*mygrid.n);
            mystate.q_g    = new_x_var(1+3*mygrid.n:4*mygrid.n);
            mystate.p_ref  = new_x_var(1+4*mygrid.n:5*mygrid.n);
            mystate.i      = abs(new_x_var(1+5*mygrid.n:2*mygrid.m+5*mygrid.n));    %because for the equation h_fix the sign of h is irrelevant
            mystate.f      = new_x_var(1+2*mygrid.m+5*mygrid.n);
            
            Physics.ctrl_angle_correction(mystate);
            newstate = mystate;
            
            %             temp = Physics.h(newstate, mygrid);
            %             newstate.i = -newstate.i;
            %             assert(norm(Physics.h(newstate, mygrid) - temp) == 0);
            
        end
        
        function ctrl_angle_correction(mystate)
            % correct voltage angles to lie in [-pi, pi]
            mystate.theta = mod(mystate.theta + pi/2, pi) - pi/2;
            mystate.theta = mystate.theta - mystate.theta(1);
        end
        
        function val = h_aug(mystate, mygrid, x_var)
            assert(isa(mystate,'State'));
            
            val = [Physics.h(x_var,mygrid); mygrid.E*(mystate.getx - x_var)];
        end
        
        
    end
end