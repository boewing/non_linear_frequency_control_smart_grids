classdef State < handle
    properties
        v
        theta
        p_g
        q_g
        p_ref
        i_re
        i_im
        f
    end
    
    methods
        function obj = State(mygrid)
            obj.v = ones(mygrid.n,1);
            obj.theta = zeros(mygrid.n,1);
            
            min_demand = -sum(min(mygrid.p_ref_upper_limit(1),0));
            max_production = sum(max(mygrid.p_ref_upper_limit(1),0));
            obj.p_ref = min_demand/max_production*max(mygrid.p_ref_upper_limit(1),0) + min(mygrid.p_ref_upper_limit(1),0); %if load: set it to the minimal load, if generator: set it such that it has proportional share of the loads
            assert(abs(sum(obj.p_ref)) <= 1e-10);
            obj.p_g = obj.p_ref;
            obj.q_g = zeros(mygrid.n,1);
            obj.i_re = zeros(2*mygrid.m,1);
            obj.i_im = zeros(2*mygrid.m,1);
            obj.f = 0;
            obj = Physics.retraction(obj,mygrid);
            assert(max(abs(Physics.h(obj,mygrid))) < 1e-3, 'retraction failed');
        end
        
        function x = getx(obj)
            x=[obj.v; obj.theta; obj.p_g; obj.q_g; obj.p_ref; obj.i_re; obj.i_im; obj.f];
        end
        
        function setx(obj, x)
            n = length(obj.v);
            m = length(obj.i_re)/2;
            obj.v = x(1:n);
            obj.theta = x(1+n:2*n);
            obj.p_g = x(1+2*n:3*n);
            obj.q_g = x(1+3*n:4*n);
            obj.p_ref = x(1+4*n:5*n);
            obj.i_re = x(1+5*n:5*n+2*m);
            obj.i_im = x(1+5*n+2*m:5*n+4*m);
            obj.f = x(1+5*n+4*m);
        end
    end
end