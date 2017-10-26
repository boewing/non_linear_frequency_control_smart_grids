classdef state < handle
   properties
      v
      theta
      p_g
      q_g
      p_ref
      i       
      f  
   end
   
   methods
       function obj = state(mygrid)
%            obj.v = ones(mygrid.n,1);
%            obj.theta = zeros(mygrid.n,1);
           obj.v = [1.02, 0.99, 0.98, 1]';
           obj.theta = 2*pi/360*[0, 0, 0 ,0]';
           
           u = obj.v.*exp(1i*obj.theta);
           temp1 = diag(u)*conj(mygrid.Y*u);
           obj.p_g = real(temp1);
           obj.q_g = imag(temp1);
           obj.p_ref = obj.p_g;
           obj.i = abs((mygrid.C*mygrid.A - mygrid.Csh*mygrid.A_t)*u);
           obj.f = 50;
       end
       
       function x = getx(obj)
           x=[obj.v; obj.theta; obj.p_g; obj.q_g; obj.p_ref; obj.f];
       end
       
       function setx(obj, x)
           n = length(obj.v);
           obj.v = x(1:n);
           obj.theta = x(1+n:2*n);
           obj.p_g = x(1+2*n:3*n);
           obj.q_g = x(1+3*n:4*n);
           obj.p_ref = x(1+4*n:5*n);
           obj.f = x(end);
       end
   end
end