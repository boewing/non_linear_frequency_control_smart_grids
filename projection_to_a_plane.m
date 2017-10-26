function projection_to_a_plane()

%solve argmin( ||v_t - v||^2 
%           s.t. A*v_t = 0

N=100;
A=rand(N);
v=rand(N,1);

v_t=v-A'*inv(A*A')*A*v;


norm(A*v)
norm(v_t)


end