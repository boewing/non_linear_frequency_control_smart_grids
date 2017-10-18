function smallest_viable_product()

%%building up the network from constants
A_onedirec=[-1 1 0 0;-1 0 1 0; 0 -1 1 0; -1 0 0 1; 0 -1 0 1];
grid.n = size(A_onedirec,2);
grid.m = size(A_onedirec,1);
A = [A_onedirec;-A_onedirec];
A_t = (A-A.^2)/2;
R_half=[0.0192 0.0452 0.057 0.0132 0.0472];
X_half=[0.0575 0.1652 0.1737 0.0379 0.1983];
C_half=diag(1./(R_half + 1i*X_half));
C = [C_half,zeros(length(C_half));zeros(length(C_half)),C_half];
grid.Y=A'*C*A;

grid.K=[0 0 0 0 0]; %stiffness of the primary frequency controller

Gsh = [0 0 0 0 0];
Bsh = [0 0 0 19 0];
Csh_half= diag(Gsh + 1i*Bsh);
grid.Csh = [Csh_half,zeros(length(Csh_half));zeros(length(Csh_half)),Csh_half]

%%state initialisation
state.v = zeros(grid.n,1);
state.theta = zeros(grid.n,1);
state.p_g = zeros(grid.n,1);
state.q_g = zeros(grid.n,1);
state.p_ref = zeros(grid.n,1);
state.i = zeros(grid.m,1);
state.f = 0;

%%equality constraints
    function h = h(state,grid)
        temp1 = diag(u)*conj(grid.Y*u));
        temp3 = (1 + K*f)*state.p_ref;
        temp4 = abs((grid.C*grid.A-grid.Csh*grid.A_t)*u).^2 - i.^2
        
        h=[real(temp1), imag(temp1), temp3, temp4]';
    end

end