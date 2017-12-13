%%building up the network from constants
mygrid = Grid();

%%state initialisation
x = State(mygrid);

%run the simulation and record
iterations=500;
myrec = Recorder(x, Controller.Jt(x,mygrid),norm(Physics.h(x,mygrid)),iterations);

for k=1:iterations
    
    %calculate a feasible gradient
    %d must be in the tangent plane defined by nabla_h(x), nabla_h(x)*d=0
    %d must be as close to the negative gradient of the cost function as possible:
    %analytic version without constraints: d = n_Jt(x,mygrid)' - n_h(x,mygrid)'*inv(n_h(x,mygrid)*n_h(x,mygrid)')*n_h(x,mygrid)*n_Jt(x,mygrid)';
    d = Controller.getStep(x, mygrid, k);
    
    %make a step
    x.setx(x.getx() + d);
    
    %retraction
    assert(abs(x.theta(1)) < 1e-3, 'reference bus has non-zero voltage angle');              %check that the angle reference stays 0
    x = Physics.retraction(x, mygrid);
    %assert(max(abs(Physics.h(x,mygrid))) < 1e-3, 'retraction failed');    %check that the retraction works makes h(x)=0
    
    %check if retraction is transverse to the tangent plane
    T = [mygrid.E; Physics.n_h(x, mygrid)]; %mygrid.E_short;
    if rank(T) < 5*mygrid.n + 2*mygrid.m + 1
        disp(['non-transverse retraction: Rank(T) = ', num2str(rank(T)), ' in iteration ', num2str(k)]);
    end
    
    %recording
    myrec.count();
    myrec.store(x, Controller.Jt(x,mygrid),norm(Physics.h(x,mygrid)));
    myrec.store_limits(x,mygrid);
    
end

myrec.plotAll();