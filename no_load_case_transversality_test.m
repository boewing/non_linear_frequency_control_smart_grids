
mygrid = Grid();

mygrid.Y*[1 1 1 1]'
x = State(mygrid);
x.setx([1 1 1 1,zeros(1,4*mygrid.n +4*mygrid.m +1)]');
T = [mygrid.E; Physics.n_h(x, mygrid)];
    ABCD = T(1+1+2*mygrid.n:1+4*mygrid.n,1:2*mygrid.n);
    A = ABCD(1:mygrid.n,1:mygrid.n)
    B = ABCD(1:mygrid.n,1+mygrid.n:2*mygrid.n)
    C = ABCD(1+mygrid.n:2*mygrid.n,1:mygrid.n);
    D = ABCD(1+mygrid.n:2*mygrid.n,1+mygrid.n:2*mygrid.n);
    T_small = [A(:,mygrid.is_PQ) B(:,2:mygrid.n) -mygrid.K; C(mygrid.is_PQ,mygrid.is_PQ) D(mygrid.is_PQ,2:mygrid.n) zeros(sum(mygrid.is_PQ),1)];