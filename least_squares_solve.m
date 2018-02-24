function x = least_squares_solve(op,F,T,E,Sb,Sp,bidx,theta_ref,N)

np  = length(op);
A   = cell(np, 1);
b   = cell(np,1);
for k = 1:length(op)
    [Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(op{k}.ids,F,T,E,op{k}.u0,Sb,Sp);
    [I,U] = pv_ref_mats(bidx,N,op{k}.u0,op{k}.ids);
    
    A{k} = [Art,  Aru, sparse(N.t,N.pv), U.ref.p,           sparse(N.t,N.ref);
            Amt,  Amu, U.pv,             sparse(N.t,N.ref), U.ref.q];
    b{k} = [br;bm];
    
    if k == 1
       Afix = [sparse(N.pv,N.t),  I.pv, sparse(N.pv,N.pv+2*N.ref);
               I.ref,            sparse(N.ref,N.t+N.pv+2*N.ref);
               sparse(N.ref,N.t), I.ref, sparse(N.ref,N.pv+2*N.ref)]; 
       bfix = [U.bu.pv; theta_ref; U.bu.ref];
    end
end

A = [vertcat(A{:}); Afix];
b = [vertcat(b{:}); bfix];

x = A\b;