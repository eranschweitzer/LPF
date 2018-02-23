function [x, A,b] = single_solve(Art,Aru,Amt,Amu,br,bm,I,U,theta_ref,N)

A = [Art,           Aru, sparse(N.t,N.pv), U.ref.p,         sparse(N.t,N.ref);
     Amt,           Amu, U.pv,         sparse(N.t,N.ref), U.ref.q;
     sparse(N.pv,N.t), I.pv,sparse(N.pv,N.pv+2*N.ref);
     I.ref,         sparse(N.ref,N.t+N.pv+2*N.ref);
     sparse(N.ref,N.t),I.ref,sparse(N.ref,N.pv+2*N.ref)];

if ~isfield(U,'bu')
    bupv = zeros(N.pv,1);
    buref = zeros(N.ref,1);
else
    bupv = U.bu.pv;
    buref= U.bu.ref;
end
b = [br;bm;bupv;theta_ref;buref];

x = A\b;
