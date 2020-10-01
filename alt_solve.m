function varstest = alt_solve(ids,F,T,E,Sb,Sp,bidx,theta_ref,N,u0,varargin)
[Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(ids,F,T,E,u0,Sb,Sp);
[matI,matU] = pv_ref_mats(bidx,N,u0,ids);
A11 = [Art(bidx.pq, bidx.pq), Art(bidx.pq, bidx.pv), Aru(bidx.pq,bidx.pq); 
       Art(bidx.pv, bidx.pq), Art(bidx.pv, bidx.pv), Aru(bidx.pv,bidx.pq);
       Amt(bidx.pq, bidx.pq), Amt(bidx.pq, bidx.pv), Amu(bidx.pq,bidx.pq)];

A12 = sparse(2*N.pq + N.pv, N.pv + 2*N.ref);
A21 = A12';

A22 = sparse(N.pv+2*N.ref,N.pv+2*N.ref);
   
C1  = [Amt(bidx.pv, bidx.pq), Amt(bidx.pv, bidx.pv), Amu(bidx.pv, bidx.pq);
       Art(bidx.ref, bidx.pq), Art(bidx.ref,bidx.pv), Aru(bidx.ref,bidx.pq);
       Amt(bidx.ref, bidx.pq), Amt(bidx.ref,bidx.pv), Amu(bidx.ref,bidx.pq)];

C2  = blkdiag(matU.pv(bidx.pv,:), matU.ref.p(bidx.ref,:), matU.ref.q(bidx.ref,:));

B1  = [Aru(bidx.pq, bidx.pv), Art(bidx.pq, bidx.ref), Aru(bidx.pq, bidx.ref);
       Aru(bidx.pv, bidx.pv), Art(bidx.pv, bidx.ref), Aru(bidx.pv, bidx.ref);
       Amu(bidx.pq, bidx.pv), Amt(bidx.pq, bidx.ref), Amu(bidx.pq, bidx.ref)];
   
B2  = speye(N.pv+2*N.ref);
% B2  = sparse(N.pv+1:N.pv+2*N.ref, N.pv+1:N.pv+2*N.ref, 1, N.pv+2*N.ref, N.pv+2*N.ref);

D   = [Amu(bidx.pv, bidx.pv), Amt(bidx.pv, bidx.ref), Amu(bidx.pv, bidx.ref);
       Aru(bidx.ref, bidx.pv), Art(bidx.ref, bidx.ref), Aru(bidx.ref, bidx.ref);
       Amu(bidx.ref, bidx.pv), Amt(bidx.ref, bidx.ref), Amu(bidx.ref, bidx.ref)];
   

f   = [br(bidx.pq); br(bidx.pv); bm(bidx.pq)]; 
g   = [bm(bidx.pv); br(bidx.ref); bm(bidx.ref)];
% u0ad= [Amu(bidx.pv, bidx.pv), Amu(bidx.pv, bidx.ref);
%        Aru(bidx.ref, bidx.pv), Aru(bidx.ref, bidx.ref);
%        Amu(bidx.ref, bidx.pv), Amu(bidx.ref, bidx.ref)];
% g   = g - u0ad*[u0(bidx.pv); u0(bidx.ref)];
y   = [u0(bidx.pv); theta_ref; u0(bidx.ref)];

%%
x1 = A11 \(f-B1*y);
varstest = struct( 'u', u0, 'theta', theta_ref*ones(N.t,1));

varstest.theta(bidx.pq) = x1(1:N.pq);
varstest.theta(bidx.pv) = x1(N.pq+1:N.pq+N.pv);
varstest.u(bidx.pq)     = x1(N.pq+N.pv+1:2*N.pq+N.pv);

varstest.v = exp(varstest.u + 1i*varstest.theta);

Ybus = myMakeYbus(F,T,Sb);
scalc = diags(varstest.v)*conj(Ybus*varstest.v) ;%+ Sp.Pd(bidx.ref) + 1i*Sp.Qd(bidx.ref);
varstest.Qpv   = imag(scalc(bidx.pv)) + Sp.Qd(bidx.pv); 
varstest.Pref = real(scalc(bidx.ref)) + Sp.Pd(bidx.ref);
varstest.Qref = imag(scalc(bidx.ref)) + Sp.Qd(bidx.ref);

% Sg   = makeSg(varstest,Sp,bidx,N);
% residual = pfresidual(varstest.v,Ybus,real(Sg) - Sp.Pd, imag(Sg) - Sp.Qd);
% %% plots
% vt_comp_plots(varstest, vtrue, E, residual,bidx) 
%%
return
rhs   = [f-B1*y; g+(speye(N.pv+2*N.ref)-B2-D)*y];
xtest = [A11, A12; A21+C1, A22+C2]\rhs;



varstest = struct( 'u', u0, 'theta', theta_ref*ones(N.t,1));

varstest.theta(bidx.pq) = xtest(1:N.pq);
varstest.theta(bidx.pv) = xtest(N.pq+1:N.pq+N.pv);
varstest.u(bidx.pq)     = xtest(N.pq+N.pv+1:2*N.pq+N.pv);
varstest.Qpv            = xtest(2*N.pq+N.pv+1:2*N.pq+2*N.pv);
varstest.Pref           = xtest(2*N.pq+2*N.pv+1:2*N.pq+2*N.pv+N.ref);
varstest.Qref           = xtest(2*N.pq+2*N.pv+N.ref+1:end);


Atest2= [A11, A12, B1; A21, A22, B2; C1, C2, D];
btest2= [f;y;g];
xtest2= Atest2\btest2;

varstest2 = struct( 'u', u0, 'theta', theta_ref*ones(N.t,1));

varstest2.theta(bidx.pq) = xtest2(1:N.pq);
varstest2.theta(bidx.pv) = xtest2(N.pq+1:N.pq+N.pv);
varstest2.u(bidx.pq)     = xtest2(N.pq+N.pv+1:2*N.pq+N.pv);
varstest2.Qpv            = xtest2(2*N.pq+N.pv+1:2*N.pq+2*N.pv);
varstest2.Pref           = xtest2(2*N.pq+2*N.pv+1:2*N.pq+2*N.pv+N.ref);
varstest2.Qref           = xtest2(2*N.pq+2*N.pv+N.ref+1:2*N.pq+2*N.pv+N.ref+2*N.ref);

%% third test

pq = find(bidx.pq);
pv = find(bidx.pv);
ref = find(bidx.ref);
% order is:
%             theta_pq, theta_pv, u_pq,  Qpv, Pref, Qref,                   upv, theta_ref, uref
permutation = [pq;       pv;      N.t+pq; ((2*N.t+1):(2*(N.t+N.ref)+N.pv))'; N.t+pv; ref; N.t+ref];
Atest3     = A(permutation,permutation);
btest3     = b(permutation);

xtest3     = Atest3\btest3;

varstest3 = struct( 'u', u0, 'theta', theta_ref*ones(N.t,1));

varstest3.theta(bidx.pq) = xtest3(1:N.pq);
varstest3.theta(bidx.pv) = xtest3(N.pq+1:N.pq+N.pv);
varstest3.u(bidx.pq)     = xtest3(N.pq+N.pv+1:2*N.pq+N.pv);
varstest3.Qpv            = xtest3(2*N.pq+N.pv+1:2*N.pq+2*N.pv);
varstest3.Pref           = xtest3(2*N.pq+2*N.pv+1:2*N.pq+2*N.pv+N.ref);
varstest3.Qref           = xtest3(2*N.pq+2*N.pv+N.ref+1:2*N.pq+2*N.pv+2*N.ref);


cidx1= 1:size(A11,2);
cidx2= size(A11,2)+1:size(A11,2)+size(A12,2);
cidx3= size(A11,2)+size(A12,2)+1:size(A11,2)+size(A12,2)+size(B1,2);
ridx1 = 1:size(A11,1);
ridx2 = size(A11,1)+1:size(A11,1)+size(A21,1);
ridx3 = size(A11,1)+size(A21,1)+1:size(A11,1)+size(A21,1)+size(C1,1);

Atest3(ridx1,cidx1) - A11
Atest3(ridx1,cidx2) - A12
Atest3(ridx1,cidx3) - B1

Atest3(ridx2,cidx1) - A21
Atest3(ridx2,cidx2) - A22 
Atest3(ridx2,cidx3) - B2  

Atest3(ridx3,cidx1) - C1
Atest3(ridx3,cidx2) - C2
Atest3(ridx3,cidx3) - D

