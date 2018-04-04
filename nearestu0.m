function out = nearestu0(F,T,Sb,N,bidx)
Ytr = myMakeYbus(F,T,Sb,'tronly',true);
%order as pq buses then pv buses
Y   = [Ytr(bidx.pq, bidx.pq), Ytr(bidx.pq, bidx.pv);
       Ytr(bidx.pv, bidx.pq), Ytr(bidx.pv, bidx.pv)];
if N.t < 1e4
    Z   = Y\speye(N.pv + N.pq);
    zdiag = diag(Z);
    Z = abs(ones(N.pv+N.pq,1)*zdiag.' + zdiag*ones(1,N.pv+N.pq) - Z - Z.');
    [zpv,out] = min(Z(1:N.pq,N.pq+1:end),[],2);
    out(abs(zdiag(1:N.pq)) < zpv) = 0;
else
    [l,u] = lu(Y);

    out = zeros(N.pq,1);
    for k = 1:N.pq
        rhs  = sparse( [k*ones(N.pv,1); (N.pq+1:N.pq+N.pv)'], [1:N.pv,1:N.pv]', ...
            [ones(N.pv,1);-ones(N.pv,1)], N.pv+N.pq,N.pv);
        Ztmp = u\(l\rhs);
        Zth  = Ztmp(k,:) - diag(Ztmp(N.pq+1:N.pq+N.pv,1:N.pv)).';
        [zpv,pvidx] = min(abs(Zth));
        Zref = u\(l\sparse(k,1,1,N.pv+N.pq,1));
        if abs(Zref(k)) > zpv
            %%% closest is pv bus number pvidx
            %%% if Zref <= zpv then the closest bus is the reference which will
            %%% be indicated by a value of 0 in the output.
            out(k) = pvidx;
        end
    end
end