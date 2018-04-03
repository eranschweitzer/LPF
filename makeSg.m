function Sg = makeSg(vars, Sp, bidx, N)

Sg = sparse([find(bidx.ref);find(bidx.pv)], 1, ...
    [vars.Pref+1i*vars.Qref; Sp.Pg(bidx.pv)+1i*vars.Qpv], N.t, 1);